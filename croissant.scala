import scala.io._
import java.io._
import org.rogach.scallop._ // Command-line parsing
import scala.annotation.tailrec
import scala.language.reflectiveCalls // Prevents warnings on compile
import scala.collection.mutable.{Seq => MSeq}
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Map
import java.util.zip._
import htsjdk.samtools._
import htsjdk.samtools.util.Log
import htsjdk.samtools.util.ProgressLogger
import java.net.InetAddress
import java.util.Arrays
import java.util.List
import java.util.stream.Collectors
import scala.jdk.CollectionConverters._

// Version 2

class HaploConf(arguments: Seq[String]) extends ScallopConf(arguments) {
    version("Croissant 0.3 (c) 2015-2021 Gabriel Keeble-Gagnere, Agriculture Victoria")
    val alignment = opt[String](default = None,
        descr="Sorted BAM alignment file")
    val verbose = opt[Boolean]()
    val mp = opt[Boolean](default = Some(false),
        descr="Aligned data is mate pair")
    val covOnly = opt[Boolean](default = Some(false),
        descr="Only report coverage (no haplotype information)")
    val window = opt[Int](default = Some(1),
        descr="Number of mismatches around target to consider in haplotype calculation")
    verify()
}


case class Mismatch(position: Int, allele: Char)

case class Extent(start: Int, end: Int)

case class Read(mismatches: Seq[Mismatch], location: Extent)

class AlignmentWindow(window: Int) {
  println("AlignmentWindow created")

  // state tracks for each pos: coverage, extent of reads, and list of reads
  type ReadBuf = ArrayBuffer[Read]
  val state = Map[Int,(Int,(Int,Int),ReadBuf)]()
  var reads = ArrayBuffer[Read]()
  var currentPosition = 1
  var haps = 0

  def addRead(read: Read) {
    for (pos <- currentPosition to read.location.end) {
      if (!(state contains pos)) {
        if (pos >= read.location.start) {
          state += ((pos,(1,(read.location.start,read.location.end),new ReadBuf)))
        }
        else {
          state += ((pos,(0,(pos,pos),new ReadBuf)))
        } 
      }
      else {
        val minEx = if (state(pos)._2._1 > read.location.start) read.location.start else state(pos)._2._1
        val maxEx = if (state(pos)._2._2 < read.location.end) read.location.end else state(pos)._2._2
        state(pos) = (state(pos)._1 + 1, (minEx,maxEx), state(pos)._3)
      }
    }
    // Only track reads with mismatches
    if (read.mismatches.size > 0) {
      reads :+= read
    }
    for (m <- read.mismatches) {
      for (r <- readsOverlappingPosition(m.position)) {
        state(m.position)._3 += r
      }
    }
    currentPosition = read.location.end
    callAndDropBefore(currentPosition - window)
  }

  def readsOverlappingPosition(pos: Int): ReadBuf = {
    reads.filter(x => x.location.start <= pos && x.location.end >= pos)    
  }

  def callAndDropBefore(position: Int) {
    val toDrop = state.keys.filter(_ < position).toSeq.sorted
    for (p <- toDrop) {
      println(s"${p}\t${callPosition(p)}") // We have maximum info at this point
      removePosition(p)
    }
    reads = reads.dropWhile(_.location.end < position)
  }

  def callPosition(position: Int): Int = {
    println(s"callPosition: position = ${position}, state(position)._3.size = ${state(position)._3.size}")
    if (state(position)._1 > 0 && state(position)._3.size == 0) {
      haps = 1
    }
    else {
      var myset = Set[String]()
      val haplotypes = if (state(position)._3.size > 0) {

        // Position has a mismatch - work out haplotype.
        // Calculate positions of all mismatches from all reads stored at this base.

        val mms = state(position)._3.map(_.mismatches).flatten.map(_.position).toSet.toVector.sorted
        println("mms",mms)

        for (r <- state(position)._3) { // Step through each read at this base.
          var mystring = ""
          for (m <- mms) {
            var here = r.mismatches.filter(_.position == m)
            // If this read has a mismatch here, add its allele to mystring
            if (here.size > 0) mystring = mystring + here.head.allele
            // If the read covers this position, add "*" to mystring
            else if (r.location.start <= m && r.location.end >= m) mystring = mystring + "*"
            // Otherwise, add " " to mystring
            else mystring = mystring + " "
          }
          println(mystring, r)
          myset = myset + mystring
        }
        println("myset",myset)
        //val hapstrings = myset.map(x => x.trim).filter(_.size == mms.size)
        haps = myset.size
        haps
      }
    }
    haps
  }

  def removePosition(position: Int) {
    state -= position
  }
}

object Croissant {
    /*
        DNA sequence helper methods
        Should they be moved to some outside package eventually?
    */

    def complement(base: Char): Char = {
        base match {
            case 'A' => 'T'
            case 'C' => 'G'
            case 'G' => 'C'
            case 'T' => 'A'
            case 'N' => 'N'
        }
    }

    def reverseComplement(seq: String): String = {
        seq.reverse.map(complement(_))
    }

    def mdTagToMismatches(md: String, startpos: Int, readSeq: String): Seq[Mismatch] = {
      val pattern = "([0-9]*)([ACGT])".r
      var currPos = 0
      var mms = Vector(): Vector[Mismatch]
      for (pattern(count,kind) <- pattern findAllIn md) {
        currPos += count.toInt
        mms = mms :+ Mismatch(currPos+startpos,readSeq(currPos))
        currPos += 1
      }
      mms
    }

    def convertSamRecordToMismatches(record: SAMRecord): Read = {
      Read(mdTagToMismatches(record.getAttribute("MD").toString,
        record.getAlignmentStart(),record.getReadString()),
        Extent(record.getAlignmentStart(), record.getAlignmentEnd()))
    }

    def main(args: Array[String]) {

        val conf = new HaploConf(args)

        val inputFile = new File(conf.alignment.toOption.get)
        val eagerDecode = true //useful to test (realistic) scenarios in which every record is always fully decoded.

        val readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
        if (eagerDecode) {
            readerFactory.enable(SamReaderFactory.Option.EAGERLY_DECODE)
        }

        // Set up dictionary for reference sequence length lookup
        val referenceSizeDict = Map[String,Int]()
        val referenceDict = readerFactory.getFileHeader(inputFile).getSequenceDictionary()
        for (i <- 0 to referenceDict.size() - 1) {
          var seqRec = referenceDict.getSequence(i)
          referenceSizeDict += ((seqRec.getSequenceName(),seqRec.getSequenceLength()))
        }

        // Set arguments
        val w = conf.window()
        val matePair = conf.mp()
        val covOnly = conf.covOnly()

        // Open BAM file
        val bamReader = readerFactory.open(inputFile);
        val header = bamReader.getFileHeader();
        // Check if sorted by coordinate (ie: by reference)
        val sorted = header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate);

        var readCount = 0

        val alnWindow = new AlignmentWindow(500)

        if (!sorted) {
          throw new Exception("SAM/BAM needs to be sorted (eg: with samtools sort)")
        }
        else {
          var targetNow = ""
          var refCount = 0
          for {record <- bamReader.asScala; if (!(record.getReadUnmappedFlag()))} {
            if (!(record.getCigarString() contains "D")) { // Ignore indels for now
              if (record.getReferenceName() != targetNow) {
                targetNow = record.getReferenceName()
                println(s"Target now ${targetNow}")
                // New reference set-up
              }
              val read = convertSamRecordToMismatches(record)
              alnWindow.addRead(read)
              readCount += 1
              if (readCount % 100 == 0) {
                println(s"Processed ${readCount} reads")
              }
                //println((record.getReferenceName(), record.getCigar(), record.getAlignmentStart(),
                //  record.getAlignmentEnd(),"->",record.getAttribute("MD").toString));
                //println(convertSamRecordToMismatches(record)) 
            }
          }
        }

    }
}
