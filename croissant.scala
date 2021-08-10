import scala.io._
import java.io._
import org.rogach.scallop._ // Command-line parsing
import scala.annotation.tailrec
import scala.language.reflectiveCalls // Prevents warnings on compile
import scala.collection.mutable.{Seq => MSeq}
import scala.collection.mutable.Map
import java.util.zip._
import htsjdk.samtools._
import htsjdk.samtools.util.Log
import htsjdk.samtools.util.ProgressLogger
import java.net.InetAddress
import java.util.Arrays
import java.util.List
import java.util.stream.Collectors
import scala.collection.JavaConversions._

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
    val start = opt[Int](default = Some(0),
        descr="Starting base position")
    val end = opt[Int](default = Some(-1),
        descr="End base position")
    verify()
}


case class Mismatch(position: Int, allele: Char)

case class Extent(start: Int, end: Int)

case class Read(mismatches: Seq[Mismatch], location: Extent)

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

    def mdTagToMismatches(md: String, startpos: Int): Seq[Mismatch] = {
      val pattern = "([0-9]*)([ACGT])".r
      var currPos = 0
      var mms = Vector(): Vector[Mismatch]
      for (pattern(count,kind) <- pattern findAllIn md) {
        currPos += count.toInt
        mms = mms :+ Mismatch(currPos+startpos,kind(0))
      }
      mms
    }

    def convertSamRecordToMismatches(record: SAMRecord): Read = {
      Read(mdTagToMismatches(record.getAttribute("MD").toString, record.getAlignmentStart()),
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
        val s = conf.start()
        val w = conf.window()
        val matePair = conf.mp()
        val covOnly = conf.covOnly()

        // Open BAM file
        val bamReader = readerFactory.open(inputFile);
        val header = bamReader.getFileHeader();
        // Check if sorted by coordinate (ie: by reference)
        val sorted = header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate);

        var readCount = 0

        if (!sorted) {
          throw new Exception("SAM/BAM needs to be sorted (eg: with samtools sort)")
        }
        else {
          var targetNow = ""
          var refCount = 0
          for {record <- bamReader; if (!(record.getReadUnmappedFlag()))} {
            if (!(record.getCigarString() contains "D")) { // Ignore indels for now
              if (record.getReferenceName() != targetNow) {
                targetNow = record.getReferenceName()
                println(s"Target now ${targetNow}")
                // New reference set-up
              }
              val read = convertSamRecordToMismatches(record)
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

        //val firstLine = in.head.split("\t").toVector
        //var targetNow = firstLine(0)
        //var targetSize = firstLine(1).toInt

        // FIRST PASS

        /*while (in.hasNext) {

          var base = MSeq.fill(targetSize + 1)((0, Vector(): Vector[(Seq[(Int, Char)], Seq[(Int, Int)])]))
          var lines = Vector(): Vector[(Seq[(Int, Char)], Seq[(Int, Int)])]

          while (in.hasNext && in.head.split("\t")(0) == targetNow) {
            val line = in.next().split("\t").toVector
            val read = convertLineToMismatches(line, sepChar) // Convert line to "read".

            // Add this "read" to all positions for which it has a mismatch.

            for (mm <- read._1.map(_._1)) {
              base(mm) = (base(mm)._1, base(mm)._2 :+ read)
            }

            lines = lines :+ read
          }

          // After first pass, every base with a mismatch has all reads with a mismatch at that base in
          // its stack.

          // SECOND PASS

          for (line <- lines) {
              for (region <- line._2) {
                  for (i <- region._1 to region._2) {
                      base(i) = (base(i)._1 + 1, base(i)._2) // Increase coverage for bases of read.

                      // If position i is a mismatch and this read (which covers this base) has a mismatch (not
                      // necessarily at this base), add read to stack.

                      if ((line._1.size > 0) && (base(i)._2.size > 0))
                          base(i) = (base(i)._1, base(i)._2 :+ line)
                  }
              }
          }

          // After second pass, the base array has correct coverage count (counting all reads) and every
          // position with a mismatch includes every read with a mismatch (not necessarily at that
          // position) covering that base.

          if (covOnly) {
            for (i <- 0 to base.size - 1) {
              println(s"${targetNow}\t${i+1}\t${base(i)._1}")
            }
          }
          else {
            var haps = 1 // Stores the last seen haplotype number.
            var until = base.size
            for (i <- 0 to base.size - 1) { // Step through each base.
                val haplotypes = if (base(i)._2.size > 0) {

                    // Position has a mismatch - work out haplotype.

                    // Calculate positions of all mismatches from all reads stored at this base.

                    val mms = base(i)._2.map(_._1).flatten.map(_._1).toSet.toVector.sorted

                    val index = mms.indexOf(i)
                    var myset = Set[String]()
                    var toCheck = Vector[Int]()
                    val min = if ((index - w) >= 0) index - w else 0
                    val max = if ((index + w) < mms.size) index + w else mms.size - 1

                    for (r <- min to max) toCheck = toCheck :+ mms(r)

                    until = toCheck.last // haps is valid until this base.

                    for (r <- base(i)._2) { // Step through each read at this base.
                        var mystring = ""
                        for (m <- toCheck) {
                            var here = r._1.filter(_._1 == m)
                            if (here.size > 0) mystring = mystring + here.head._2
                            else if (r._2.filter(x => x._1 <= m && x._2 >= m).size > 0) mystring = mystring + "*"
                            else mystring = mystring + " "
                        }
                        myset = myset + mystring
                    }
                    val hapstrings = myset.map(x => x.trim).filter(_.size == toCheck.size)
                    haps = hapstrings.size
                    if (!(hapstrings.contains("*"*toCheck.size)) && (base(i)._2.map(_._1).flatten.filter(_._1 == i).size < base(i)._1))
                        haps += 1 // Add the reference haplotype.
                    haps
                }

                // The logic here:
                // If the coverage at this base is greater than 1, then if haps=0 (eg: if no mismatches
                // here), then number of haplotypes is 1. Otherwise, haps is the right number.
                // If coverage at this base is exactly 1, then the number of haplotypes is 1.
                // Finally, if none of the above, coverage must be 0, so haps = 0.

                else if (base(i)._1 > 1) { // If there are no reads above this base with mismatches, but
                                           // coverage > 0.
                    if ((haps == 0) || (i > until)) { haps = 1 }
                    haps // This is equal to last seen haplotype if there are no reads with mismatches above
                         // this base.
                }
                else if (base(i)._1 == 1) 1
                else { haps = 0; 0 }
                println(s"$targetNow\t$i\t${base(i)._1}\t$haplotypes")
            }
          }
          if (in.hasNext) {
            targetNow = in.head.split("\t")(0)
            targetSize = in.head.split("\t")(1).toInt
          }
        }*/
    }
}

//Croissant.main(args)
