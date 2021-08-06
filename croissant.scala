import scala.io._
import java.io._
import org.rogach.scallop._; // Command-line parsing
import scala.annotation.tailrec
import scala.language.reflectiveCalls // Prevents warnings on compile
import scala.collection.mutable.{Seq => MSeq}
import java.util.zip._

// Version 2

class HaploConf(arguments: Seq[String]) extends ScallopConf(arguments) {
    version("Croissant 0.2 (c) 2015-2016 Gabriel Keeble-Gagnere, DEDJTR")
    val alignment = opt[String](default = None,
        descr="Gydle alignment file")
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

    def convertLineToMismatches(l: Vector[String], sepChar: String): (Seq[(Int, Char)], Seq[(Int, Int)]) = {

        def convertStrand(in: String): String = {
            in match {
                case "1+" | "+" => "+1"
                case "1-" | "-" => "-1"
                case "2+" => "+2"
                case "2-" => "-2"
            }
        }

        def cigarToOffsets(cigar: String): Seq[Int] = {
            val pattern = "([0-9]*)([im])".r
            var currPos = 0
            var mms = Vector(): Vector[Int]
            for (pattern(count,kind) <- pattern findAllIn cigar) {
                currPos += count.toInt
                if (kind == "m") {
                    mms = mms :+ currPos
                }
            }
            mms
        }

        def alignmentStringToMismatches(a: Vector[String], seq: String): (Seq[(Int,Char)], (Int, Int)) = {
            val strand = convertStrand(a(0))(0)
            val pair = convertStrand(a(0))(1)
            val qstart = a(7).toInt // NOTE: query/subject can change with ordering of gydle output
            val qend = a(8).toInt
            val sstart = a(9).toInt
            val send = a(10).toInt
            val cigar = a(11)
            if (strand == '+')
                (cigarToOffsets(cigar).map(x => (x + sstart, seq(x + qstart - 1))), (sstart, send))
            else
                (cigarToOffsets(cigar).map(x => (x + sstart, (reverseComplement(seq.slice(qstart-1,qend)))(x))), (sstart, send))
        }

        val id = l(2)
        val hits = l(7).toString.split("-").map(x => x.toInt).sum
        val alStrings = l.slice(8, 8 + hits).map(x => x.split(" ").toVector)
        val seqs = l(8 + hits).split(sepChar)
        val read1Als = alStrings.filter(x => convertStrand(x(0))(1).toInt == '1')
        val read2Als = alStrings.filter(x => convertStrand(x(0))(1).toInt == '2')
        val mm1 = read1Als.map(x => alignmentStringToMismatches(x, seqs(0)))
        val mm2 = read2Als.map(x => alignmentStringToMismatches(x, seqs(1)))
        ((mm1.map(x => x._1) ++ mm2.map(x => x._1)).flatten.sortBy(_._1), (mm1.map(x => x._2) ++ mm2.map(x => x._2)))

    }

    def main(args: Array[String]) {

        val conf = new HaploConf(args)

        def reader(file: Option[String]): BufferedSource = {
          file match {
            case Some(file) =>
              val stream = (if (file.endsWith("gz"))
                new GZIPInputStream(new FileInputStream(file))
              else
                new FileInputStream(file)
              )
              new BufferedSource(stream)
            case None => Source.stdin
          }
        }

        val s = conf.start()
        val w = conf.window()
        val matePair = conf.mp()
        val covOnly = conf.covOnly()

        val sepChar = if (matePair) "~" else "-"

        val in = reader(conf.alignment.toOption).getLines.buffered

        val firstLine = in.head.split("\t").toVector
        var targetNow = firstLine(0)
        var targetSize = firstLine(1).toInt

        // FIRST PASS

        while (in.hasNext) {

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
                /*
                
                The logic here:

                If the coverage at this base is greater than 1, then if haps=0 (eg: if no mismatches
                here), then number of haplotypes is 1. Otherwise, haps is the right number.

                If coverage at this base is exactly 1, then the number of haplotypes is 1.

                Finally, if none of the above, coverage must be 0, so haps = 0.

                */

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
        }
    }
}

//Croissant.main(args)
