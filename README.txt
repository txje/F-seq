Welcome to F-Seq.

-----------------------------------
This fork is intended to modify F-seq to accept real-valued read weights
(instead of 1) and is a work in progress.

In BED format, the score field will be parsed as a floating-point value for
the read's weight. It need NOT be in the range 0-1000.

In SAM/BAM format, the floating-point score should be in an optional tag named
"XW".

Remaining points of confusion:
Main.java:~320 Is the random cut distribution for computing standard deviation
still accurate? Same at ~350 for wgThreshold

The current version of HTSJDK (samtools/htsjdk commit 9c01b696ab) does not
play nice with the current version of Snappy (1.1.1.4), so I basically
disabled Snappy in HTSJDK by changing
htsjdk/samtools/util/SnappyLoader.java:86 to just ask "if(true)"

Another token of knowledge --
pysam will write BAM files with optional tags with the type 'd' (double),
however this is not allowed in the "official" spec, so HTSJDK will not read a
BAM file with 'd' typed tags. It is trivial to add this functionality to
HTSJDK. Do so. Remember to add a getDoubleAttribute() function to the SamReader.

See txje/sequence-bias-adjustment for more information on read weights
-----------------------------------

This software requires Java version 1.5 or greater.  
To see your version of java (or if it is installed), type 'java -version'

If java is not installed or you do not have the correct version, download at 
http://java.sun.com/javase/downloads/index.jsp

Current accepted input formats:
Bed, http://genome.ucsc.edu/FAQ/FAQformat#format1
Sam/Bam, http://samtools.github.io/hts-specs/SAMv1.pdf

Currently accepted output formats:
Wiggle, http://genome.ucsc.edu/FAQ/FAQformat#format6
Bed (see above)

-----------------------------------
Unix (Linux)

Build with ant from the base directory:
~$ git clone https://github.com/aboyle/F-seq.git
~$ cd F-seq/
~/F-seq$ ant

This will build F-seq and package it in the dist~ folder. To then run F-seq:
~/F-seq$ cd dist~/
~/F-seq/dist~$ tar -xvf fseq.tgz
~/F-seq/dist~$ cd fseq/bin/
~/F-seq/dist~/fseq/bin$ ./fseq
F-Seq Version 1.85
usage: fseq [options]... [file(s)]...
 -b <background dir>     background directory (default=none)
 -c <arg>                genomic count of sequence reads (defualt =
                         calculated)
 -d <input dir>          input directory (default=current directory)
 -f <arg>                fragment size (default=estimated from data)
 -h                      print usage
 -l <arg>                feature length (default=600)
 -o <output dir>         output directory (default=current directory)
 -of <wig | bed | npf>   output format (default wig)
 -p <ploidy dir>         ploidy/input directory (default=none)
 -s <arg>                wiggle track step (default=1)
 -t <arg>                threshold (standard deviations) (default=4.0)
 -v                      verbose output
 -wg <arg>               wg threshold set (defualt = calculated)
 

Make sure 'bin/fseq' is executable (chmod 0755 bin/fseq)

For a list of options, type 'fseq -h'

Example: fseq -v -of wig chr1.bed chr2.bed

This takes as input the chr1.bed and chr2.bed files. 
Will use verbose output and outputs the densities in the wiggle format.

-----------------------------------
Troubleshooting

A likely cause for errors is an "OutOfMemory" exception.  
To increase the available memory to the java virtual machine, edit
'bin/fseq' file and change the JAVA_OPTS property to increase the heap size.

Users also typically complain about an 'ArrayOutOfBounds' exception. This is due to a low number of reads
and is fixed in the latest version of F-seq. You can also get around this error by using the -f option
to set your fragment size.

-----------------------------------
License

Licensed under the GNU General Public License 3.0 license.

F-seq is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

F-seq is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with F-seq.  If not, see <http://www.gnu.org/licenses/>.
