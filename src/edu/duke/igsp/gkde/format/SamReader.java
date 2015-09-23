
/*****************************************************************************
  SamReader.java

  (c) 2014 - Jeremy Wang
  Department of Genetics
  University of North Carolina at Chapel Hill
  jeremy@unc.edu
  
  Modified from BedReader.java by Alan Boyle
  aboyle@gmail.com

  Licensed under the GNU General Public License 3.0 license.

  This file is part of F-seq.

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

******************************************************************************/

package edu.duke.igsp.gkde.format;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import edu.duke.igsp.gkde.KDEChromosome;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * Sam/Bam files are not necessarily sorted. Must load everything before sort.
 * Assumes all come from one chromosome
 * 
 * @author jhg9
 * 
 * Modified to keep strand information
 * @author apb9
 * 
 * Modified from BedReader to read SAM/BAM files using HTSJDK
 * @author jrwang
 * 
 */
public class SamReader {

  public static KDEChromosome[] read(File[] files, int weight_clip) throws IOException {

    HashMap<String, ArrayList<KDEChromosome.Sequence>> chrMap = new HashMap<String, ArrayList<KDEChromosome.Sequence>>();

    String currentChr = null;
    ArrayList<KDEChromosome.Sequence> currentCuts = null;
    boolean lengthSet = false;
    int sequenceLength = 0;
    
    for (int i = 0; i < files.length; ++i) {

      //final htsjdk.samtools.SamReader reader = SamReaderFactory.makeDefault().open(files[i]);
      final htsjdk.samtools.SamReaderFactory factory = SamReaderFactory.make();
      factory.validationStringency(ValidationStringency.SILENT); // STRICT, LENIENT, SILENT
      final htsjdk.samtools.SamReader reader = factory.open(files[i]);

      for (final SAMRecord samRecord : reader) {
        String chrom = samRecord.getReferenceName();
        if (chrom != currentChr) {
          if (!chrMap.containsKey(chrom)) {
            chrMap.put(chrom, new ArrayList<KDEChromosome.Sequence>());
          }
          currentChr = chrom;
          currentCuts = chrMap.get(chrom);
        }
        
        KDEChromosome.Sequence seq = new KDEChromosome.Sequence();
        
        try {
          long s = samRecord.getAlignmentStart();
          long e = samRecord.getAlignmentEnd();
          long diff = e - s;
          if(!lengthSet) { // assume all reads are the same length
        	  sequenceLength = (int)diff;
        	  lengthSet = true;
          }
          double weight;
          try {
            weight = samRecord.getDoubleAttribute("XW"); // read weight
            if(weight_clip != 0) {
              if(weight > weight_clip) {
                weight = weight_clip;
              } else if(weight < 1.0f/weight_clip) {
                weight = 1.0f/weight_clip;
              }
            }
          } catch (Exception ex) {
            weight = 1.0f;
          }
          if(samRecord.getReadNegativeStrandFlag()) {
            seq = new KDEChromosome.Sequence(s, false, weight);
          } else {
            seq = new KDEChromosome.Sequence(e, true, weight);
          }
          currentCuts.add(seq);
        } catch (NumberFormatException e) {
          badFile(files[i]);
        }
      }

      reader.close();
    }

    currentCuts = null;

    KDEChromosome[] chrs = new KDEChromosome[chrMap.size()];
    int i = 0;
    String[] chrnames = chrMap.keySet().toArray(new String[0]);
    for (String chr : chrnames) {
      currentCuts = chrMap.remove(chr);
      KDEChromosome.Sequence[] cuts = new KDEChromosome.Sequence[currentCuts.size()];
      for (int j = 0; j < cuts.length; ++j)
        cuts[j] = currentCuts.get(j);

      Arrays.sort(cuts);
      chrs[i++] = new KDEChromosome(chr, cuts, sequenceLength);
    }
    return chrs;
  }

  private static void badFile(File f) throws IOException {
    throw new IOException("Bad '.bed' format for file " + f.getAbsolutePath());
  }
}
