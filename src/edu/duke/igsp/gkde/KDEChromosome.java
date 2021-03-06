/*****************************************************************************
  KDEChromosome.java

  (c) 2008-2013 - Alan Boyle
  Department of Genetics
  Stanford University
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

package edu.duke.igsp.gkde;

import java.util.Arrays;
import java.util.Random;
import java.io.File;

import com.sun.corba.se.impl.logging.UtilSystemException;

import edu.duke.igsp.gkde.background.BffReader;
import edu.duke.igsp.gkde.background.IffReader;
import edu.duke.igsp.gkde.background.WigChromosome;
import edu.duke.igsp.gkde.format.DensityWriter;


public class KDEChromosome {
  
  static int BATCH_SIZE = 1024 * 10;
  
  private long _firstCut;
  private long _lastCut;
  private Sequence[] _cuts;
  private String _chromosome;
  private float _threshold;
  private int _sequenceLength;
  private double _totalWeight;
  
  public KDEChromosome(String chromosome, Sequence[] cuts, int sequenceLength){
    _chromosome = chromosome;
    _cuts = cuts;
    _firstCut = cuts[0].getPosition();
    _lastCut = cuts[cuts.length-1].getPosition();
    _sequenceLength = sequenceLength;

    _totalWeight = 0.0d;
    for(int i = 0; i < _cuts.length; i++) {
      _totalWeight = _totalWeight + _cuts[i].getWeight();
    }
  }
  
  public int getSequenceLength() {
	  return _sequenceLength;
  }
  
  public Sequence[] getCuts() {
	  	return _cuts;
  }
  
  public String getChromosome(){
    return _chromosome;
  }
  
  public long getFirstPos(){
    return _firstCut;
  }
  
  public long getLastPos(){
    return _lastCut;
  }
  
  public long getLength(){
    return _cuts.length;
  }

  public double getTotalWeight(){
	    return _totalWeight;
  }
  
  public void run(Settings settings, DensityWriter dw, boolean verboseFlag, float wg_threshold) throws Exception {
    
//    _threshold = computeThreshold(settings);
    _threshold = wg_threshold;
    dw.setThreshold(_threshold);
    
//    if(verboseFlag){
//    	System.out.println(_chromosome + ": Threshold: " + _threshold);
//    }
    
    Sequence[] cuts = _cuts;
//    long[] cuts = _cuts;
    
    float[] density = new float[BATCH_SIZE];
    
    if(verboseFlag){
      System.out.println(_chromosome + ": first=" + _firstCut + ", last=" + _lastCut);
      for(int i = 0; i < 20; ++i){
        System.out.print(".");
      }
      System.out.println();
    }
    int numBases = (int)Math.abs(_lastCut - _firstCut);
    
    int incr =  numBases / 20;
    int mod = BATCH_SIZE -1;
    int peaks = 0;
    
    long start = System.currentTimeMillis();
    int cutIdx = 0;
    long currentChromPos = 0;
    int arrPos = 0;
    boolean aboveThreshold = false;
    for(int i = 0; i < numBases; ++i){
      currentChromPos = i + _firstCut;
      arrPos = i % BATCH_SIZE;
      density[arrPos] = (float)density(settings, currentChromPos, cutIdx, cuts);
      
      if(!aboveThreshold && density[arrPos] > _threshold){
        aboveThreshold = true;
        ++peaks;
      }else if(aboveThreshold && density[arrPos] < _threshold){
        aboveThreshold = false;
      }
      while(cutIdx < cuts.length && currentChromPos > cuts[cutIdx].getPosition())
        ++cutIdx;
      
      if(verboseFlag && i % incr == 0)
        System.out.print(".");
      if(i % BATCH_SIZE == mod){
        dw.writeDensity(density, 0, BATCH_SIZE); 
      }
    }
    
    int len = numBases % BATCH_SIZE;
    dw.writeDensity(density, 0, len);
    
    if(verboseFlag){
      System.out.println();
      System.out.println(_chromosome + ": Completed in " + (System.currentTimeMillis() - start)/1000d + " seconds.");
      System.out.println(_chromosome + ": Found " + peaks + " peaks.");
      System.out.println("-----------------------");
    }
  }
  
  /*
   * This function is the same as above but reads in a background uniqueness file and adjusts the
   * density estimate based on this.
   */
  public void runBG(Settings settings, DensityWriter dw, boolean verboseFlag, float wg_threshold, File bgfile) throws Exception {
//    _threshold = computeThreshold(settings);
	    _threshold = wg_threshold;
	    dw.setThreshold(_threshold);
	    
		WigChromosome bgchr = null;
		try {
			bgchr = BffReader.read(bgfile);
		} catch (Exception e){
			e.printStackTrace();
		}
//		System.out.println(chr.getChromosome());
//		System.out.println(chr.getLength());
//		System.out.println(chr.getStart());
//	    if(verboseFlag){
//	    	System.out.println(_chromosome + ": Threshold: " + _threshold);
//	    }
	    
	    Sequence[] cuts = _cuts;
//	    long[] cuts = _cuts;
	    
	    float[] density = new float[BATCH_SIZE];
	    
	    if(verboseFlag){
	      System.out.println(_chromosome + ": first=" + _firstCut + ", last=" + _lastCut);
	      for(int i = 0; i < 20; ++i){
	        System.out.print(".");
	      }
	      System.out.println();
	    }
	    int numBases = (int)Math.abs(_lastCut - _firstCut);
	    
	    int incr =  numBases / 20;
	    int mod = BATCH_SIZE -1;
	    int peaks = 0;
	    
	    long start = System.currentTimeMillis();
	    int cutIdx = 0;
	    long currentChromPos = 0;
	    int arrPos = 0;
	    boolean aboveThreshold = false;
	    for(int i = 0; i < numBases; ++i){
	      currentChromPos = i + _firstCut;
	      arrPos = i % BATCH_SIZE;
	      density[arrPos] = (float)bgdensity(settings, currentChromPos, cutIdx, cuts, bgchr);
	      
	      if(!aboveThreshold && density[arrPos] > _threshold){
	        aboveThreshold = true;
	        ++peaks;
	      }else if(aboveThreshold && density[arrPos] < _threshold){
	        aboveThreshold = false;
	      }
	      while(cutIdx < cuts.length && currentChromPos > cuts[cutIdx].getPosition())
	        ++cutIdx;
	      
	      if(verboseFlag && i % incr == 0)
	        System.out.print(".");
	      if(i % BATCH_SIZE == mod){
	        dw.writeDensity(density, 0, BATCH_SIZE); 
	      }
	    }
	    
	    int len = numBases % BATCH_SIZE;
	    dw.writeDensity(density, 0, len);
	    
	    if(verboseFlag){
	      System.out.println();
	      System.out.println(_chromosome + ": Completed in " + (System.currentTimeMillis() - start)/1000d + " seconds.");
	      System.out.println(_chromosome + ": Found " + peaks + " peaks.");
	      System.out.println("-----------------------");
	    }
  }

  /*
   * This function is the same as above but reads in a background uniqueness file and adjusts the
   * density estimate based on this.
   */
  public void run(Settings settings, DensityWriter dw, boolean verboseFlag, float wg_threshold, File[] bgfile, File[] ipfile) throws Exception {
//    _threshold = computeThreshold(settings);
	    _threshold = wg_threshold;
	    dw.setThreshold(_threshold);
	    
	    boolean bg_hit = false;
	    boolean ip_hit = false;
	    boolean bg_used = true;
	    boolean ip_used = true;
	    File bg_file_used = null;
	    File ip_file_used = null;
	    
	    if(bgfile.length == 0) {
	    	// No background used
	    	bg_used = false;
	    } else {
	    	for(int j = 0; j < bgfile.length; ++j) {
	    		if(bgfile[j].getName().equals(_chromosome + ".bff")) {
	    			bg_file_used = bgfile[j];
	    			//System.out.println("Running background on Chromosome " + _chromosome);
	    			bg_hit = true;
	    		}
	    	}
	    	if(!bg_hit) {
	    		System.out.println("No background for Chromosome " + _chromosome);
	    		return;
	    	}
	    }
	    
	    if(ipfile.length == 0) {
	    	// No ploidy used
	    	ip_used = false;
	    } else {
	    	for(int j = 0; j < ipfile.length; ++j) {
	    		if(ipfile[j].getName().equals(_chromosome + ".iff")) {
	    			ip_file_used = ipfile[j];
	    			//System.out.println("Running input on Chromosome " + _chromosome);
	    			ip_hit = true;
	    		}
	    	}
	    	if(!ip_hit) {
	    		System.out.println("No input for Chromosome " + _chromosome);
	    		return;
	    	}
	    }
	    
	    //Load files
    	WigChromosome bgchr = null;
    	WigChromosome ipchr = null;
	    if(bg_used) {
	    	try {
	    		bgchr = BffReader.read(bg_file_used);
	    	} catch (Exception e){
	    		e.printStackTrace();
	    	}
	    }
	    if(ip_used) {
	    	try {
	    		ipchr = IffReader.read(ip_file_used);
	    	} catch (Exception e){
	    		e.printStackTrace();
	    	}
	    }
//		System.out.println(chr.getChromosome());
//		System.out.println(chr.getLength());
//		System.out.println(chr.getStart());
//	    if(verboseFlag){
//	    	System.out.println(_chromosome + ": Threshold: " + _threshold);
//	    }
	    
	    Sequence[] cuts = _cuts;
//	    long[] cuts = _cuts;
	    
	    float[] density = new float[BATCH_SIZE];
	    
	    if(verboseFlag){
	      System.out.println(_chromosome + ": first=" + _firstCut + ", last=" + _lastCut);
	      for(int i = 0; i < 20; ++i){
	        System.out.print(".");
	      }
	      System.out.println();
	    }
	    int numBases = (int)Math.abs(_lastCut - _firstCut);
	    
	    int incr =  numBases / 20;
	    int mod = BATCH_SIZE -1;
	    int peaks = 0;
	    
	    long start = System.currentTimeMillis();
	    int cutIdx = 0;
	    long currentChromPos = 0;
	    int arrPos = 0;
	    boolean aboveThreshold = false;
	    for(int i = 0; i < numBases; ++i){
	      currentChromPos = i + _firstCut;
	      arrPos = i % BATCH_SIZE;
	      
	      if(!bg_used && !ip_used) {
	    	  density[arrPos] = (float)density(settings, currentChromPos, cutIdx, cuts);
	      } else {
	    	  if(bg_used && !ip_used) {
	    		  density[arrPos] = (float)bgdensity(settings, currentChromPos, cutIdx, cuts, bgchr);
	    	  } else {
	    		  if(!bg_used && ip_used) {
	    			  density[arrPos] = (float)ipdensity(settings, currentChromPos, cutIdx, cuts, ipchr);
	    		  } else {
	    			  density[arrPos] = (float)fulldensity(settings, currentChromPos, cutIdx, cuts, bgchr, ipchr);
	    		  }
	    	  }
	      }
//	    	  density[arrPos] = (float)bgdensity(settings, currentChromPos, cutIdx, cuts, bgchr);
	      
	      if(!aboveThreshold && density[arrPos] > _threshold){
	        aboveThreshold = true;
	        ++peaks;
	      }else if(aboveThreshold && density[arrPos] < _threshold){
	        aboveThreshold = false;
	      }
	      while(cutIdx < cuts.length && currentChromPos > cuts[cutIdx].getPosition())
	        ++cutIdx;
	      
	      if(verboseFlag && i % incr == 0)
	        System.out.print(".");
	      if(i % BATCH_SIZE == mod){
	        dw.writeDensity(density, 0, BATCH_SIZE); 
	      }
	    }
	    
	    int len = numBases % BATCH_SIZE;
	    dw.writeDensity(density, 0, len);
	    
	    if(verboseFlag){
	      System.out.println();
	      System.out.println(_chromosome + ": Completed in " + (System.currentTimeMillis() - start)/1000d + " seconds.");
	      System.out.println(_chromosome + ": Found " + peaks + " peaks.");
	      System.out.println("-----------------------");
	    }
  }
  /*  
  private float computeThreshold(Settings settings){
    Random r = new Random();
    
    double size = (int)Math.abs(_lastCut - _firstCut);
    double ncuts = getTotalWeight();
    int totalWindow = 1 + (int)(settings.window * 2);
    
    int cutDensity = (int)((ncuts / size) * totalWindow);
    int thresholdIterations = 1000;
    long[] cuts = new long[cutDensity];
    double[] densities = new double[thresholdIterations];
    
    
    for(int i = 0; i < thresholdIterations; ++i){
      for(int j = 0; j < cuts.length; ++j)
        cuts[j] = r.nextInt(totalWindow);
      Arrays.sort(cuts);
      //densities[i] = density(settings, (long)settings.window, (int)(cuts.length/2.0), cuts);
    }
    
    double mean = Util.mean(densities);
    double std = Util.std(densities);
    
    return (float)(mean + settings.threshold * std);
  }
*/
  
  private static float density(Settings settings, long chromPos, int cutIdx, Sequence[] cuts){
    
    long minPos = chromPos - settings.window;
    long maxPos = chromPos + settings.window;
    
    double[] PRECOMPUTE = settings.precompute;
    
    double sum = 0.0;
    for(int i = cutIdx-1; i > -1; --i){
      if (cuts[i].getPosition() < minPos) break;
      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
      if(settings.experimentType.equals("chip")) {
    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
    		  d = Math.abs((int)(cuts[i].getPosition() + settings.offset - chromPos));
    	  } else if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
          d = Math.abs((int)(cuts[i].getPosition() - settings.offset - chromPos));
    	  } else {
          continue;
        }
      }
    	sum += settings.precompute[d] * cuts[i].getWeight();
    }
    
    for(int i = cutIdx; i < cuts.length; ++i){
      if (cuts[i].getPosition() > maxPos) break;
      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
      if(d > PRECOMPUTE.length-1)
        throw new IllegalStateException();
      if(settings.experimentType.equals("chip")) {
    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
    		  d = Math.abs((int)(cuts[i].getPosition() + settings.offset - chromPos));
    	  } else if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
          d = Math.abs((int)(cuts[i].getPosition() - settings.offset - chromPos));
        } else {
          continue;
        }
      }
    	sum += settings.precompute[d] * cuts[i].getWeight();
    }
    
    return (float)(sum / (double)settings.bandwidth);
  }

  private float bgdensity(Settings settings, long chromPos, int cutIdx, Sequence[] cuts, WigChromosome bgdata){
    
    long minPos = chromPos - settings.window;
    long maxPos = chromPos + settings.window;
    
    double[] PRECOMPUTE = settings.precompute;
   
    int bgOffset = (int)_firstCut - bgdata.getStart(); //convert to 0 based
    
    double sum = 0.0;
    int b = 0;
    
    for(int i = cutIdx-1; i > -1; --i){
      if (cuts[i].getPosition() < minPos) 
        break;
      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
      
      if(settings.experimentType.equals("chip")) {
    		d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
    	  } else if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
          b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
    	  } else {
          continue;
        }
      } else {
    	  if(cuts[i].getStrand()) {
    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
    	  } else {
    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
    	  }
      }
		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
			  sum += settings.precompute[d] * (double)bgdata.getValues()[b] * cuts[i].getWeight();
		  }
    }
    
    for(int i = cutIdx; i < cuts.length; ++i){
      if (cuts[i].getPosition() > maxPos) break;
      
      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
      
      if(d > PRECOMPUTE.length-1)
        throw new IllegalStateException();
      
      if(settings.experimentType.equals("chip")) {
    		d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
    		  b = (int)cuts[i].getPosition() -bgdata.getStart();
    	  } else if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
          b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
        } else {
          continue;
        }
      } else {
    	  if(cuts[i].getStrand()) {
    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
    	  } else {
    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
    	  }
      }
		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
			  sum += settings.precompute[d] * (double)bgdata.getValues()[b] * cuts[i].getWeight();
		  }
    }
    
    return (float)(sum / (double)settings.bandwidth);
  }
  
  private float ipdensity(Settings settings, long chromPos, int cutIdx, Sequence[] cuts, WigChromosome bgdata){
	    
	    long minPos = chromPos - settings.window;
	    long maxPos = chromPos + settings.window;
	    
	    double[] PRECOMPUTE = settings.precompute;
	   
	    int bgOffset = (int)_firstCut - bgdata.getStart(); //convert to 0 based
	    
	    double sum = 0.0;
	    int b = 0;
	    
	    for(int i = cutIdx-1; i > -1; --i){
	      if (cuts[i].getPosition() < minPos) 
	        break;
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      
	      
	      if(settings.experimentType.equals("chip")) {
	    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d) * cuts[i].getWeight();
	    		  }
	    	  } else {
	    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
	    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
	        		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	        			  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d) * cuts[i].getWeight();
	        		  }  		  
	    		  }
	    	  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    	  }
			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
				  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d) * cuts[i].getWeight();
			  }
	      }
	    }
	    
	    for(int i = cutIdx; i < cuts.length; ++i){
	      if (cuts[i].getPosition() > maxPos) break;
	      
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      
	      if(d > PRECOMPUTE.length-1)
	        throw new IllegalStateException();
	      
	      if(settings.experimentType.equals("chip")) {
	    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
	    		  b = (int)cuts[i].getPosition() -bgdata.getStart();
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d) * cuts[i].getWeight();
	    		  }
	    	  } else {
	    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
	    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
	    			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	        			  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d) * cuts[i].getWeight();
	        		  }		  
	    		  }
	    	  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    	  }
			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
				  sum += settings.precompute[d] / ((double)bgdata.getValues()[b]/1000d) * cuts[i].getWeight();
			  }
	      }
	    }
	    
	    return (float)(sum / (double)settings.bandwidth);
	  }
  
  private float fulldensity(Settings settings, long chromPos, int cutIdx, Sequence[] cuts, WigChromosome bgdata, WigChromosome ipdata){
	    
	    long minPos = chromPos - settings.window;
	    long maxPos = chromPos + settings.window;
	    
	    double[] PRECOMPUTE = settings.precompute;
	   
	    int bgOffset = (int)_firstCut - bgdata.getStart(); //convert to 0 based
	    
	    double sum = 0.0;
	    int b;
	    int c;
	    
	    for(int i = cutIdx-1; i > -1; --i){
	      if (cuts[i].getPosition() < minPos) 
	        break;
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      
	      
	      if(settings.experimentType.equals("chip")) {
	    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart(); //index of ip for particular sequence i
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
	    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d) * cuts[i].getWeight();
	    			  }
	    		  }
	    	  } else {
	    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
	    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
		    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
	        		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
		    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
		    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d) * cuts[i].getWeight();
		    			  }
	        		  }  		  
	    		  }
	    	  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart();
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
	    	  }
    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d) * cuts[i].getWeight();
    			  }
    		  }  		
	    	  //sum += settings.precompute[d] * cuts[i].getWeight();
	      }
	    }
	    
	    for(int i = cutIdx; i < cuts.length; ++i){
	      if (cuts[i].getPosition() > maxPos) break;
	      
	      int d = Math.abs((int)(cuts[i].getPosition() - chromPos));
	      
	      if(d > PRECOMPUTE.length-1)
	        throw new IllegalStateException();
	      
	      if(settings.experimentType.equals("chip")) {
	    	  if(cuts[i].getStrand() && cuts[i].getPosition() <= chromPos) {
	    		  d = Math.abs((int)(cuts[i].getPosition() + (int)settings.offset - chromPos));
	    		  b = (int)cuts[i].getPosition() -bgdata.getStart();
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart();
	    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
	    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
	    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000) * cuts[i].getWeight();
	    			  }	    		  }
	    	  } else {
	    		  if(!cuts[i].getStrand() && cuts[i].getPosition() >= chromPos) {
	    			  d = Math.abs((int)(cuts[i].getPosition() - (int)settings.offset - chromPos));
	    			  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength;
		    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
	    			  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
		    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
		    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d) * cuts[i].getWeight();
		    			  }	        		  }		  
	    		  }
	    	  }
	      } else {
	    	  if(cuts[i].getStrand()) {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart(); //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart();
	    	  } else {
	    		  b = (int)cuts[i].getPosition() - bgdata.getStart() - _sequenceLength; //index of bg for particular sequence i
	    		  c = (int)cuts[i].getPosition() - ipdata.getStart() - _sequenceLength;
	    	  }
    		  if(b >= 0 && b < (int)bgdata.getLength() && bgdata.getValues()[b] > 0) {
    			  if(c >= 0 && c < (int)ipdata.getLength() && ipdata.getValues()[c] > 0) {
    				  sum += settings.precompute[d] * (double)bgdata.getValues()[b] / ((double)ipdata.getValues()[c]/1000d) * cuts[i].getWeight();
    			  }
    		  }  		
	    	  //sum += settings.precompute[d] * cuts[i].getWeight();
	      }
	    }
	    
	    return (float)(sum / (double)settings.bandwidth);
	  }
  
  public static class Sequence implements Comparable
  {
      long position;
      boolean strand;
      double weight;
      
      public long getPosition() {
    	  	return position;
      }
      
      public boolean getStrand() {
    	  	return strand;
      }

      public double getWeight() {
          return weight;
      }
      
      public void setPosition(long pos) {
    	  this.position = pos;
      }
    
      public void setStrand(boolean str) {
    	  this.strand = str;
      }

      public void setWeight(double wt) {
    	  this.weight = wt;
      }
      
      public Sequence(long x, boolean y, double w) {
          this.position = x;
          this.strand = y;
          this.weight = w;
      }
      
      public Sequence() {
    	  this.position = 0;
    	  this.strand = false;
        this.weight = 1.0f;
      }
      
      //CompareTo overload
      public int compareTo(Object obj) {
    	  Sequence tmp = (Sequence)obj;
    	  if(this.position < tmp.position) {
    		  return -1;
    	  }else if(this.position > tmp.position) {
    		  return 1;
    	  }
    	  return 0;
      }
  }

  public static class Settings {
    private static final double PI2 = Math.sqrt(Math.PI * 2);
    static final float DEFAULT_THRESHOLD = 4.0f;
    
    public final double[] precompute;
    public final long window;
    public final long bandwidth;
    public final int step;
    public final float threshold;
    public final int offset;
    public final String experimentType;
    public final double ncuts;
    
    public Settings(long bandwidth, long window, float threshold, int offset, double ncuts, String inputDataType){
      this.bandwidth = bandwidth;
      this.window = window;
      this.threshold = threshold;
      this.offset = offset;
      this.step = 0;
      this.ncuts = ncuts;
      
      experimentType = inputDataType;//isDNase(offset);
      precompute = precompute(window, bandwidth, ncuts);
    }
    
    public Settings(long featureLength, float threshold, int offset, double ncuts, String inputDataType){
      this.bandwidth = computeBandwidth(featureLength);
      this.window = computeOptimalWindow(bandwidth);
      this.threshold = threshold;
      this.offset = offset;
      this.step = 0;
      this.ncuts = ncuts;
      
      experimentType = inputDataType;//isDNase(offset);
      precompute = precompute(window, bandwidth, ncuts);
    }
    
    private static double sequenceNormalize(double value, double ncuts) {
  	  return ((value * 20000000d)/ncuts);
    }//Add this to all density outputs.
    
    private static double[] precompute(long window, long bandwidth, double ncuts){
      double[] precompute = new double[(int)(window + 1)];
      for(int i = 0; i < precompute.length; ++i){
        double x = i / (double)bandwidth;
        precompute[i] = sequenceNormalize( Math.exp(-(x * x) / 2) / PI2 , ncuts);
      }
      return precompute;
    }
    
    private static boolean isDNase(int offset){
    	if(offset == 0) {
    		return true;
    	}
    	return false;
    }
    
    public static long computeBandwidth(long featureLength){
      return (long)((featureLength / 2.0) / 3.0);  // 3 standard deviations
    }

    public static int computeOptimalWindow(long bandwidth){
      int i = 1400;
      
      double bw = (double)bandwidth;
      while(true){
        double x = ++i / bw;
        double v = Math.exp(-(x * x) / 2) / PI2;
        if(v < Float.MIN_VALUE)
          break;
      }
      return (i-1);
    }
  }
  
}
