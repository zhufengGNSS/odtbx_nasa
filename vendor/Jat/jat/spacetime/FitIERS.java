/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2005 United States Government as represented by the
 * administrator of the National Aeronautics and Space Administration.
 * All rights reserved.
 *
 * This file is part of JAT. JAT is free software; you can
 * redistribute it and/or modify it under the terms of the
 * NASA Open Source Agreement, version 1.3 or later.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * NASA Open Source Agreement for more details.
 *
 * You should have received a copy of the NASA Open Source Agreement
 * along with this program; if not, write to the NASA Goddard
 * Space Flight Center at opensource@gsfc.nasa.gov.
* 
*/
package jat.spacetime;

import jat.util.FileUtil;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.StringTokenizer;

/**
 * International Earth Rotation and Reference Service parameters.  Parses and calculates
 * the polar motion for the Celestial Ephemeris Pole and the offset UT1-UTC.<hr>
 * File format for the polar motion datafile is as follows:<br>
 * First line contains a single number indicating the number of lines in the file
 * <ol>Each subsequent line contains the following items and their units, whitespace-delimited:
 * <li>time (Modified Julian Date in UTC)</li>
 * <li>xpole (arcsec)</li>
 * <li>xpole_error (arcsec, <i>currently unused</i>)</li>
 * <li>ypole (arcsec)</li>
 * <li>ypole_error (arcsec, <i>currently unused</i>)</li>
 * <li>UT1-UTC (seconds)</li>
 * <li>time_error (seconds, <i>currently unused</i>)</li>
 * </ol>
 * <br>
 * <b>Each time value should be exactly one day apart.</b>  This is very important
 * because the internal algorithms compare and search on this assumption.  Be sure to
 * set the first and last entries in the datafile to appropriate values given the
 * behavior of search() to clamp to the nearest values for times outside the data set.<br>
 * <hr>
 * Note: The error values are not currently used.  This class only processes
 * the time, xpole, ypole, and UT1-UTC items.  If the error are not available, 
 * format their locations with zeros.
 * <hr>
 * <b>WARNING: This class holds the loaded data as static data - it applies to
 * all instances.  Specifying a different data file anywhere in the JVM will
 * cause <i>all instances of FitIERS</i> to behave differently.</b>
 * 
 * @author Richard C. Page III
 * @author Allen Brown
 */
public class FitIERS implements Serializable {

	private static final long serialVersionUID = 1L;

	/* 
	 * Developer Notes:
	 * 
	 * TODO: Implement proper thread-safe data caching mechanism which handles multiple data sources.
	 * 
	 * This class has two data caching mechanisms in its design.
	 * 
	 * First, it remembers the last set of data calculated in order to avoid
	 * a search through the IERS dataset.  Usually multiple calls are made
	 * to search(double mjd) with the same or similar mjd arguments.  This is
	 * appropriate at the instance level.
	 * 
	 * Second, the data set itself should be cached from the filesystem.  Many
	 * clients of this class create and destroy instances frequently.  Caching
	 * the filesystem data for all instances would reduce the data read hit to
	 * performance.  However, as the to do mentions above, this class currently
	 * does not have a proper filesystem caching implementation.  Currently
	 * any instance which reads from a different file will overwrite all data
	 * used by all instances.  See jat.eph.DE405 for an example of a safer
	 * mechanism.
	 * 
	 * The minimum effort to fix this was: Make process() non-public and
	 * make it static & final since its called from constructors;  Coordinate
	 * all cached data and its data source (filename, prev_filename); make all
	 * modifier and accessor methods synchronize on the static data.
	 * 
	 * TODO: Implement constructor argument checks and data checks...
	 * properly throw exceptions from the constructor and handle in calling 
	 * classes.  Ensure that a bad parse doesn't impact other instances/threads.
	 * 
	 * First, if the data file, assumed or specified, can't be found or read, 
	 * then any prior data cached by this class is silently preserved.  No 
	 * exceptions are thrown by the constructors so the user doesn't know 
	 * they're not using the data set they intended.  Unfortunately this is 
	 * kept right now to avoid changing client code.
	 * 
	 * Secondly, all data points in the data file are assumed by search() to be
	 * exactly 1.0 JD apart.  This assumption was in the original algorithm and 
	 * replicated in the current implementation so that behavior isn't changed.
	 * However, this assumption is never checked when reading in the data file (!).
	 * An earlier version of the jat/spacetime/EOP.dat file violated this 
	 * assumption with its last data point and that caused problems in 
	 * search()'s "quick search" algorihtm using a linear fit across the entire
	 * data set.  Later this same violation caused problems with times around
	 * the last valid data point.  Ultimately this "hole" in the input data
	 * has to be ultimately disallowed - there are too many cases for this logic
	 * to break.  If someone needs EOP data in different increments then the
	 * data storage, search logic, and any caching mechanisms should be re-
	 * designed from scratch.
	 * 
	 */
	
	/*
	 * Static Members
	 */
	
	/**
	 *  Comparator for performing the binary search.  This has a special
	 *  compare() operator for searching the FitIERS.mdj array.
	 */
	protected static class MDJTimeComparator implements Comparator<Double> {
		
		/**
		 * This special compare operator is for finding a time in the
		 * FitIERS.mjd array. <br>
		 * The two times are "equal" if t1 and t2 are not more than 1.0 apart. <br>
		 * t1 is greater than t2 if t1 >= t2 + 1.0 <br>
		 * t2 is greater than t1 if t2 >= t1 + 1.0
		 * @param t1 time to compare
		 * @param t2 time to compare
		 * @return A negative integer, zero, or a positive integer as the first argument is less than, equal to, or greater than the second. 
		 */
		public int compare(Double t1, Double t2) {
			final double tmp = t1-t2;
			if (tmp <= -1.0) {
				return -1;
			}
			else if (tmp >= 1.0) {
				return 1;
			} else {
				return 0;
			}
		}
	}
	
	/**
	 * The TimeComparator instance used by this class to compare times when
	 * searching the mjd array.
	 */
	protected static final MDJTimeComparator timecomp = new MDJTimeComparator();
	
	/**
	 * A convenient non-null Object that can be used for synchronization before
	 * any static data is populated and cached or accessed.
	 */
	protected static final Object syncObject;
	
	/**
	 * Cached data:
	 * Modified Julian Date Doubles array of the cached time data points. 
	 * This applies to all instances.
	 */
	/* 
	 * We'll use a Double[] here in order to be compatible with java.util's
	 * binarSearch routines.  Its built once, when we read the file, so there
	 * is low overhead when searching.
	 */
	protected static Double[] mjd;
	
	/**
	 * Cached data:
	 * The value of the first element of mjd[].
	 * This applies to all instances.
	 * Legacy: only used for debugging, left in case the debugging statements
	 * are re-enabled.
	 */
	protected static double mjd_1=0.0;
	
    /**
     * Cached data:
	 * The IERS EOP file data that is cached in this instance:
	 * xpole, ypole, and UT1-UTC.
	 * This applies to all instances.
     */
    protected static double[][] eop;
    
    /**
     * Cached data:
	 * Whether the cache has been initialized (true) or not (false).
	 * This applies to all instances.
     */
    protected static boolean initialized = false;
	
    /**
     * Cached data:
     * Filename of the current IERS EOP file (Bulletin A) cached data.
	 * This applies to all instances.
     */
    protected static String filename;
    
    /**
     * Cached data:
     * Filename of the previous IERS EOP file (Bulletin A) cached data, for 
     * comparison.
	 * This applies to all instances.
     */
    protected static String prev_filename = "UNINITIALIZED";
    
    /*
     * Static initialization block:
     */
    static {
    	// the syncObject just needs a valid non-null Object
    	syncObject = new String("FitIERS");
    }
    

	/**
	 * Returns the earliest time in the cached dataset.  This method is static
	 * so the data can be queried without disturbing it with a constructor.
	 * @return The earliest time in the dataset (Modified Julian Date UTC),
	 * or 0.0 if uninitialized.
	 */
	public static double getEarliestTime() {

		synchronized (syncObject) {
			if (initialized) {
				
				// Search forwards for a non-zero eop element, we shouldn't
				// have to go far.
				for (int i=0; i<mjd.length-1; i++) {
					if (eop[i][0] + eop[i][1] + eop[i][2] != 0.0) {
						return mjd[i];
					}
				}
			}
			// default, either not initialized or all data is zeros
			return 0.0; 
		} // synchronized
	}

	/**
	 * Returns the latest valid time in the cached dataset.  This method is 
	 * static so the data can be queried without disturbing it with a 
	 * constructor.
	 * <br>
	 * A valid time is a time that contains useful (non-zero) data for any EOP 
	 * parameter.  Invalid times may be present at the end of the data file
	 * but they contain EOP data that is all zeros.
	 * @return The latest time in the dataset (Modified Julian Date UTC),
	 * or 0.0 if uninitialized.
	 */
	public static double getLatestTime() {

		synchronized (syncObject) {
			if (initialized) {
				
				// Search backwards for a non-zero eop element, we shouldn't
				// have to go far.
				for (int i=mjd.length-1; i>=0 ;i--) {
					if (eop[i][0] + eop[i][1] + eop[i][2] != 0.0) {
						return mjd[i];
					}
				}
			}
			// default, either not initialized or all data is zeros
			return 0.0; 
		} // synchronized
	}
	
	/*
	 * Instance Members:
	 */
	
	/**
	 * When true, enables the checking of the instance cache data, the
	 * last_mjd, last_out, and last_ind private members that can help
	 * eliminate a full search based on the last answer.
	 * Used inside search().
	 */
	private boolean instance_cache_enabled = true;
    
    /**
     * Instance cache of the last result (to avoid searching).
     * This is the last MJD time presented to search(). 
     */
	private double last_mjd = 0.0;
	
	/**
	 * Instance cache of the last result (to avoid searching).
	 *  The last returned results from search() based on last_mjd.
	 */
	private double[] last_out = new double[3];
	
	/**
	 * Instance cache of the last result (to avoid searching).
	 *  The index that supported last_out.
	 */
	private int last_ind = 0;
	
	/**
	 * Running total of how many times the instance cache has been
	 * hit during search() for this instance's lifetime.  This
	 * initializes to zero and only increments when the cache is
	 * used to return a result in search().  Not incremented when
	 * the cache is disabled.  NOTE: This value will wrap without
	 * checking.  
	 */
	private int cache_hits = 0;
	
	/**
	 * Running total of how many times the instance cache has been
	 * missed during search() for this instance's lifetime.  This
	 * initializes to zero and only increments when the cache is
	 * can't return a result in search().  Not incremented when
	 * the cache is disabled.  NOTE: This value will wrap without
	 * checking.  
	 */
	private int cache_misses = 0;

    /**
     * Default constructor 
     */
    public FitIERS() {
        String fs = FileUtil.file_separator();
		String directory;
		try{
		    directory = FileUtil.getClassFilePath("jat.spacetime","FitIERS")+fs;
		}catch(Exception ne){
		    directory = "C:/Code/Jat/jat/spacetime/";
		}
		//filename = directory+"iers.dat";
		
		try {
			// this may throw a IOException but since we don't yet
			// want to update every single class that uses FitIERS we'll just
			// dump a stack trace here.
			
			process(directory+"EOP.dat");
			
		// We have two catches because we should differentiate between a missing file
		// and an unparsable file. 
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.err.println("The file for the FitIERS default constructor is missing.  Check your installation!");
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("The file for the FitIERS default constructor can't be parsed.  Check your installation!");
		}
		
		// initialize the cached instance data to the first data point
		last_mjd = mjd[0];
		last_out[0] = eop[0][0];
		last_out[1] = eop[0][1];
		last_out[2] = eop[0][2];
		last_ind = 0;
    }
    
    /**
     * Constructor
     * @param fname Filename of the IERS reference data (eg "iers.dat").
     */
	public FitIERS(final String fname) {
		
		try {
			// this may throw an IOException but since we don't yet
			// want to update every single class that uses FitIERS we'll just
			// dump a stack trace here.
			
			process(fname);
			
		// We have two catches because we should differentiate between a missing file
		// and an unparsable file. 
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.err.println("The file for the FitIERS String constructor is missing.  Check your installation!");
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("The file for the FitIERS String constructor can't be parsed.  Check your installation!");
		}
		
		// initialize the cached instance data to the first data point
		last_mjd = mjd[0];
		last_out[0] = eop[0][0];
		last_out[1] = eop[0][1];
		last_out[2] = eop[0][2];
		last_ind = 0;
	}
	
	/**
	 * Parse the file into data.
	 * @throws IOException If the given file is not found or if there is a 
	 * problem reading the file.
	 */
	final static protected void process(final String newFileName) throws IOException {
		
		synchronized (syncObject) {
			
			// did we actually update the data?
			boolean updated = false;

			//System.out.println("Processing File: " + filename);
			if(!prev_filename.equals(newFileName) || !initialized){
				
				FileReader fr;
				BufferedReader in = null;
				
				// This may throw a FileNotFoundException:
				fr = new FileReader(newFileName);
				in = new BufferedReader(fr);			

				ArrayList<String> lineCollection = new ArrayList<String>();
	
				// read in the file
				boolean eof = false;
				while (!eof) {
					String line;
					try {
						// the readLine may throw an IOException
						if ((line = in.readLine()) == null) {
							eof = true;
						} else {
							// add to the collection
							lineCollection.add(line);
						}
					} catch (IOException e1) {
						try {
							in.close();
						} catch (IOException e) {
							// ignore an exception on our clean-up attempt
						}
						try {
							fr.close();
						} catch (IOException e) {
							// ignore an exception on our clean-up attempt
						}
						e1.printStackTrace();
						throw e1; // intentional re-throw
					}
				}
				// these may throw an IO Exception:
				in.close();
				fr.close();
	
				int n = lineCollection.size() - 1;
	//			int n = 28;
	//			int n = 48;
				//	int n = 2191;
	//			double[] mjd = new double[n];
				mjd = new Double[n];
				double[] y = new double[n];
	//			Matrix hth = new Matrix(3, 3);
	//			Matrix ht = new Matrix(3, n);
	//			Matrix ata = new Matrix(5, 5);
	//			VectorN aty = new VectorN(5);
	//			VectorN atx = new VectorN(5);
	
				String str;
				StringTokenizer tok;
	
				eop = new double[n][3];
	
				for (int i = 0; i < n; i++) {
	//				str = (String) lineCollection.get(i);
	//				tok = new StringTokenizer(str, " ");
	//				mjd[i] = Double.parseDouble(tok.nextToken());
	//				double pmx = Double.parseDouble(tok.nextToken());
	//				double pmy = Double.parseDouble(tok.nextToken());
	//				y[i] = Double.parseDouble(tok.nextToken());
	
					//***Read EOP.dat
					//System.out.println("parse:  "+filename);
	
					//str = (String) lineCollection.get(i+1860);//1870);
					str = (String) lineCollection.get(i+1);//all
	//				System.out.println("line "+i+":  "+str);
					tok = new StringTokenizer(str, " ");
					mjd[i] = Double.parseDouble(tok.nextToken());
					double pmx = Double.parseDouble(tok.nextToken());
					eop[i][0] = pmx;
					tok.nextToken();
					double pmy = Double.parseDouble(tok.nextToken());
					eop[i][1] = pmy;
					tok.nextToken();
					y[i] = Double.parseDouble(tok.nextToken());
					eop[i][2] = y[i];
					tok.nextToken();
	//				}//*DEBUG
					//System.out.println(""+mjd[i]+"  "+pmx+"  "+pmy+"  "+y[i]);
	
	//				double x = 0.0;
	//				double dtleap = 0.0;
	//				if (i > 0) {
	//				//x = mjd[i] - mjd[i - 1];
	//				//x += (mjd[i]-mjd[i-1]);
	//				x = i * 1.0; // one day increment in data file
	//				if (Math.abs(y[i] - y[i - 1]) > 0.5) {
	//				dtleap = 1.0;
	//				}
	//				y[i] = y[i] - dtleap;
	//				}
	
	//				double x2 = x * x;
	//				double x3 = x2 * x;
	//				double x4 = x3 * x;
	//				hth.set(2, 2, (hth.get(2, 2) + 1.0));
	//				hth.set(2, 1, (hth.get(2, 1) + x));
	//				hth.set(1, 2, (hth.get(2, 1)));
	//				hth.set(2, 0, (hth.get(2, 0) + x2));
	//				hth.set(1, 1, (hth.get(2, 0)));
	//				hth.set(0, 2, (hth.get(2, 0)));
	//				hth.set(1, 0, (hth.get(1, 0) + x3));
	//				hth.set(0, 1, (hth.get(1, 0)));
	//				hth.set(0, 0, (hth.get(0, 0) + x4));
	
	//				ht.set(0, i, x2);
	//				ht.set(1, i, x);
	//				ht.set(2, i, 1.0);
	//				// process polar motion data
	//				tp = mjd[0] - 1.0;
	//				double num = 2.0 * (mjd[i] - tp) * MathUtils.PI;
	//				double ca = num / 365.25;
	//				double cc = num / 435.0;
	//				double cosa = Math.cos(ca);
	//				double sina = Math.sin(ca);
	//				double cosc = Math.cos(cc);
	//				double sinc = Math.sin(cc);
	
	//				ata.set(0, 0, (ata.get(0, 0) + 1.0));
	//				ata.set(1, 0, (ata.get(1, 0) + cosa));
	//				ata.set(2, 0, (ata.get(2, 0) + sina));
	//				ata.set(3, 0, (ata.get(3, 0) + cosc));
	//				ata.set(4, 0, (ata.get(4, 0) + sinc));
	//				ata.set(0, 1, (ata.get(1, 0)));
	//				ata.set(0, 2, (ata.get(2, 0)));
	//				ata.set(0, 3, (ata.get(3, 0)));
	//				ata.set(0, 4, (ata.get(4, 0)));
	//				ata.set(1, 1, (ata.get(1, 1) + cosa * cosa));
	//				ata.set(2, 1, (ata.get(2, 1) + cosa * sina));
	//				ata.set(3, 1, (ata.get(3, 1) + cosa * cosc));
	//				ata.set(4, 1, (ata.get(4, 1) + cosa * sinc));
	//				ata.set(1, 2, (ata.get(2, 1)));
	//				ata.set(1, 3, (ata.get(3, 1)));
	//				ata.set(1, 4, (ata.get(4, 1)));
	//				ata.set(2, 2, (ata.get(2, 2) + sina * sina));
	//				ata.set(3, 2, (ata.get(3, 2) + sina * cosc));
	//				ata.set(4, 2, (ata.get(4, 2) + sina * sinc));
	//				ata.set(2, 3, (ata.get(3, 2)));
	//				ata.set(2, 4, (ata.get(4, 2)));
	//				ata.set(3, 3, (ata.get(3, 3) + cosc * cosc));
	//				ata.set(3, 4, (ata.get(3, 4) + cosc * sinc));
	//				ata.set(4, 3, (ata.get(3, 4)));
	//				ata.set(4, 4, (ata.get(4, 4) + sinc * sinc));
	
	//				atx.set(0, (atx.get(0) + pmx));
	//				atx.set(1, (atx.get(1) + pmx * cosa));
	//				atx.set(2, (atx.get(2) + pmx * sina));
	//				atx.set(3, (atx.get(3) + pmx * cosc));
	//				atx.set(4, (atx.get(4) + pmx * sinc));
	//				aty.set(0, (aty.get(0) + pmy));
	//				aty.set(1, (aty.get(1) + pmy * cosa));
	//				aty.set(2, (aty.get(2) + pmy * sina));
	//				aty.set(3, (aty.get(3) + pmy * cosc));
	//				aty.set(4, (aty.get(4) + pmy * sinc));
	
	//				}
	//				//* Get inverse
	//				Matrix hthinv = hth.inverse();
	//				//* Get hty
	//				VectorN ymat = new VectorN(y);
	//				VectorN hty = ht.times(ymat);
	
	//				//* solve
	//				u = hthinv.times(hty);
	//				//* invert ata (least squares fit for x- and y- parameters
	//				Matrix atainv = ata.inverse();
	//				//* multiply the inverses
	//				ax = atainv.times(atx);
	//				ay = atainv.times(aty);
				}
				mjd_1 = mjd[0];
				updated = true;
			}
		
			if (updated) {
				// note this cached data has been initialized
				initialized = true;
				
				// update the cached filename since we're reading new data
				filename = newFileName;
				prev_filename = filename;		
			}
			
		} // synchronized
	}

	/**
	 * This method scans the EOP dataset to determine the corresponding 
	 * parameters of the given MJD.  This method linearly interpolates
	 * the data to find the parameters. <br>
	 * This method tries to avoid searching by comparing to the last
	 * answer returned by this instance's search() call.  If that doesn't
	 * match then a limited search is attempted.  Finally, a full binary
	 * searh is performed to find the matching time.  <br>
	 * Out-of-range data is clamped to the first or last data point in the set
	 * and bypasses cache checks.  Updates the cache hit and miss counters
	 * if the cache is enabled.
	 * @param mjd_arg The Modified Julian Date of interest
	 * @return [xpole ypole UT1-UTC].  Clamps to first or last dataset if 
	 * mjd_arg is out-of-range for the dataset.
	 * @throws RuntimeException Deliberately thrown when the underlying data 
	 * doesn't support the search.
	 */
	public double[] search(double mjd_arg) throws RuntimeException {

		synchronized (syncObject) {
			//		System.out.println("***SEARCH***");

			final int n = mjd.length;
			int out = n; // n is out of range for an index so here it indicates "not found"
			double[] output = new double[3];

			/*
			 * Argument checks for out-of-range conditions:
			 */

			// Clamp for out-of-range conditions: before data
			if( mjd_arg < mjd[0]) {
				output[0] = eop[0][0];
				output[1] = eop[0][1];
				output[2] = eop[0][2];
				return output;
			}

			// Clamp for out-of-range conditions: after data
			if( mjd_arg >= mjd[n-1]) {
				output[0] = eop[n-1][0];
				output[1] = eop[n-1][1];
				output[2] = eop[n-1][2];
				return output;
			}

			// instance cache logic checks
			if (this.instance_cache_enabled) {

				// check for exact match from last call
				if(mjd_arg == last_mjd){
					output[0] = last_out[0];
					output[1] = last_out[1];
					output[2] = last_out[2];
					this.cache_hits++;
					return output;
				}

				// check for same data range from last call:
				final double tmp = mjd_arg - mjd[last_ind];
				if ( (tmp >= 0.0) && (last_ind + 1 < n) && (mjd_arg < mjd[last_ind+1]) ) {
					out = last_ind;
					this.cache_hits++;
				}
				else {
					this.cache_misses++;
				}
			} // instance cache checks
			

			/*
			 * No match to prior data then do a search.
			 * First try a quick search using a linear fit approximation.
			 */
			if (out == n) {

				//System.out.println("mjd0: "+mjd[0]+"   mjdarg: "+mjd_arg);
				//System.out.println("mjd2: "+mjd[n/2]+"   mjdf: "+mjd[n-1]);
				//System.out.println("mjd_arg: "+mjd_arg);
				//System.out.println("length: "+mjd.length);
				//System.out.println("den:    "+(mjd[mjd.length-1]-mjd[0]));

				final double timespan;  // the approx. timespan of the mjd dataset
				final double timecount;  // the corresponding count of the mjd data
				if (n > 10) {
					// If possible, ignore the first few and last few data points,
					// which could be outliers from a linear progression.
					// This used to be a significant problem when the data file
					// was loosely allowed to violate the 1.0 day separation
					// near the end of the file.  This caused too many problems 
					// but we'll keep this "safety check" just in case for a while.
					timespan = mjd[mjd.length-5]-mjd[4];
					timecount = mjd.length-8;
				}
				else {
					timespan = mjd[mjd.length-1]-mjd[0];
					timecount = mjd.length;
				}

				// Estimate the approximate index, ai, by a linear fit over the
				// first and last times.
				final double ref = (timecount/timespan);
				int ai = (int)( (mjd_arg-mjd[0])*ref-5); // fudge by backing up 5 spots
				
				// sometimes the approximations and truncation can cause 
				// ai to be negative, fix it here
				if (ai < 0) {
					ai = 0;
				}

				if (ai >= 0) {
					// Check from ai forward 10 positions, only if we don't run off the 
					// end of the array
					int af = ((ai + 10) > n-1) ? n-1 : (ai + 10);

					for (int i=ai; i<af ;i++) {
						//System.out.println("MATLAB: FitIERS search: "+i+"  ai: "+ai+"  ref: "+ref);
						if((mjd_arg - mjd[i])<1.0){
							out = i;
							//System.out.println("mjd1: "+mjd[i]+"  "+eop[i][0]);
							break;
						}
					}
				}
			} // data range check
			
			/*
			 * If the quick search didn't work.... bring out the big guns.
			 */
			if (out == n) {
				/* 
				 * Implement a binary search using our own Comparator to search
				 * for the appropriate closest data point. 
				 */
				out = Collections.binarySearch(Arrays.asList(mjd), mjd_arg, timecomp);
			}
			
			// Check the solution before using
			if ((out == n) || (out < 0)) {
				// This should only occur if there is a problem with the data...
				// usually with the 1.0 day increment assumption.
				
				throw new RuntimeException("FitIERS.search() could not search its dataset (" 
						+ FitIERS.filename + ") for " + mjd_arg 
						+ " MJD.\n Perhaps the dataset violates the 1.0 day increment assumption?");
			}

			//* linear interpolation
			output[0] = (eop[out+1][0]-eop[out][0])/(mjd[out+1]-mjd[out])*(mjd_arg-mjd[out])+eop[out][0];
			output[1] = (eop[out+1][1]-eop[out][1])/(mjd[out+1]-mjd[out])*(mjd_arg-mjd[out])+eop[out][1];
			output[2] = (eop[out+1][2]-eop[out][2])/(mjd[out+1]-mjd[out])*(mjd_arg-mjd[out])+eop[out][2];

			// Cache this result:
			last_mjd = mjd[out];
			last_out[0] = output[0];
			last_out[1] = output[1];
			last_out[2] = output[2];
			last_ind = out;

			return output;

		} // synchronized
	}
	
	/**
	 * Enables the instance cache checking during the search().
	 * This may eliminate a full search of the IERS data or reduce
	 * the search length based on the prior call to this instance.
	 */
	public void enableInstanceCache() {
		this.instance_cache_enabled = true;
	}
	
	/**
	 * Disables the instance cache checking during the search().
	 * This makes search() always perform a full search.
	 */
	public void disableInstanceCache() {
		this.instance_cache_enabled = false;
	}
	
	/**
	 * Returns the status of the instance cache check during search().
	 * If enabled, this may eliminate a full search of the IERS data or reduce
	 * the search length based on the prior call to this instance.  If disabled
	 * then search() always performs a full search.
	 * @return True if the instance cache is enabled, false if the check is
	 * disabled.
	 */
	public boolean isInstanceCacheEnabled() {
		return this.instance_cache_enabled;
	}
	
	/**
	 * @return Running total of how many times the instance cache has been
	 * hit during search() for this instance's lifetime.  This
	 * initializes to zero and only increments when the cache is
	 * used to return a result in search().  Not incremented when
	 * the cache is disabled.  NOTE: This value will wrap without
	 * checking.  
	 */
	public int getCacheHits() {
		return cache_hits;
	}

	/**
	 * @return Running total of how many times the instance cache has been
	 * missed during search() for this instance's lifetime.  This
	 * initializes to zero and only increments when the cache is
	 * can't return a result in search().  Not incremented when
	 * the cache is disabled.  NOTE: This value will wrap without
	 * checking.  
	 */
	public int getCacheMisses() {
		return cache_misses;
	}
	
//	/**
//	 * Fit the x and y poles to the given Julian Date
//	 * @param mjd Modified Julian Date
//	 */
//	public void fit(double mjd){
//	    // Find the nearest UT1 value
////	    for(int i=0; i<this.mjd.length; i++){
////	        if(this.mjd[i] < mjd){
////	            UT1 = this.mjd[i];
////System.out.println("fit - UT1: "+UT1+"  mjd: "+mjd);	        
////	        }else
////	            break;
////	    }
//	    double timeDiff = mjd-tp;
//	    //double timeDiff = 5.000006;
//	    double A1 = 2*MathUtils.PI/365.25*(timeDiff);
//	    double C1 = 2*MathUtils.PI/435.0*(timeDiff);
//	    pmx = ax.get(0)+ax.get(1)*Math.cos(A1)+ax.get(2)*Math.sin(A1)+
//	    	ax.get(3)*Math.cos(C1)+ax.get(4)*Math.sin(C1);
//	    pmy = ay.get(0)+ay.get(1)*Math.cos(A1)+ay.get(2)*Math.sin(A1)+
//    		ay.get(3)*Math.cos(C1)+ay.get(4)*Math.sin(C1);
////	    pmx = ax.get(4)+ax.get(3)*Math.cos(A1)+ax.get(2)*Math.sin(A1)+
////	    	ax.get(1)*Math.cos(C1)+ax.get(0)*Math.sin(C1);
////	    pmy = ay.get(4)+ay.get(3)*Math.cos(A1)+ay.get(2)*Math.sin(A1)+
////			ay.get(1)*Math.cos(C1)+ay.get(0)*Math.sin(C1);
//	    diffUT1UTC = u.get(2)+u.get(1)*(mjd-mjd_1)+u.get(0)*(mjd-mjd_1)*(mjd-mjd_1); 
//	}

		
//	public void fit_value(double mjd_arg){
//	    double[] fit = search(mjd_arg);
//	    pmx = fit[0];
//	    pmy = fit[1];
//	    diffUT1UTC = fit[2];
//	}
	
//	public void print(){
//	    if(u!=null && ax !=null && ay!=null){
//	        System.out.println(""+mjd_1+"   "+u.get(0)+"   "+u.get(1)+"   "+u.get(2));
//	        System.out.println(""+tp+"   "+ax.get(0)+"   "+ax.get(1)+"   "+ax.get(2)+"   "+ax.get(3)+"   "+ax.get(4));
//	        System.out.println(""+tp+"   "+ay.get(0)+"   "+ay.get(1)+"   "+ay.get(2)+"   "+ay.get(3)+"   "+ay.get(4));
//	    }
//	    System.out.println("***");
//	    System.out.println("x: "+pmx+"  y: "+pmy+"  diff: "+diffUT1UTC);
//	}
//	
//	public void set_update_interval(double val){
//	    this.update_interval = val;
//	}
//	
//	public double get_update_interval(){return this.update_interval;}
//	
//	/**
//	 * Returns the Earth's x pole;
//	 * @return
//	 */
//	public double getX(){ return pmx;}
//	/**
//	 * Returns the Earth's y pole;
//	 * @return
//	 */
//	public double getY(){ return pmy;}
//	/**
//	 * Returns the time difference between UT1 and UTC;
//	 * @return
//	 */
//	public double getDUT1(){ return diffUT1UTC;}
	
//	private void setCoeff(double[] u, double[] ax, double[] ay){
//	    this.u = new VectorN(u);
//	    this.ax = new VectorN(ax);
//	    this.ay = new VectorN(ay);
//	}
	
//	public static void main(String[] args) {
//
//	}
//	
//	public void printDiff(double[] u, double[] a){
//	    if(u!=null && a!=null){
//	        System.out.println("DIFFERENCE");
//	        System.out.println(""+(u[0]-this.u.get(2))+"  "+(u[1]-this.u.get(1))+"  "+(u[2]-this.u.get(0)));
//	        System.out.println(""+(a[0]-this.ax.get(0))+"  "+(a[1]-this.ax.get(1))+"  "+(a[2]-this.ax.get(2))+"  "+(a[3]-this.ax.get(3))+"  "+(a[4]-this.ax.get(4)));
//	        System.out.println(""+(a[5]-this.ay.get(0))+"  "+(a[6]-this.ay.get(1))+"  "+(a[7]-this.ay.get(2))+"  "+(a[8]-this.ay.get(3))+"  "+(a[9]-this.ay.get(4)));
//	    }
//	}
}
