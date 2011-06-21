/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2002 United States Government as represented by the
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

package jat.alg.integrators;
import jat.util.FileUtil;

import java.io.*;

/** <P>
 * The LinePrinter Class provides a way to print out integrator output data
 * to a tab delimited ASCII file. Note: remember to close the LinePrinter
 * when you are done using it.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class LinePrinter implements Printable {

	private final int[] indices;
	private PrintWriter pw;
	private final boolean printall;
	private int multiple = 1;
	private int counter = 0;
	private boolean isadaptive = false;

	/** Default constructor. Prints the entire state array to System.out.
	 */
	public LinePrinter() {
		this(null, true, null);
	}

	/** Prints a subset of the state array to System.out.
	 * @param i Integer array containing the indices of the elements of the state array to be printed.
	 */
	public LinePrinter(int[] i) {
		this(null, false, (i == null ? new int[0] : i));
	}

	/** Prints the entire state array to the user supplied directory and file.
	 * @param fname String containing the filename.
	 */
	public LinePrinter(String fname) {
        // We want to make sure we throw NullPointerException if null is passed in,
        // but still want to use "this(...)" syntax which requires there be only
        // one line.  So we sneakily call toString() on the string.
		this(fname.toString(), true, null);
	}

	/** Prints the entire state array to the user supplied directory and file.
	 * @param dir String containing the directory.
	 * @param fname String containing the filename.
	 * @param i Integer array containing the indices of the elements of the state array to be printed.
	 */
    public LinePrinter(String fname, int[] i) {
      // We want to make sure we throw NullPointerException if null is passed in,
      // but still want to use "this(...)" syntax which requires there be only
      // one line.  So we sneakily call toString() on the string.
      this(fname.toString(), false, (i == null ? new int[0] : i));
    }
    
    
    /** Prints the entire state array to the user supplied directory and file.
     * @param dir String containing the directory.
     * @param fname String containing the filename.
     * @param i Integer array containing the indices of the elements of the state array to be printed.
     * @throws IOException thrown if it can't open the file.
     */
    private LinePrinter(String fname, boolean printAllIndexes, int[] i) {
		this.printall = printAllIndexes;
		this.indices = (printAllIndexes ? null : new int[i.length]);
        if (!printAllIndexes) {
          System.arraycopy(i, 0, this.indices, 0, i.length);
        }
        
        // We open a writer to the output media.
        // If fname is null, we open a writer to System.out
        // If fname is a filename with no directory names at all
        //   (or if it specifies the default directory)
        //   we put it in the default directory.
        // Otherwise, we open a writer to the file and if not
        //   an absolute file name will create it relative to
        //   the working directory.
        
        try {
          if (fname == null) {
            this.pw = new PrintWriter(System.out, true);
          }
          else {
            OutputStream ostrm = FileUtil.openOutputFile(fname);
            this.pw = new PrintWriter(ostrm);
          }
        } catch (IOException e) {
            System.err.println("LinePrinter error opening file: " + e);
            System.err.println("All output to line printer will be lost.");
        }
	}

    /**
     * This method overrides the call to finalize() in order to close the LinePrinter when
     * the object dies.
     */
    protected void finalize() throws Throwable{
    	close();
    	super.finalize();
    }
    
	/** Closes the LinePrinter, the PrintWriter and the output file. Always remember to call this method when you are done printing!
	 */
	public void close() {
	  // close the PrintWriter
      if (pw != null) {
		pw.close();
      }
	}
	
	/** Set whether the print calls are fixed step (false) or adaptive (true).
	 * If true, the print call will disregard thinning procedures.
	 * @param a boolean parameter
	 */
	public void setIsAdaptive(boolean a){
	    this.isadaptive = a;
	}
	
	/** Print out data only every n points
	 * @param n  how often to print
	 */
	
	public void setThinning(int n) {
		multiple = n;
	}
	
	/**
	 * Return the value of the thinning.
	 * @see setThinning(double)
	 * @return value of thinning
	 */
	public int getThinning(){
	    return multiple;
	}

	/** Implements the Printable interface. This method is called once per integration step by the integrator.
	 * @param t time or independent variable.
	 * @param y state or dependent variable array.
	 */
	public void print(double t, double[] y) {
      if (pw != null) {
//		double tprint = counter * multiple;
		if (counter >= multiple || isadaptive) {
			counter = 0;
			// print the time variable
			pw.print(t + "\t");

			// print the state array
			if (printall) {
				// print all of the y array
				for (int j = 0; j < y.length; j++) {
					pw.print(y[j] + "\t");
				}
			} else {
				// check to make sure there is enough to print
				if (indices.length > y.length) {
					System.out
							.println("LinePrinter: too many elements to print");
					return;
				}
				// print the requested parts of the y array
				for (int j = 0; j < indices.length; j++) {
					pw.print(y[indices[j]] + "\t");
				}
			}

			// add the linefeed for the next line
			pw.println();
		}
		counter = counter + 1;
      }
	}

	/** Implements the Printable interface. This method is called once per integration step by the integrator.
	 * @param y double array.
	 */
	public void print(double[] y) {
	  if (pw != null) {
		// print the state array
		if (printall) {
			// print all of the y array
			for (int j = 0; j < y.length; j++) {
				pw.print(y[j] + "\t");
			}
		} else {
			// check to make sure there is enough to print
			if (indices.length > y.length) {
				System.out.println("LinePrinter: too many elements to print");
				return;
			}
			// print the requested parts of the y array
			for (int j = 0; j < indices.length; j++) {
				pw.print(y[indices[j]] + "\t");
			}
		}

		// add the linefeed for the next line
		pw.println();
      }
	}

	/** Print a string to a line
	 * @param str string to be printed
	 */
	public void println(String str) {
		pw.println(str);
	}

	/** Print a string array to a line
	 * @param str string array to be printed
	 */
	public void println(String[] str) {
      if (pw != null) {
		for (int i = 0; i < str.length; i++) {
			pw.print(str[i] + "\t");
		}
		pw.println();
      }
	}

	/** Print a string array to a line
	 * @param title string containing a title
	 * @param str string array to be printed
	 */
	public void println(String title, String[] str) {
      if (pw != null) {
		pw.print(title + "\t");
		for (int i = 0; i < str.length; i++) {
			pw.print(str[i] + "\t");
		}
		pw.println();
      }
	}

}
