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
package jat.forces.gravity;

import jat.math.MathUtils;
import jat.matvec.data.Matrix;
import jat.spacetime.*;
import jat.util.FileUtil;
import jat.eph.*;
import jat.forces.gravity.earth.*;
import jat.forces.gravity.moon.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;

/**
 * This class alows access to the gravity model data files (e.g. JGM2, JGM3, etc)
 * 
 * @author David Gaylor
 *
 */
public class GravityModel extends SphericalHarmonicGravity {

	private static final long serialVersionUID = 1L;

	private String file_separator = System.getProperty("file.separator");

    /** The GRV file to use */
	private URL grvFile = null;
	
    /** A file we output intermediate data to for test purposes. */
	private File outputFile = null;

	private int size = 0;
	
	private double gm = 0.0;
	
	private double rref = 0.0;
	
	private int nmax = 0;
	
	private int mmax = 0;

	private boolean foundSize = false;
	
	private boolean normalized;
	
	private boolean monteCarlo = false;
	
	private double gmUncertainty = 0.0;
	
	private Random rand = new Random(System.currentTimeMillis());
	
    /** Construct a gravity model from a STK gravity file included in JAT
     * This constructor assumes Earth.
     * @param n Desired degree.
     * @param m Desired order.
     * @param typ GravityModelType
     */
	public GravityModel(int n, int m, EarthGravityType typ) {
		super(n, m, new EarthFixedRef());		
		// We assume the gravity model is in the jat.forces.gravity.earth
        // package
        ClassLoader loader = getClass().getClassLoader();
        grvFile = loader.getResource("jat/forces/gravity/earth/" + typ.toString() +
            ".grv");
        System.out.println(grvFile.toString());
        // Put the output file in the earthGravity/ directory.
        String path = FileUtil.getClassFilePath("jat.forces.gravity.earth", "EarthGravityType");
		String filename = path + file_separator + "test.out";
        outputFile = new File(filename);
		initialize();
	}
	
    /** Construct a gravity model from a STK gravity file included in JAT
     * This constructor assumes Earth.
     * @param n Desired degree.
     * @param m Desired order.
     * @param typ GravityModelType
     */
	public GravityModel(int n, int m, MoonGravityType typ) {
		super(n, m, new LunaFixedRef());		
		// We assume the gravity model is in the jat.forces.gravity.earth
        // package
        ClassLoader loader = getClass().getClassLoader();
        grvFile = loader.getResource("jat/forces/gravity/moon/" + typ.toString() +
        ".grv");
        System.out.println(grvFile.toString());
        // Put the output file in the earthGravity/ directory.
        String path = FileUtil.getClassFilePath("jat.forces.gravity.moon", "MoonGravityType");
		String filename = path + file_separator + "test.out";
        outputFile = new File(filename);
		initialize();
	}	
	
    /** Construct a gravity model from any STK gravity file
     * @param n Desired degree.
     * @param m Desired order.
     * @param filepath path and filename of gravity file
     */
	public GravityModel(int n, int m, ReferenceFrame bodyFixedFrame, 
          File filepath) {
		super(n, m, bodyFixedFrame);
        try {
          grvFile = filepath.toURI().toURL(); 
          // We do the toURI() to handle spaces
        }
        catch (MalformedURLException e) {
          throw new RuntimeException("Error decoding grv filename \"" +
              filepath + "\"");
        }
        System.out.println(grvFile.toString());
        outputFile = new File(filepath.getParentFile(), "test.out");
		initialize();
	}

	/** Construct a gravity model from any STK gravity file
     * @param n Desired degree.
     * @param m Desired order.
     * @param filepath path and filename of gravity file
     */
	public GravityModel(int n, int m, ReferenceFrame bodyFixedFrame, File filepath, boolean monte) {
		super(n, m, bodyFixedFrame);
		monteCarlo = monte;
        try {
          grvFile = filepath.toURI().toURL(); 
          // We do the toURI() to handle spaces
        }
        catch (MalformedURLException e) {
          throw new RuntimeException("Error decoding grv filename \"" +
              filepath + "\"");
        }
        System.out.println(grvFile.toString());
        outputFile = new File(filepath.getParentFile(), "test.out");
		initialize();
		
	}
	

    /** Construct a gravity model from a STK gravity file in the classpath
     * @param n Desired degree.
     * @param m Desired order.
     * @param resouceName the name of a resource file 
     *  (e.g. jat/forces/moonGravity/LP165P.grv).
     *  @param monte Boolean to turn on gravity model variation for Monte Carlo
     */
    public GravityModel(int n, int m, ReferenceFrame bodyFixedRef, 
      String resourceName) {
        super(n, m, bodyFixedRef);        
        // We assume the gravity model is in the jat.forces.gravity.earth
        // package
        ClassLoader loader = getClass().getClassLoader();
        grvFile = loader.getResource(resourceName);
        System.out.println(grvFile.toString());
        // Put the output file in a temporary directory and put the 
        // resource into the name.
        int index = resourceName.lastIndexOf("/");
        String shortName = (index < 0 ?
            resourceName : resourceName.substring(index+1));
        String path = System.getProperty("java.io.tmpdir");
        outputFile = new File(path, shortName + "test.out");
        initialize();
    }
    
    /** Construct a gravity model from a STK gravity file in the classpath
     * @param n Desired degree.
     * @param m Desired order.
     * @param resouceName the name of a resource file 
     *  (e.g. jat/forces/moonGravity/LP165P.grv).
     *  @param monte Boolean to turn on gravity model variation for Monte Carlo
     */
    public GravityModel(int n, int m, ReferenceFrame bodyFixedRef, 
      String resourceName, boolean monte) {
        super(n, m, bodyFixedRef);        
        // We assume the gravity model is in the jat.forces.gravity.earth
        // package
        ClassLoader loader = getClass().getClassLoader();
        grvFile = loader.getResource(resourceName);
        System.out.println(grvFile.toString());
        // Put the output file in a temporary directory and put the 
        // resource into the name.
        int index = resourceName.lastIndexOf("/");
        String shortName = (index < 0 ?
            resourceName : resourceName.substring(index+1));
        String path = System.getProperty("java.io.tmpdir");
        outputFile = new File(path, shortName + "test.out");
        monteCarlo = monte;
        initialize();
    }    
    
	public void initialize() {
		// process the gravity file
		Matrix cs = readGrvFile(grvFile);
		// check desired gravity model size vs what we have
		if ((this.n_desired > nmax)||(this.m_desired > mmax)) {
			System.err.println("GravityModel.initialize: desired size is greater than max size from gravity file");
			System.exit(-99);
		}
		// initialize
        this.initializeGM(gm, rref);
        this.initializeCS(nmax, mmax, cs.A);
	}
	
	private Matrix readGrvFile(URL grvFile) {
		
		PrintWriter pw = null;
        if (outputFile != null)
        {
          try {
            pw = new PrintWriter(outputFile);
          } catch (FileNotFoundException e1) {
            e1.printStackTrace();
          }
        }
		
		
		Reader rdr = null;
		try {
			rdr = new InputStreamReader(grvFile.openStream());
		} catch (IOException e) {
			e.printStackTrace();
		}
		BufferedReader in = new BufferedReader(rdr);
		String line = null;
		boolean beginFound = false;
		boolean endFound = false;

		while (!beginFound) {
			try {
				line = in.readLine();
			} catch (IOException e1) {
				e1.printStackTrace();
				break;
			}
			// break line into space-delimited tokens
			String[] tokens = line.trim().split("\\s");
			
			String first = tokens[0];			

			// look for size of gravity field
			if (first.matches("Degree")) {
				nmax = Integer.parseInt(tokens[tokens.length-1]);
				size = nmax + 1;
				foundSize = true;
				System.out.println("GravityModel: Max Degree = " + nmax);
			}

			if (first.matches("Order")) {
				mmax = Integer.parseInt(tokens[tokens.length-1]);
				System.out.println("GravityModel: Max Order = " + mmax);
			}

			if (first.matches("Gm")) {
				gm = Double.parseDouble(tokens[tokens.length-1]) + 1;
				System.out.println("GravityModel: GM = " + gm);
			}

			if (first.matches("RefDistance")) {
				rref = Double.parseDouble(tokens[tokens.length-1]) + 1;
				System.out.println("GravityModel: Rref = " + rref);
			}

			if (first.matches("Normalized")) {
				if (tokens[tokens.length-1].matches("Yes")) {
					normalized = true;
				} else {
					normalized = false;
				}
				System.out.println("GravityModel: normalized = " + normalized);
			}
			
			if (first.matches("GmUncertainty")) {
				gmUncertainty = Double.parseDouble(tokens[tokens.length-1]) + 1;
				System.out.println("GravityModel: GM Uncertainty = " + gmUncertainty);
				
				// if doing Monte Carlo, apply GM uncertainty
				if (monteCarlo) {
					gm = gm + rand.nextDouble()*gmUncertainty;
					System.out.println("GravityModel: Adjusted GM = " + gm);					
				}
				
			}

			
			if (first.matches("BEGIN")) {
				beginFound = true;
				System.out.println("found BEGIN");
				if (!foundSize) {
					System.err
							.println("GravityModel: Size (Degree) of Gravity Field Not Found\n");
					System.exit(-99);
				}
			}
		}

		Matrix out = new Matrix(size, size);
		out.set(0,0,1.0);
		int nlines = 0;
		
		while (!endFound) {
			try {
				line = in.readLine();
			} catch (IOException e1) {
				e1.printStackTrace();
				break;
			}
			// break line into tokens
			StringTokenizer tokens;
			StringTokenizer tmp1 = new StringTokenizer(line," ");
			StringTokenizer tmp2 = new StringTokenizer(line,"\t");
			if(tmp1.countTokens() > tmp2.countTokens()){
				tokens = tmp1;
			}else{
				tokens = tmp2;
			}
			int ntokens = tokens.countTokens();
			nlines = nlines + 1;
//			System.out.println("number of tokens = "+ntokens);
			
			String first = " ";			
			if (ntokens > 0) {
				first = tokens.nextToken();
			} else {
				first = line.trim();
			}
			
			// look for the end
			if (first.matches("END")) {
				endFound = true;
				System.out.println("found END");
				out.print(pw);
				pw.close();
				System.out.println("Lines Processed: "+nlines);
				return out;
			}
			
			// parse the data records
			
			if (ntokens == 4) {
				if (monteCarlo){
					System.err.println("You did not supply variances on the gravity coefficients");					
				}
				int n = Integer.parseInt(first);
				int m = Integer.parseInt(tokens.nextToken());
				
				// compute unnormalization factor
				double nf = 1.0;				
				if (normalized) nf = normFactor(n, m);
				
				// obtain the c and s values
				double c = nf * Double.parseDouble(tokens.nextToken());
				double s = nf * Double.parseDouble(tokens.nextToken());
//				System.out.println(ntokens+":"+n+":"+m+":"+c+":"+s);
				
				// put c and s into the correct places in the cs matrix
				out.set(n, m, c);
				if (m > 0) out.set(m-1, n, s);				
			}
			// parse the data records
			if (ntokens == 6) {
				int n = Integer.parseInt(first);
				int m = Integer.parseInt(tokens.nextToken());
				
				// compute unnormalization factor
				double nf = 1.0;				
				if (normalized) nf = normFactor(n, m);
				
				// obtain the c and s values
				double c = nf * Double.parseDouble(tokens.nextToken());
				double s = nf * Double.parseDouble(tokens.nextToken());
//				System.out.println(ntokens+":"+n+":"+m+":"+c+":"+s);
				
				// obtain the c and s uncertainties
				double cSigma = nf * Double.parseDouble(tokens.nextToken());
				double sSigma = nf * Double.parseDouble(tokens.nextToken());
				
				// if doing Monte Carlo, apply the uncertainties
				double adjustC = rand.nextDouble()*cSigma;
				double adjustS = rand.nextDouble()*sSigma;
				c = c + adjustC;
				s = s + adjustS;
//				System.out.println(ntokens+":"+n+":"+m+":"+adjustC+":"+adjustS);
				
				// put c and s into the correct places in the cs matrix
				out.set(n, m, c);
				if (m > 0) out.set(m-1, n, s);				
			}
			
		}
		
		out.print(pw);
        pw.close();
        
        try {
          in.close();
          rdr.close();
        }
        catch (IOException ioe) {
          ioe.printStackTrace();
        }
		return out;
	}
	
	private double normFactor (int n, int m) {
		double nn = new Integer(n).doubleValue();
		double mm = new Integer(m).doubleValue();
		double nmmfact = MathUtils.factorial(nn - mm);
		double npmfact = MathUtils.factorial(nn + mm);
		double  delta = 0.0;
		if (m == 0) delta = 1.0;
		double num = nmmfact*(2.0*nn+1.0)*(2.0-delta);
		double out = 0.0;
		if (npmfact == 0.0) {
			System.out.println("GravityModel.normFactor: denominator = 0, n= "+n+" m = "+m+" npmfact = "+npmfact);
		} else {
			out = Math.sqrt(num/npmfact);
		}
		return out;
	}

	public static void main(String args[]) {
		BodyCenteredInertialRef ref = new BodyCenteredInertialRef(DE405_Body.GEOCENTRIC_MOON);
		GravityModel x = new GravityModel(165, 165, ref, "jat/forces/gravity.moon/LP165P-Cov.grv", true);
		x.printParameters();
	}



}

