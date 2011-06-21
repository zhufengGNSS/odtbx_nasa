/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2007 United States Government as represented by the
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
 * 
 * File Created on May 20, 2007
 */
 
package jat.groundstations;
 
import jat.alg.integrators.*;
import jat.math.*;
import jat.util.FileUtil;

import java.io.*;
import java.util.*;

/**
* The GroundStationList.java Class provides a way to deal with
* a list of ground stations read from a file.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class NASA_GroundStationList implements Serializable {

	private static final long serialVersionUID = 1L;

	private ArrayList<NASA_GroundStation> list = new ArrayList<NASA_GroundStation>();
	
	/** Constructor
	 */
	public NASA_GroundStationList() {
	}

	/** Constructor
	 * @param filename String containing directory and filename where the data resides
	 */
	public NASA_GroundStationList(String filename) {
		this.readFromFile(filename);
	}

	/** Add a GroundStation to the collection
	 * @param gs GroundStation object
	 */
	public void add(NASA_GroundStation gs) {
		list.add(gs);
	}
	
	/** Add a GroundStation to the collection
	 * @param index Index of GroundStation object
	 */	public void remove(int index) {
		list.remove(index);
	}

	/** Get a GroundStation out of the collection
	 * @param index Index of the Ground Station
	 * @return the FiniteBurn
	 */
	public NASA_GroundStation get(int index) {
		return (NASA_GroundStation) list.get(index);
	}

	/** Return the size of the list
	 * @return the number of ground stations in the list
	 */
	public int size() {
		return list.size();
	}

	public void print() {
		for (NASA_GroundStation next : list) {
			System.out.println(next.getSTDN()+"\t"+next.getLocation()+"\t"+next.getName()+"\t"
					+next.getLatitude()+"\t"+next.getLongitude()+"\t"+next.getHAE());
		}
	}

	/**
	 * Returns whether there is more data
	 * @param index measurement index
	 * @return true if there is more data to be read
	 */
	public boolean hasNext(int index) {
		boolean out = false;
		if (index < (this.size() - 1)) {
			out = true;
		}
		return out;
	}

	/** Read the ground station data from a tab-delimited ASCII text file.
	 * @param file filename and directory
	 */
	public void readFromFile(String file) {
		try {
			FileReader fr = new FileReader(file);
			BufferedReader in = new BufferedReader(fr);
			String line;

			// loop through the file, one line at a time
			while ((line = in.readLine()) != null) {
				StringTokenizer tok = new StringTokenizer(line, "\t");
				int total = tok.countTokens();

				// check for consistent number of columns
				if (total != 6) {
					System.out.println(
						"GroundStationList.readFromFile: Number of columns do not match");
					System.exit(-99);
				}
				// get name
				tok.nextToken();
				String stdn = tok.nextToken();
				String loc = tok.nextToken();
				String nasa = tok.nextToken();

				// get lat, lon, alt
				double[] temp = new double[3];
				for (int i = 0; i < 3; i++) {
					String token = tok.nextToken();
					temp[i] = Double.parseDouble(token)*MathUtils.DEG2RAD;
				}
				NASA_GroundStation gs = new NASA_GroundStation(stdn, loc, nasa, temp[0], temp[1], temp[2]);
				this.add(gs);
			}
			in.close();
			fr.close();
		} catch (IOException e) {
			System.err.println("Error opening:" + file);
			return;
		}
	}
	/** Read the ground station data from a NASA NDOSL ASCII text file.
	 * @param file filename and directory
	 */
	public void readFromNDOSLFile(String file) {
		try {
			FileReader fr = new FileReader(file);
			BufferedReader in = new BufferedReader(fr);
			String line;

			// loop through the file, one line at a time
			while ((line = in.readLine()) != null) {
				StringTokenizer tok = new StringTokenizer(line, "|");
				int total = tok.countTokens();

				// check for correcy number of columns
				if (total == 11) {
					// get stdn
					tok.nextToken();
					String stdn = tok.nextToken().trim();
					if (!stdn.contains("%") && !stdn.contains("&") && !stdn.contains("$") && !stdn.contentEquals("")) {
						String loc = tok.nextToken().trim();
						String nasa = tok.nextToken().trim();
						tok.nextToken();
						tok.nextToken();
						String lat_str = tok.nextToken().trim();
						double lat = this.parseDMS(lat_str);
						String lon_str = tok.nextToken().trim();
						double lon = this.parseDMS(lon_str);
						String alt_str = tok.nextToken().trim();
						double alt = Double.parseDouble(alt_str);
						NASA_GroundStation gs = new NASA_GroundStation(stdn, loc, nasa, lat, lon, alt);
						this.add(gs);
					}
				}
			}
			in.close();
			fr.close();
		} catch (IOException e) {
			System.err.println("Error opening:" + file);
			return;
		}
	}	
	
	
	/** Write the ground station data out to tab-delimited ASCII text file.
	 * @param file filename and directory
	 */
	public void sendToFile(String file) {
		LinePrinter lp = new LinePrinter(file);
		int index = 0;

		// loop through the file, one line at a time
		while (this.hasNext(index)) {
			NASA_GroundStation gs = this.get(index);
			double x = MathUtils.RAD2DEG;
			String out = gs.getSTDN()+"\t"+gs.getLocation()+"\t"+gs.getNASA()+"\t"+x*gs.getLatitude()+"\t"+x*gs.getLongitude()+"\t"+x*gs.getHAE();
			lp.println(out);
			index = index + 1;
		}
		lp.close();
	}
	
	private double parseDMS(String str){
		double radians = 0.0;
		try {
			StringTokenizer tok = new StringTokenizer(str, " ");
			double degrees = Double.parseDouble(tok.nextToken());
			degrees = degrees + Double.parseDouble(tok.nextToken()) / 60.0;
			degrees = degrees + Double.parseDouble(tok.nextToken()) / 3600.0;
			radians = degrees * MathUtils.DEG2RAD;
		} catch (NumberFormatException e) {
		}
		return radians;
	}
	
	public static void main(String [] args){
		NASA_GroundStationList gs = new NASA_GroundStationList();
    	String fs = FileUtil.file_separator();		
    	String dir = FileUtil.getClassFilePath("jat.groundstations","GroundStation");
    	String filename = "DBS_NDOSL_WGS84.txt";
    	String file = dir+fs+filename;		
		gs.readFromNDOSLFile(file);
	}


}
