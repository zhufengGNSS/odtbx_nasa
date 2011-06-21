/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2003 United States Government as represented by the
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
 * File Created on May 19, 2003
 */

package jat.gps;
import jat.alg.estimators.*;

import java.io.*;
import java.util.*;

/**
* The GPS_MeasurementList.java Class provides a way to deal with
* a list of GPS measurements read from a file.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class GPS_MeasurementList implements MeasurementData, Serializable {

	private ArrayList list = new ArrayList();
	
	private boolean carrier = false;

	/** Constructor
	 */
	public GPS_MeasurementList() {
	}

	/** Constructor
	 * @param filename String containing directory and filename where the data resides
	 */
	public GPS_MeasurementList(String filename) {
		this.readFromFile(filename);
	}
	
	/**
	 * Set Measurement Flag. Determines which type of measurement is returned
	 * by the function z.
	 * @param boolean true = carrier phase measurement, false = pseudorange
	 */
	public void setMeasFlag(boolean flag){
		this.carrier = flag;
	}

	/** Add a GPS Measurement to the collection
	 * @param meas GPS_Measurement object
	 */
	public void add(GPS_Measurement meas) {
		list.add(meas);
	}

	/** Get a GPS_Measurement out of the collection
	 * @param index Index of the measurement
	 * @return the GPS Measurement
	 */
	public GPS_Measurement get(int index) {
		return (GPS_Measurement) list.get(index);
	}

	/** Return the size of the list
	 * @return the number of measurements in the list
	 */
	public int size() {
		return list.size();
	}

	/** Returns the range measurement
	 * @param index int containing the measurement index
	 * @return double containing the measurement corresponding to the index.
	 */
	public double z(int index) {
		GPS_Measurement meas = this.get(index);
		if (carrier) {
			return meas.phaseRange();
		} else {
			return meas.range();
		}
	}

	/**
	 * Returns the time of the measurement
	 * @param index int containing the measurement index
	 * @return double containing time of the measurement corresponding to the index.
	 */
	public double time(int index) {
		GPS_Measurement meas = this.get(index);
		return meas.t();
	}

	/**
	 * Returns the predicted measurement based on the current state
	 * @param index measurement index
	 * @param t time of the measurement
	 * @param xref VectorN with the current state at the measurement time
	 */
	public boolean hasNext(int index) {
		boolean out = false;
		if (index < (this.size())) {
			out = true;
		}
		return out;
	}


	/** Read the measurement data from a tab-delimited ASCII text file.
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
				if (total != 5) {
					System.out.println(
						"GPS_MeasurementList.readFromFile: Number of columns do not match");
					System.exit(-99);
				}

				double[] temp = new double[4];
				for (int i = 0; i < 4; i++) {
					String token = tok.nextToken();
					temp[i] = Double.parseDouble(token);
				}
				String toke = tok.nextToken();
				int sv = Integer.parseInt(toke);
				GPS_Measurement meas =
					new GPS_Measurement(temp[0], temp[1], temp[2], temp[3], sv);
				this.add(meas);
			}
			in.close();
			fr.close();
		} catch (IOException e) {
			System.err.println("Error opening:" + file);
			return;
		}
	}

	/** Write the measurement data out to tab-delimited ASCII text file.
	 * @param file filename and directory
	 */
	public void sendToFile(String file) {
		try {
			FileOutputStream outfile = new FileOutputStream(file);
			PrintWriter pw = new PrintWriter(outfile);
			int index = 0;

			// loop through the file, one line at a time
			while (this.hasNext(index)) {
				GPS_Measurement meas = this.get(index);
				pw.println(
					meas.t()
						+ "\t"
						+ meas.t_mjd()
						+ "\t"
						+ meas.range()
						+ "\t"
						+ meas.phaseRange()
						+ "\t"
						+ meas.svid());
				index = index + 1;
			}
			pw.close();
			outfile.close();
		} catch (IOException e) {
			System.err.println("Error opening:" + file);
			return;
		}
	}

	/** Recover a serialized RGPS_MeasurementList File
	 * @param filename string containing the directory and filename.
	 * @return the trajectory
	 */
	public static GPS_MeasurementList recover(String filename) {
		GPS_MeasurementList out = new GPS_MeasurementList();
		try {

			FileInputStream file = new FileInputStream(filename);
			ObjectInputStream in = new ObjectInputStream(file);
			out = (GPS_MeasurementList) in.readObject();
			in.close();
			file.close();
		} catch (Exception e) {
			System.err.println("recover: " + e);
		}
		return out;
	}

	/** Write the trajectory out to a GPS_Measurement File
	 * @param filename string containing the directory and filename.
	 */
	public void serialize(String filename) {
		try {
			FileOutputStream file = new FileOutputStream(filename);
			ObjectOutputStream out = new ObjectOutputStream(file);
			out.writeObject(this);
			out.close();
			file.close();
		} catch (Exception e) {
			System.err.println("serialize: " + e);
		}
	}
	
	public static void main(String[] args){
		String dir = "C:\\Jat\\jat\\input\\";
		String gpsmeasfile = "gpsmeas.jat";
		String outfile = "gpsmeasblk.txt";
		GPS_MeasurementList gps = GPS_MeasurementList.recover(dir+gpsmeasfile);
		gps.sendToFile(dir+outfile);
		
	}
}
