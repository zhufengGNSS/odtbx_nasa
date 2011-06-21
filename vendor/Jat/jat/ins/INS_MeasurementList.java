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

package jat.ins;
//import jat.alg.estimators.*;
//import jat.ins.*;
import jat.matvec.data.*;
import jat.alg.integrators.*;

import java.io.*;
import java.util.*;

/**
* The INS_MeasurementList.java Class provides a way to deal with
* a list of INS measurements read from a file.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class INS_MeasurementList implements Serializable {

	private ArrayList list = new ArrayList();
	
	/** Constructor
	 */
	public INS_MeasurementList() {
	}

	/** Constructor
	 * @param filename String containing directory and filename where the data resides
	 */
	public INS_MeasurementList(String filename) {
		this.readFromFile(filename);
	}

	/** Add a INS Measurement to the collection
	 * @param meas INS_Measurement object
	 */
	public void add(INS_Measurement meas) {
		list.add(meas);
	}

	/** Get a GPS_Measurement out of the collection
	 * @param index Index of the measurement
	 * @return the GPS Measurement
	 */
	public INS_Measurement get(int index) {
		return (INS_Measurement) list.get(index);
	}

	/** Return the size of the list
	 * @return the number of measurements in the list
	 */
	public int size() {
		return list.size();
	}

	/** Returns the accelerometer measurement
	 * @param index int containing the measurement index
	 * @return double containing the measurement corresponding to the index.
	 */
	public VectorN f(int index) {
		INS_Measurement meas = this.get(index);
		return meas.f;
	}

	/**
	 * Returns the time of the measurement
	 * @param index int containing the measurement index
	 * @return double containing time of the measurement corresponding to the index.
	 */
	public double time(int index) {
		INS_Measurement meas = this.get(index);
		return meas.t;
	}

	/**
	 * Returns whether there is more data
	 * @param index measurement index
	 * @return true if there is more data to be read
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
				if (total != 7) {
					System.out.println(
						"INS_MeasurementList.readFromFile: Number of columns do not match");
					System.exit(-99);
				}

				double[] temp = new double[7];
				for (int i = 0; i < 7; i++) {
					String token = tok.nextToken();
					temp[i] = Double.parseDouble(token);
				}
				VectorN sf = new VectorN(temp[1], temp[2], temp[3]);
				VectorN w = new VectorN(temp[4], temp[5], temp[6]);
				INS_Measurement meas = new INS_Measurement(temp[0], sf, w);
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
		LinePrinter lp = new LinePrinter(file);
		int index = 0;

		// loop through the file, one line at a time
		while (this.hasNext(index)) {
			INS_Measurement meas = this.get(index);
			VectorN sf = meas.f;
			VectorN w = meas.omega;
			VectorN out = new VectorN(sf, w);
			lp.print(meas.t, out.x);
			index = index + 1;
		}
		lp.close();
	}

	/** Recover a serialized INS_MeasurementList File
	 * @param filename string containing the directory and filename.
	 * @return the trajectory
	 */
	public static INS_MeasurementList recover(String filename) {
		INS_MeasurementList out = new INS_MeasurementList();
		try {

			FileInputStream file = new FileInputStream(filename);
			ObjectInputStream in = new ObjectInputStream(file);
			out = (INS_MeasurementList) in.readObject();
			in.close();
			file.close();
		} catch (Exception e) {
			System.err.println("recover: " + e);
		}
		return out;
	}

	/** Write the trajectory out to a INS_Measurement File
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

}
