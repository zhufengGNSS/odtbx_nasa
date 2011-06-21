package jat.cm.rendezvous;

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
 * File Created on Aug 26, 2003
 */
 
import jat.matvec.data.*;
import jat.alg.integrators.*;

import java.io.*;
import java.util.*;

/**
 * <P>
 * The CW_Target Class contains a list of guidance targets.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */ 
public class CW_TargetList {
	private ArrayList list = new ArrayList();
	
	/** Constructor
	 */
	public CW_TargetList() {
	}

	/** Constructor
	 * @param filename String containing directory and filename where the data resides
	 */
	public CW_TargetList(String filename) {
		this.readFromFile(filename);
	}

	/** Add a DeltaV to the collection
	 * @param dv DeltaV object
	 */
	public void add(CW_Target tgt) {
		list.add(tgt);
	}

	/** Get a GPS_Measurement out of the collection
	 * @param index Index of the measurement
	 * @return the GPS Measurement
	 */
	public CW_Target get(int index) {
		return (CW_Target) list.get(index);
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
	public VectorN rtgt(int index) {
		CW_Target tgt = this.get(index);
		return tgt.rtgt;
	}

	/**
	 * Returns the time of the measurement
	 * @param index int containing the measurement index
	 * @return double containing time of the measurement corresponding to the index.
	 */
	public double time(int index) {
		CW_Target tgt = this.get(index);
		return tgt.t;
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
				if (total != 4) {
					System.out.println(
						"CW_TargetList.readFromFile: Number of columns do not match");
					System.exit(-99);
				}

				double[] temp = new double[4];
				for (int i = 0; i < 4; i++) {
					String token = tok.nextToken();
					temp[i] = Double.parseDouble(token);
				}
				VectorN rtgt = new VectorN(temp[1], temp[2], temp[3]);
				CW_Target tgt = new CW_Target(temp[0], rtgt);
				this.add(tgt);
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
			CW_Target tgt = this.get(index);
			VectorN rtgt = tgt.rtgt;
			lp.print(tgt.t, rtgt.x);
			index = index + 1;
		}
		lp.close();
	}

	/** Recover a serialized DeltaV_List File
	 * @param filename string containing the directory and filename.
	 * @return the trajectory
	 */
	public static CW_TargetList recover(String filename) {
		CW_TargetList out = new CW_TargetList();
		try {

			FileInputStream file = new FileInputStream(filename);
			ObjectInputStream in = new ObjectInputStream(file);
			out = (CW_TargetList) in.readObject();
			in.close();
			file.close();
		} catch (Exception e) {
			System.err.println("recover: " + e);
		}
		return out;
	}

	/** Write the trajectory out to a DeltaV File
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
