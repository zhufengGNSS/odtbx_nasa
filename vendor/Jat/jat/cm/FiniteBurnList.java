
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
 * File Created on Aug 28, 2003
 */
 
package jat.cm;
 
import jat.matvec.data.*;
import jat.alg.integrators.*;

import java.io.*;
import java.util.*;

/**
* The GroundStationList.java Class provides a way to deal with
* a list of finite burns read from a file.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class FiniteBurnList implements Serializable {

	private ArrayList list = new ArrayList();
	
	/** Constructor
	 */
	public FiniteBurnList() {
	}

	/** Constructor
	 * @param filename String containing directory and filename where the data resides
	 */
	public FiniteBurnList(String filename) {
		this.readFromFile(filename);
	}

	/** Add a burn to the collection
	 * @param burn FiniteBurn object
	 */
	public void add(FiniteBurn burn) {
		list.add(burn);
	}

	/** Get a FiniteBurn out of the collection
	 * @param index Index of the measurement
	 * @return the FiniteBurn
	 */
	public FiniteBurn get(int index) {
		return (FiniteBurn) list.get(index);
	}

	/** Return the size of the list
	 * @return the number of burns in the list
	 */
	public int size() {
		return list.size();
	}

	/** Returns the thrust direction unit vector
	 * @param index int containing the burn index
	 * @return VectorN containing the thrust direction unit vector.
	 */
	public VectorN unitVector(int index) {
		FiniteBurn burn = this.get(index);
		return burn.unitVector;
	}

	/**
	 * Returns the start time of the burn
	 * @param index int containing the burn index
	 * @return double containing the start time.
	 */
	public double startTime(int index) {
		FiniteBurn burn = this.get(index);
		return burn.tstart;
	}

	/**
	 * Returns the stop time of the burn
	 * @param index int containing the burn index
	 * @return double containing the stop time.
	 */
	public double stopTime(int index) {
		FiniteBurn burn = this.get(index);
		return burn.tstop;
	}

	/**
	 * Returns the acceleration of the burn
	 * @param index int containing the burn index
	 * @return double containing time of the acceleration.
	 */
	public double accel(int index) {
		FiniteBurn burn = this.get(index);
		return burn.accel;
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

	/** Read the burn data from a tab-delimited ASCII text file.
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
						"DeltaV_List.readFromFile: Number of columns do not match");
					System.exit(-99);
				}

				double[] temp = new double[6];
				for (int i = 0; i < 6; i++) {
					String token = tok.nextToken();
					temp[i] = Double.parseDouble(token);
				}
				VectorN unit = new VectorN(temp[3], temp[4], temp[5]);
				FiniteBurn burn = new FiniteBurn(temp[0], temp[1], temp[2], unit);
				this.add(burn);
			}
			in.close();
			fr.close();
		} catch (IOException e) {
			System.err.println("Error opening:" + file);
			return;
		}
	}
	/** Write the burn data out to tab-delimited ASCII text file.
	 * @param file filename and directory
	 */
	public void sendToFile(String file) {
		LinePrinter lp = new LinePrinter(file);
		int index = 0;

		// loop through the file, one line at a time
		while (this.hasNext(index)) {
			FiniteBurn burn = this.get(index);
			String out = burn.tstart+"\t"+burn.tstop+"\t"+burn.accel+"\t"+burn.unitVector.toString();
			lp.println(out);
			index = index + 1;
		}
		lp.close();
	}

	/** Recover a serialized GroundStationList File
	 * @param filename string containing the directory and filename.
	 * @return the trajectory
	 */
	public static FiniteBurnList recover(String filename) {
		FiniteBurnList out = new FiniteBurnList();
		try {

			FileInputStream file = new FileInputStream(filename);
			ObjectInputStream in = new ObjectInputStream(file);
			out = (FiniteBurnList) in.readObject();
			in.close();
			file.close();
		} catch (Exception e) {
			System.err.println("recover: " + e);
		}
		return out;
	}

	/** Write the burn data out to a FiniteBurn File
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
	
    public static void main(String[] args)
    {
        FiniteBurnList x = new FiniteBurnList();
        x.readFromFile("C:\\Jat\\jat\\input\\burns\\vbar_burns.txt");

		
        
        
        // Serialize the burns
        x.serialize("C:\\Jat\\jat\\input\\burns\\vbar_burns.jat");
        System.out.println("burn list serialized");
        
        // Recover the trajectory and print all to screen       
        FiniteBurnList bl = FiniteBurnList.recover("C:\\Jat\\jat\\input\\burns\\vbar_burns.jat");
        System.out.println("Printing Recovered DeltaV's");
        int index = 0;
        for (int i = 0; i < bl.size(); i++){
	        double tstart = bl.startTime(i);
	        double tstop = bl.stopTime(i);
	        double acc = bl.accel(i);
	        VectorN unit = bl.unitVector(i);
	        System.out.println("burn: "+i+" "+tstart+" "+tstop+" "+acc+" "+unit);
        }
    }	


}
