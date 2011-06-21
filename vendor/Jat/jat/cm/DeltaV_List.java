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

package jat.cm;
import jat.matvec.data.*;
import jat.alg.integrators.*;

import java.io.*;
import java.util.*;

/**
* The DeltaV_List.java Class provides a way to deal with
* a list of impulsive delta_V's read from a file.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class DeltaV_List implements Serializable {

	private ArrayList list = new ArrayList();
	
	/** Constructor
	 */
	public DeltaV_List() {
	}

	/** Constructor
	 * @param filename String containing directory and filename where the data resides
	 */
	public DeltaV_List(String filename) {
		this.readFromFile(filename);
	}

	/** Add a DeltaV to the collection
	 * @param dv DeltaV object
	 */
	public void add(DeltaV dv) {
		list.add(dv);
	}

	/** Get a GPS_Measurement out of the collection
	 * @param index Index of the measurement
	 * @return the GPS Measurement
	 */
	public DeltaV get(int index) {
		return (DeltaV) list.get(index);
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
	public VectorN dv(int index) {
		DeltaV deltav = this.get(index);
		return deltav.dv;
	}

	/**
	 * Returns the time of the measurement
	 * @param index int containing the measurement index
	 * @return double containing time of the measurement corresponding to the index.
	 */
	public double time(int index) {
		DeltaV deltav = this.get(index);
		return deltav.t;
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
						"DeltaV_List.readFromFile: Number of columns do not match");
					System.exit(-99);
				}

				double[] temp = new double[4];
				for (int i = 0; i < 4; i++) {
					String token = tok.nextToken();
					temp[i] = Double.parseDouble(token);
				}
				VectorN dv = new VectorN(temp[1], temp[2], temp[3]);
				DeltaV deltav = new DeltaV(temp[0], dv);
				this.add(deltav);
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
			DeltaV deltav = this.get(index);
			VectorN dv = deltav.dv;
			lp.print(deltav.t, dv.x);
			index = index + 1;
		}
		lp.close();
	}

	/** Recover a serialized DeltaV_List File
	 * @param filename string containing the directory and filename.
	 * @return the trajectory
	 */
	public static DeltaV_List recover(String filename) {
		DeltaV_List out = new DeltaV_List();
		try {

			FileInputStream file = new FileInputStream(filename);
			ObjectInputStream in = new ObjectInputStream(file);
			out = (DeltaV_List) in.readObject();
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
	
    public static void main(String[] args)
    {
        DeltaV_List x = new DeltaV_List();

		// delta-v from CW equations
		double [] deltav1 = new double[3];
		deltav1[0] = -0.19094164920359352;
//		deltav1[1] = -0.8776834149116475;
		deltav1[1] = -0.878;
		deltav1[2] = 0.0;
		VectorN dv1 = new VectorN(deltav1);
		
		DeltaV deltav = new DeltaV(0.0, dv1);
		
		x.add(deltav);
		
		double [] deltav2 = new double[3];
		deltav2[0] = -0.5603757896639145;
		deltav2[1] = 0.6249070361006656;
		deltav2[2] = 0.0;
		VectorN dv2 = new VectorN(deltav2);
		
		deltav = new DeltaV(5000.0, dv2);
		
		x.add(deltav);

		double [] deltav3 = new double[3];
		deltav3[0] = 0.1606911732628235;
		deltav3[1] = 0.14140915221166075;
		deltav3[2] = 0.0;
		VectorN dv3 = new VectorN(deltav3);
		
		deltav = new DeltaV(5750.0, dv3);
		
		x.add(deltav);

		double [] deltav4 = new double[3];
		deltav4[0] = 0.0707966875285207;
		deltav4[1] = 0.06230149024070955;
		deltav4[2] = 0.0;
		VectorN dv4 = new VectorN(deltav4);
		
		deltav = new DeltaV(6500.0, dv4);
		
		x.add(deltav);

		double [] deltav5 = new double[3];
		deltav5[0] = 0.03560401104777747;
		deltav5[1] = 0.034059802166368884;
		deltav5[2] = 0.0;
		VectorN dv5 = new VectorN(deltav5);
		
		deltav = new DeltaV(7250.0, dv5);
		
		x.add(deltav);
		

        // Print out the trajectory to the screen
//        LinePrinter lp1 = new LinePrinter();        
//        x.traj.printAll(lp1);
//        lp1.close();
        
        
        // Serialize the trajectory
        x.serialize("C:\\Jat\\jat\\input\\deltav.jat");
        System.out.println("deltav list serialized");
        
        // Recover the trajectory and print all to screen       
        DeltaV_List dvl = DeltaV_List.recover("C:\\Jat\\jat\\input\\deltav.jat");
        System.out.println("Printing Recovered DeltaV's");
        int index = 0;
        while (dvl.hasNext(index)){
	        double t = dvl.time(index);
	        VectorN delv = dvl.dv(index);
	        System.out.println("deltav: "+index+" "+t+" "+delv);
	        index = index + 1;
        }
    }
        
	

}
