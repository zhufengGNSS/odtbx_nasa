package jat.traj;

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
 */
import java.io.*;
import java.util.*;

import jat.matvec.data.VectorN;
import jat.timeRef.*;
import jat.alg.integrators.*;
/**
* <P>
* The Trajectory.java Class provides the means for storing and accessing trajectory data.
* Trajectory data includes at least time, position and velocity and may include other
* variables such as mass, thrust, specific force, etc.
* 
* Note: time can be in sim time or mjd, user needs to be consistent.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/

public class Trajectory implements Serializable, Printable {

	/** The trajectory data is stored in a TrajectoryList
	 */
	private TrajectoryList traj;

	/** The trajectory title */
	private String title;
	
	/** The column labels */
	private String[] labels;

	/** The number of columns of data */
	private int ncol;
	
	/** The trajectory epoch time (UTC) */
	private CalDate epoch;

	/** The trajectory central body */
	private CentralBody cb;

	/** The trajectory coordinate system */
	private CoordinateSystem cs;

	/** The trajectory distance units */
	private DistanceUnits du;

	/** The trajectory time units */
	private TimeUnits tu;

	/** False if the trajectory contains data */
	private boolean firstData = true;

	/** Construct a new Trajectory.
	 */
	public Trajectory() {
		traj = new TrajectoryList();
	}

	//*************METHODS TO READ OR ADD TRAJECTORY DATA	

	/** Add data to the trajectory. The first call to add sets the number of columns.
	 * @param in double array containing the next row of data to be added.
	 */
	public void add(double[] in) {
		if (this.firstData) {
			this.ncol = in.length;
			traj.add(in);
			this.firstData = false;
		} else {
			if (in.length != this.ncol) {
				System.out.println(
					"Trajectory.add: Number of columns do not match");
				System.exit(-99);
			} else {
				traj.add(in);
			}
		}
	}

	/** Add data to the trajectory. The first call to add sets the number of columns.
	 * @param t double containing time.
	 * @param in double array containing the next row of data to be added.
	 */
	public void add(double t, double[] in) {
		if (this.firstData) {
			this.ncol = in.length + 1;
			traj.add(t, in);
			this.firstData = false;
		} else {
			if (in.length != (this.ncol - 1)) {
				System.out.println(
					"Trajectory.add: Number of columns do not match");
				System.exit(-99);
			} else {
				traj.add(t, in);
			}
		}
	}

	/** Add data to the trajectory. The first call to add sets the number of columns.
	 * @param t double containing time.
	 * @param y1 double array containing part of data to be added.
	 * @param y2 double array containing the remainder of the data to be added.
	 */
	public void add(double t, double[] y1, double[] y2) {
		if (this.firstData) {
			this.ncol = y1.length + y2.length + 1;
			traj.add(t, y1, y2);
			this.firstData = false;
		} else {
			if ((y1.length + y2.length) != (this.ncol - 1)) {
				System.out.println("Trajectory.add: Number of columns do not match");
				System.exit(-99);
			} else {
				traj.add(t, y1, y2);
			}
		}
	}
	
    /** Implements the Printable interface to get the data out of the propagator and pass it to the trajectory.
     *  This method is executed by the propagator at each integration step.
     * @param t Time.
     * @param y Data array.
     */
    public void print(double t, double [] y)
    {
    	traj.add(t, y);
    }
	

	/** Recover a serialized JAT Trajectory File
	 * @param filename string containing the directory and filename.
	 * @return the trajectory
	 */
	public static Trajectory recover(String filename) {
		Trajectory out = new Trajectory();
		try {

			FileInputStream file = new FileInputStream(filename);
			ObjectInputStream in = new ObjectInputStream(file);
			out = (Trajectory) in.readObject();
			in.close();
			file.close();
		} catch (Exception e) {
			System.err.println("recover: " + e);
		}
		return out;
	}

	/** Read the trajectory data from a tab-delimited ASCII text file.
	 * The first line read in sets the number of columns.
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
				if (this.firstData) {
					this.ncol = total;
					this.firstData = false;
				} else {
					if (total != this.ncol) {
						System.out.println(
							"Trajectory.readFromFile: Number of columns do not match");
						System.exit(-99);
					}
				}

				double[] temp = new double[total];
				int i = 0;
				while (tok.hasMoreTokens()) {
					String token = tok.nextToken();
					//					System.out.println("token= "+token);
					temp[i] = Double.parseDouble(token);
					i = i + 1;
				}
				traj.add(temp);
			}
			in.close();
			fr.close();
		} catch (IOException e) {
			System.err.println("Error opening:" + file);
			return;
		}
	}
	/** Read the trajectory data from a tab-delimited ASCII text file.
	 * The first line read in sets the number of columns.
	 * @param file filename and directory
	 * @param delim the delimiter character
	 */
	public void readFromFile(String file, String delim, String units) {
		try {
			FileReader fr = new FileReader(file);
			BufferedReader in = new BufferedReader(fr);
			String line;

			// loop through the file, one line at a time
			while ((line = in.readLine()) != null) {
				StringTokenizer tok = new StringTokenizer(line, delim);
				int total = tok.countTokens();

				// check for consistent number of columns
				if (this.firstData) {
					this.ncol = total;
					this.firstData = false;
				} else {
					if (total != this.ncol) {
						System.out.println(
							"Trajectory.readFromFile: Number of columns do not match");
						System.exit(-99);
					}
				}

				double[] temp = new double[total];
				int i = 0;
				while (tok.hasMoreTokens()) {
					String token = tok.nextToken();
					//					System.out.println("token= "+token);
					temp[i] = Double.parseDouble(token);
					if(units.equalsIgnoreCase("km") && i>0){
						temp[i] = temp[i]*1000;
					}
					i = i + 1;
				}
				traj.add(temp);
			}
			in.close();
			fr.close();
		} catch (IOException e) {
			System.err.println("Error opening:" + file);
			return;
		}
	}

	//*************METHODS TO WRITE OUT OR STORE TRAJECTORY DATA

	/** Write the trajectory out to a JAT Trajectory File
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

	/** Write the trajectory out using a LinePrinter
	 * @param LinePrinter string containing the directory
	 */
	public void sendToLinePrinter(LinePrinter lp) {
		for (int i = 0; i < this.npts(); i++) {
			double[] temp = traj.get(i);
			lp.print(temp);
		}
		lp.close();
	}

	/** Write the trajectory data to a tab-delimited ASCII text file.
	 * The first line writen sets the number of columns.
	 * @param file filename and directory
	 */
	public void writeToFile(String file) {
			LinePrinter out = new LinePrinter(file);
			String line;

			// loop through the file, one line at a time
			double[] r = new double[7];
			for(int i=0; i<traj.size(); i++){
			    r = traj.get(i);
			    out.println(""+r[0]+"\t"+r[1]+"\t"+r[2]+"\t"+r[3]+"\t"+r[4]+"\t"+r[5]+"\t"+r[6]);
			}
			out.close();
			System.err.println("Error opening:" + file);
			return;
	}

	//*************METHODS TO ACCESS (GET) TRAJECTORY DATA

	/** Return the number of trajectory points
	 * @return the number of trajectory points
	 */
	public int npts() {
		return this.traj.size();
	}

	/** Returns the next row of trajectory data
	 * @return double array containing the next row of trajectory data 
	 */
	public double[] next() {
		double[] out = this.traj.next();
		return out;
	}
	
	/** Returns the row of trajectory data at the given index
	 * @return double array containing the row of trajectory data
	 */
	public double[] get(int i){
	    return this.traj.get(i);
	}
	/**
	 * Returns the state vector at the given index.
	 * Identical to get(int i) but excluding the first collumn
	 * @param i index
	 * @return double array containing the row of data skipping the first entry
	 */
	public double[] getState(int i){
		double[] out = new double[ncol-1];
		double[] data = this.traj.get(i);
		for(int j=1; j<ncol; j++){
			out[j-1] = data[j];
		}
		return out;
	}
	
	/** Determines if there is any more trajectory data left
	 * @return true if there is more trajectory data left
	 */
	public boolean hasNext() {
		boolean out = this.traj.hasNext();
		return out;
	}
	
	public double[] getNext() {
		return this.traj.next();
	}

	/** Return a double array containing the trajectory data in Java3D format,
	 *  which is: [x1, y1, z1, x2, y2, z2,...]
	 * @return the trajectory data in Java3D array format
	 */
	public double[] j3dArray() {
		int n = this.npts();
		double[] out = new double[3 * n];
		for (int i = 0; i < n; i++) {
			double[] temp = traj.get(i);
			out[3 * i] = temp[1];
			out[3 * i + 1] = temp[2];
			out[3 * i + 2] = temp[3];
		}
		return out;
	}
	/** Return a double array containing only times.
	 * @return double array containing times.
	 */
	public double[] timeArray() {
		int n = this.npts();
		double[] out = new double[n];
		for (int i = 0; i < n; i++) {
			double[] temp = traj.get(i);
			out[i] = temp[0];
		}
		return out;
	}
	
	public double[][] positionArray(){
	    int n = this.npts();
	    double[][] out = new double[n][3];
	    for (int i=0; i<n; i++){
	        double[] tmp = traj.get(i);
	        out[i][0] = tmp[1];
	        out[i][1] = tmp[2];
	        out[i][2] = tmp[3];
	    }
	    return out;
	}

	//*************METHODS TO GET TRAJECTORY ATTRIBUTES

	/** Return the time at the given index
	 * @return The time (first element of the row at index i)
	 */
	public double getTimeAt(int i){
		double out = 0.0;
		if ((i >= 0)&&(i < this.size())) {
			double[] tmp = traj.get(i);
			out = tmp[0];
		} else {
			System.err.println("Trajectory.getTimeAt: i out of bounds");
		}
		return out;
	}
	
	/** Return the title of the trajectory
	 * @return String containing the title of the trajectory
	 */
	public String getTitle() {
		return this.title;
	}

	/** Return the column labels of the trajectory
	 * @return String containing the column labels of the trajectory
	 */
	public String[] getColumnLabels() {
		return this.labels;
	}

	/** Return the number of columns in the trajectory
	 * @return int containing the number of columns in the trajectory
	 */
	public int numberOfColumns() {
		return this.ncol;
	}

	/** Return the epoch time of the trajectory in MJD
	 * @return double containing the epoch time of the trajectory in MJD
	 */
	public double epochMJD() {
		return this.epoch.mjd();
	}

	/** Return the epoch time of the trajectory
	 * @return String containing the epoch time of the trajectory
	 */
	public String getEpoch() {
		return this.epoch.toString();
	}

	/** Return the CentralBody for the trajectory
	 * @return the CentralBody for the trajectory
	 */
	public CentralBody centralBody() {
		return this.cb;
	}

	/** Return the CoordinateSystem for the trajectory
	 * @return the CoordinateSystem for the trajectory
	 */
	public CoordinateSystem coordinateSystem() {
		return this.cs;
	}

	/** Return the DistanceUnits used in the trajectory
	 * @return the DistanceUnits used in the trajectory
	 */
	public DistanceUnits distanceUnits() {
		return this.du;
	}

	/** Return the TimeUnits used in the trajectory
	 * @return the TimeUnits used in the trajectory
	 */
	public TimeUnits timeUnits() {
		return this.tu;
	}

	/** Print the trajectory including attributes to a LinePrinter
	 * @param LinePrinter to use for printing.
	 */
	public void printAll(LinePrinter lp) {
		lp.println("Title = " + this.title);
		lp.println("NumberOfRows = " + this.npts());
		lp.println("NumberOfColumns = " + this.ncol);
		lp.println("ScenarioEpoch = " + this.epoch);
		lp.println("CentralBody = " + this.cb);
		lp.println("CoordinateSystem = " + this.cs);
		lp.println("DistanceUnits = " + this.du);
		lp.println("TimeUnits = " + this.tu);
		if(this.labels==null){
			String[] tmp = {"t ","x [m]","y [m]","z [m]","xdot [m]","ydot [m]","zdot [m]"};
			labels = tmp;
		}
		lp.println(this.labels);
		this.sendToLinePrinter(lp);
	}

	//*************METHODS TO SET TRAJECTORY ATTRIBUTES
	/** Resets the index of the trajectory list.
	 */
	public void reset(){
		this.traj.reset();
	}
	/** Set the trajectory title
	 * @param t String containing the title
	 */
	public void setTitle(String t) {
		this.title = t;
	}

	/** Set the column labels
	 * @param l String array containing the labels
	 */
	public void setLabels(String[] l) {
		this.labels = new String[l.length];
		System.arraycopy(l, 0, this.labels, 0, l.length);
	}

	/** Set the epoch time
	 * @param Yr int containing the year, 4 digits
	 * @param Mon int containing the month
	 * @param D int containing the day of the month
	 * @param Hr int containing the hour of the day
	 * @param Mn int containing the minutes
	 * @param S double containing the seconds
	 */
	public void setEpoch(int Yr, int Mon, int D, int Hr, int Mn, double S) {
		this.epoch = new CalDate(Yr, Mon, D, Hr, Mn, S);
	}

	/** Set the epoch time
	 * @param mjd double containing the MJD time of epoch
	 */
	public void setEpoch(double mjd) {
		this.epoch = new CalDate(mjd);
	}

	/** Set the Central Body
	 * @param m a CentralBody
	 */
	public void setCentralBody(CentralBody m) {
		this.cb = m;
	}

	/** Set the Coordinate System
	 * @param c a CoordinateSystem
	 */
	public void setCoordinateSystem(CoordinateSystem c) {
		this.cs = c;
	}

	/** Set the Distance Units
	 * @param d the DistanceUnits
	 */
	public void setDistanceUnits(DistanceUnits d) {
		this.du = d;
	}

	/** Set the Time Units
	 * @param t the TimeUnits
	 */
	public void setTimeUnits(TimeUnits t) {
		this.tu = t;
	}
	
	/** Return the size (number of elements) of the trajectory.
	 * @return size of the TrajectoryList
	 */
	public int size(){
		return this.traj.size();
	}

	/**
	 * Returns the position at the given modified julian date or sim time
	 * @param time 
	 * @return VectorN position 
	 */
	public VectorN getPositionAt(double time) {
		int i=0;
		while(i<this.size() && this.getTimeAt(i)<=(time+1e-5)){
			if(Math.abs(this.getTimeAt(i)-time)<1e-5){
				double[] data = this.traj.get(i);
				return new VectorN(data[1],data[2],data[3]);
			}
			i++;
		}
		System.err.println("Trajectory.getStateAt: unable to find position data at: "+time+"  nearest: "+this.getTimeAt(i));
		return new VectorN(3);
	}
	/**
	 * Returns the full state at the given modified julian date or sim time
	 * @param time
	 * @return VectorN state
	 */
	public VectorN getStateAt(double time){
		int i=0;
		while(i<this.size() && this.getTimeAt(i)<=time){
			if(Math.abs(this.getTimeAt(i)-time)<1e-5){
				double[] data = this.traj.get(i);
				double[] out = new double[ncol-1];
				for(int j=1; j<this.ncol; j++) out[j-1] = data[j];
				return new VectorN(out);
			}
			i++;
		}
		System.err.println("Trajectory.getStateAt: unable to find position data at: "+time+"  nearest: "+this.getTimeAt(i));
		return new VectorN(this.ncol-1);
	}
}
