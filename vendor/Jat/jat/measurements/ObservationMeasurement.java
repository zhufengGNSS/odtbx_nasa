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
 * File Created on 3 October 2005
 * Created by Kathryn Bradley, Emergent Space Technologies
 * */
package jat.measurements;

import jat.alg.estimators.MeasurementFileModel;
import jat.alg.estimators.MeasurementModel;
import jat.matvec.data.RotationMatrix;
import jat.matvec.data.VectorN;
import jat.spacecraft.SpacecraftModel;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Vector;

public class ObservationMeasurement implements Serializable{
	
	//** Flags for measurement data type
	public static final int DATA_UNKNOWN = 0;
	public static final int DATA_PSEUDORANGE = 1;
	public static final int DATA_DOPLER = 2;
	public static final int DATA_UTILITY = 3;
	public static final int DATA_CARRIERPHASE = 4;
	public static final int DATA_STATE_3 = 5;
	public static final int DATA_STATE_6 = 6;
	private int[] data_index = new int[7];

	//** Flags for observation type
	public static final int TYPE_GPS = 1;
	public static final int TYPE_GPSSTATE = 2;
	public static final int TYPE_RANGE = 3;
	public static final int TYPE_STATE = 4;
	public static final int TYPE_STATEUPDATE = 5;
	
	/** Epoch time of measurement in MJD */
	private double time_mjd;
	
	/** List of Measurement Types from Header */
	private Vector measurementType = new Vector();
	
	/** Observation Type Flag */
	private int obs_type;
	/** Observation Flag indicating which state is being observed **/
	private int whichState;
	
	/** Observation Measurements from File */
	private Vector measurement = new Vector();	
	private VectorN obs_state;
	private double pseudorange;
	private boolean data_initialized=false;
	
	/** Satellite PRN number */
	public Object prn;
	
	private RotationMatrix ECF2ECI;
	
	private int stateNum = 0;  //????? from EKF.stateNum
	
	/** Measurement Model */
	private MeasurementFileModel measModel;
	
	/** Constructor
	 * @param timeMjd Time of the measurement in MJD
	 * @param observationType an array with the types of measurements
	 * @param measurement actual observation measurements.
	 * @param  PRN satellite ID number
	 */
	public ObservationMeasurement(double timeMjd, Vector observationType, 
			Vector measurement, Object prn, HashMap input, MeasurementFileModel mfm,int type) {
		this.time_mjd = timeMjd;
		this.measurementType = observationType;
		this.measurement = measurement;
		this.prn = prn;
		this.obs_type = type; //chooseType((String)prn);
		this.measModel = mfm;
		initializeMeasType();
	}
	
//	public ObservationMeasurement(double timeMjd, double x,double y, double z, 
//			String id, HashMap input, MeasurementFileModel mfm){
//		this.time_mjd = timeMjd;
//		this.measurementType = new Vector();
//		this.obs_state = new VectorN(x,y,z);
//		this.measurement = new Vector();
//		measurement.add(new Double(x));
//		measurement.add(new Double(y));
//		measurement.add(new Double(z));
//		this.prn = "0";
//		this.obs_type = ObservationMeasurement.TYPE_GPSSTATE;
//		this.measModel = mfm;
//		data_index[DATA_STATE_3] = 0;
//	}
	
	/** Override toString()
	 * @return Observation data on a single line
	 */
	public String toString() {
		String out = time_mjd+"\t"+measurementType+"\t"+measurement+"\t"+prn;
		return out;
	}
	
	protected void set_ECF2ECI(RotationMatrix m){
		this.ECF2ECI = m;
	}
	public RotationMatrix get_ECF2ECI(){
		return this.ECF2ECI;
	}
	
	/**
	 * The ID of the spacecraft in the simulation to estimate.
	 * @return int ID number
	 */
	public int get_sc_id(){
		return get_PRN();
		//return this.spacecraft_id;
	}
	public int get_stateNum(){
		return this.stateNum;
	}
	public double get_mjd(){
		return this.time_mjd;
	}
	
	private void initializeMeasType(){
		try{
		for(int i=0; i<measurementType.size(); i++){
			if(measurementType.get(i).toString().equalsIgnoreCase("C1")){
				data_index[DATA_PSEUDORANGE]=i;
				pseudorange = Double.parseDouble((String)measurement.get(i));
			}else if(measurementType.get(i).toString().equalsIgnoreCase("D1")){
				data_index[DATA_DOPLER]=i;
			}else if(measurementType.get(i).toString().equalsIgnoreCase("P1")){
				data_index[DATA_CARRIERPHASE]=i;
			}else if(measurementType.get(i).toString().equalsIgnoreCase("S1")){			
				data_index[DATA_UTILITY]=i;
			}else if(measurementType.get(i).toString().equalsIgnoreCase("S3")){
				data_index[DATA_STATE_3]=i;
				obs_state = new VectorN(3);
				obs_state.x[0]= Double.parseDouble((String)measurement.get(0));
				obs_state.x[1]= Double.parseDouble((String)measurement.get(1));
				obs_state.x[2]= Double.parseDouble((String)measurement.get(2));
			}else if(measurementType.get(i).toString().equalsIgnoreCase("S6")){
				data_index[DATA_STATE_6]=i;
				obs_state = new VectorN(6);
				obs_state.x[0]= Double.parseDouble((String)measurement.get(0));
				obs_state.x[1]= Double.parseDouble((String)measurement.get(1));
				obs_state.x[2]= Double.parseDouble((String)measurement.get(2));
				obs_state.x[3]= Double.parseDouble((String)measurement.get(3));
				obs_state.x[4]= Double.parseDouble((String)measurement.get(4));
				obs_state.x[5]= Double.parseDouble((String)measurement.get(5));
			}else{
				data_index[DATA_UNKNOWN]=i;
			}
		}
		this.data_initialized = true;
		}catch(NumberFormatException e){
			System.err.println("Error occured when initializing observation data");
		}
	}
	
	/**
	 * Return the measurement time
	 * @return the measurement time in MJD
	 */
	public double time_mjd(){
		return this.time_mjd;
	}
	
	/**
	 * Return the measurement type
	 * @return double array containing the measurement types
	 */
	public Vector measurementType() {
		return this.measurementType;
	}
	
	/**
	 * Return the carrier phase range measurement
	 * @return double containing the carrier phase range measurement
	 */
	public Vector measurement(){
		return this.measurement;
	}
	
	/**
	 * Return the SV index
	 * @return int containing the SV index
	 */
	public Object prn() {
		return this.prn;
	}
	
	public static int chooseType(String s){
//		TODO Enumerate the different types
		if(s.contains("c")){
			return ObservationMeasurement.TYPE_RANGE;
		}else if(s.contains("u")){
			return ObservationMeasurement.TYPE_STATEUPDATE;
		}else{
			return ObservationMeasurement.TYPE_GPS;
		}
	}
	
	/**
	 * Returns the measurement data
	 * @param type the measurement type - use static values DATA_MEASTYPE
	 * @return the value of the measurement
	 */
	public double get_obs_data(int type){
		int i=0;
		return Double.parseDouble((String)measurement.get(data_index[type]));
	}
	
	public double get_range(){
		if(this.data_initialized){
			return this.pseudorange;
		}else{
			System.err.println("Error when retrieving range - not initialized");
			System.out.println("Error when retrieving range - not initialized");
			System.exit(0);
			return 0;
		}
	}
	
	public VectorN get_state(int size){
		if(this.data_initialized){
			if(this.obs_state.x.length==size)
				return this.obs_state;
			else{
				System.err.println("Error when retrieving state - size_data "
						+this.obs_state.x.length+"size_in "+size);
				System.out.println("Error when retrieving state - size_data "
						+this.obs_state.x.length+"size_in "+size);
				System.exit(0);
				return null;
			}
		}else{
			System.err.println("Error when retrieving state - not initialized");
			System.out.println("Error when retrieving state - not initialized");
			System.exit(0);
			return null;
		}
	}
	
	/**
	 * Returns the residual evaluated through the spacecraft computer
	 * @param state Current spacecraft state (to be adjusted)
	 * @return The observation (observed minus predicted)
	 */
	public double get_residual(VectorN state) {
		return this.measModel.zPred(this,get_PRN(),time_mjd,state);
	}
	
	public double get_noise() {
		return this.measModel.R(this);
	}
	
	public double get_noise(SpacecraftModel sc){
		switch(obs_type){
		case TYPE_GPS:
			return sc.get_GPS_noise(0,true);
		case TYPE_GPSSTATE:
			return sc.get_GPS_noise(whichState,false);
		case TYPE_RANGE:
			return sc.get_GPS_noise(0,true);
		case TYPE_STATE:
			return sc.get_GPS_noise(whichState,false);
		case TYPE_STATEUPDATE:
			return sc.get_GPS_noise(whichState,false);
		default:
			return 0;
		}
	}
	
	public String get_measurementType() {
		switch(obs_type){
		case TYPE_GPS:
			return "GPS";
		case TYPE_GPSSTATE:
			return "GPS_State";
		case TYPE_RANGE:
			return "Range";
		case TYPE_STATE:
			return "State";
		case TYPE_STATEUPDATE:
			return "State_Update";
		default:
			return null;
		}
	}
	
	public int get_type(){
		return this.obs_type;
	}
	
	public VectorN get_H(VectorN xref) {
		return this.measModel.H(this,xref);
	}

	public int get_PRN() {
		//System.out.println((String)prn);
		String id = (String)prn;
		int out;
		try{
			out = Integer.parseInt((String)prn);
		}catch(NumberFormatException e){
			char[] array = id.toCharArray();
			int i=0;
			while(array[i]==' ' || array[i]=='c' || array[i]=='G') i++; 
			id = new String(array,i,array.length-i);
			out= Integer.parseInt(id);
		}
		return out;
	}
	
	public int get_whichState(){
		return whichState;
	}
	public void set_whichState(int arg){
		whichState = arg;
	}
}

