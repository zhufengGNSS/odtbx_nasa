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
 */

package jat.measurements;


import java.util.HashMap;

import jat.sim.*;
import jat.util.FileUtil;
import jat.gps.Visible;
import jat.matvec.data.*;
import java.util.Random;
import jat.alg.estimators.*;

public class GPSStateMeasurementModel implements MeasurementFileModel,MeasurementModel{
	
	public static VectorN R;
	public static int numStates;
	HashMap hm;// = closedLoopSim.hm;  
	//* *NOTE* Added argument to default constructor instead of call to static var at instantiation
	Random generator;
	private int FILTER_states;
	
	public GPSStateMeasurementModel(HashMap h) {
		/*Add a sleep in here to insure that the Random Number
		 * Seeds don't allign with any other random number generator
		 */
		hm = h;
		this.FILTER_states = initializer.parseInt(hm,"FILTER.states");
		try {
			Thread.sleep(20);
			generator = new Random();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			
		}
	}
	
	public VectorN getMeasurement()
	{
		String tmp;
		VectorN pointSolution = new VectorN(3);
		
		tmp = "MEAS."+EKF.measNum+".satellite";
		int sat = initializer.parseInt(hm,tmp);
		
		//double[] truth = closedLoopSim.truth[sat].sc.get_spacecraft().toStateVector();
		double[] truth = EstimatorSimModel.truth[sat].get_spacecraft().toStateVector();
		
		//Add in the measurement noise read out of the file
		for(int j = 0; j < 3; j++)
		{
			tmp = "MEAS."+EKF.measNum+".R."+j;
			double R = initializer.parseDouble(hm,tmp);
			
			/*Scale the error as Gaussian noise times the 
			 square of the measurement noise*/
			//* TODO .nextGaussian() is already zero mean
			pointSolution.x[j] =  truth[j] +(2*(generator.nextGaussian()-0.5))*R*R;
		}	
		return pointSolution;
	}
	
	public VectorN predictMeasurement(VectorN state){
		
			String tmp = "MEAS."+EKF.measNum+".satellite";
			int sat = initializer.parseInt(hm,tmp);
			
			VectorN range = new VectorN(3);
			range.set(0,state.get(0+(6*sat)));
			range.set(1,state.get(1+(6*sat)));
			range.set(2,state.get(2+(6*sat)));
			
			return range;
	}
	public VectorN predictMeasurement(int sat,VectorN state){
		
		//String tmp = "MEAS."+EKF.measNum+".satellite";
		//int sat = initializer.parseInt(hm,tmp);
		
		VectorN range = new VectorN(3);
		range.set(0,state.get(0+(6*sat)));
		range.set(1,state.get(1+(6*sat)));
		range.set(2,state.get(2+(6*sat)));
		
		return range;
}
	
	public double  zPred(int i, double time, VectorN state){
		String tmp = "MEAS."+EKF.measNum+".satellite";
		int sat = initializer.parseInt(hm,tmp);
		
		VectorN oMinusC;
		VectorN pred = predictMeasurement(state);
		VectorN obs = getMeasurement();
		oMinusC      = obs.minus(pred);
		
		//Ensure we are returning the correct state when there is more than
		//one satellite
		int j = i - 6*sat; 
		return oMinusC.get(j);
	}
	public double  zPred(ObservationMeasurement om,int i, double time, VectorN state){
		//String tmp = "MEAS."+om.get_sc_id()+".satellite";
		//int sat = initializer.parseInt(hm,tmp);
		int whichState = om.get_whichState();		//* Note i = GPS_PRN  om.get_sc_id() = Sim_SC_ID
		VectorN oMinusC;
		VectorN pred = predictMeasurement(om.get_sc_id(),state);
		VectorN obs = om.get_state(3);
		oMinusC      = obs.minus(pred);
		
		//Ensure we are returning the correct state when there is more than
		//one satellite
		//int j = i - 6*sat; 
		//return oMinusC.get(j);
		return oMinusC.get(whichState);
	}
	
	/** Return the measurement noise value for this measurement
	 * 
	 *   
	 */
	public double R()
	{
		int whichState = EKF.stateNum;
		int measNum    = EKF.measNum;
		
		
		//Ensure we are returning the correct state when there is more than
		//one satellite
		String tmp = "MEAS."+EKF.measNum+".satellite";
		int sat = initializer.parseInt(hm,tmp);
		int j = whichState - 6*sat; 
		
		
		tmp = "MEAS."+measNum+".R."+j;
		double R = initializer.parseDouble(hm,tmp);
		return R;
	}
	/** Return the measurement noise value for this measurement
	 * 
	 *   
	 */
	public double R(ObservationMeasurement om)
	{
		int whichState = om.get_whichState();//EKF.stateNum;
		int measNum    = 0;//EKF.measNum;
		double R=0.0;
		
		try{
		//Ensure we are returning the correct state when there is more than
		//one satellite
		String tmp = "MEAS."+measNum+".satellite";
		int sat = initializer.parseInt(hm,tmp);
		int j = whichState - 6*sat; 
		
		
		tmp = "MEAS."+measNum+".R."+j;
		R = initializer.parseDouble(hm,tmp);
		}catch(NumberFormatException e){
			System.err.println("Warning: Number format exception - return zero noise");
			//System.exit(0);
		}catch(NullPointerException ne){
			System.err.println("Warning: Null pointer exception - return zero noise");
			//System.exit(0);
		}
		return R;
	}
	
	public VectorN H(VectorN state)
	{
		/*Determine the number of states*/
		int whichState = EKF.stateNum;
		int numStates = initializer.parseInt(hm,"FILTER.states");
		
		
		/*for a Range measurement, the current state has H = 1, all other states H = 0 */
		VectorN H = new VectorN(numStates);
		H.set(0.0);
		H.set(whichState,1.0);
		return H;
	}
	public VectorN H(ObservationMeasurement om, VectorN state)
	{
		/*Determine the number of states*/
		int whichState = om.get_whichState();//EKF.stateNum;
		int numStates = this.FILTER_states;//initializer.parseInt(hm,"FILTER.states");
		
		
		/*for a Range measurement, the current state has H = 1, all other states H = 0 */
		VectorN H = new VectorN(numStates);
		H.set(0.0);
		H.set(whichState,1.0);
		return H;
	}
	
}