package jat.gps;

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
 * File Created on Jul 13, 2003
 */
 
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.StringTokenizer;
import jat.gps.*;
import jat.alg.integrators.LinePrinter;
import jat.cm.FiniteBurn;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;
import jat.math.*;
import jat.timeRef.RSW_Frame;
import jat.util.FileUtil;
//import jat.math.*;

/**
 * <P>
 * The Geo_Blockage_Model Class provides a model of GPS signal blockage due to 
 * a spherical earth as well as the ability to include signal strength
 * effects for sidelobe signals
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class GEO_Blockage_Models implements ExpandedVisible {
	
	private double elevationMask;
	private static final double earthRadius = 6478140.0;
	private static final double earthRadiusWithIono = 7478140.0;
	private double L1freq = 1575.42e6;
	private double C      = 2.99792458e8;
	
	double Ae = 0;  //Attenuation due to the atmosthpere
	double As = 0 ;  //To account for signal loss in front of LNA
	double PowerSV = 14.9;  //Average transmit power of block IIA sats
	double L  = -1.1;  //Receiver implimentation loss
	double Nf = -3.4;  //Noise figure for Receiver/LNA
	double Ts =  290.0;  //System Noise Temperature (K) earth pointing antenna
	
	
	double pointingBias = 0;
	//The receiver performance figures
	double AcquisitionThreshold = 25;// (dB-Hz)
	double TrackingThreshold    = 22;// (dB-Hz)
	static boolean [] AcquisitionFlag = new boolean[33];
	
	static double [] receiverAntennaDeg = new double[15];
	static double [] receiverAntennaGain = new double[15];
	static double [] GPSAntennaDeg = new double[66];
	static double [] GPSAntennaGain = new double[66];
	public static double Cn0;
	public static double elevation;
	public static double trackedCn0;
	boolean firsttime = true;

	
    /**
     * Determine if the GPS satellite is visible, including earth blockage
     * Used by GPS Measurement Generator.
     * @param losu Line of sight unit vector
     * @param r    current position vector
     * @param rGPS current position vector of GPS Satellite
     * @return boolean true if the GPS satellite is visible
     */
	public boolean visible(VectorN losu, VectorN r, VectorN rGPS) {

		//Load the antenna gain maps.  Note:  0 = receiver, 1 = GPS Sv
		//NOTE: I believe that both GPS satellite and the receiver antenna
		//are currently pointed at the center of the Earth.   (ie, zero degree
		//is pointed at the center of the Earth)
		if(firsttime == true)
		{
            String fs, dir_in;
            fs = FileUtil.file_separator();
            try{
                dir_in = FileUtil.getClassFilePath("jat.gps","GEO_Blockage_Models")+"antennas"+fs;
            }catch(Exception e){
                dir_in = "C:/Code/Jat/jat/eph/DE405data/";
            }	
//          
			//readFromFile("C:\\GOESR\\omni_antenna.txt",receiverAntennaDeg,receiverAntennaGain);
			readFromFile(dir_in+"ballhybrid_10db_60deg.txt",receiverAntennaDeg,receiverAntennaGain);
			//readFromFile(dir_in+"patch.txt",receiverAntennaDeg,receiverAntennaGain);
			readFromFile(dir_in+"GPSIIA_L1MEAN.txt",GPSAntennaDeg,GPSAntennaGain);
			//readFromFile("C:\\GOESR\\LMantenna.txt",receiverAntennaDeg,receiverAntennaGain);
			firsttime = false;
		}
		Interpolator interpR = new Interpolator(receiverAntennaDeg,receiverAntennaGain);
		Interpolator interpG = new Interpolator(GPSAntennaDeg,GPSAntennaGain);
		// check elevation mask
		boolean visible = true;
		
		// First check for Earth Obscuration
		double dist = r.mag();
		
		//  Cut off the ionosphere by adding 1000 Km to the 
		//  Mask.  This negates having to model/correct for
		//  the ionosphere
		
		
		//double cone = Math.atan2(earthRadius, dist);
		double cone = Math.atan2(earthRadiusWithIono, dist);
		
		VectorN r_unit = r.unitVector().times(-1.0);
		double cos_delta = r_unit.dotProduct(losu);
		double delta = Math.acos(cos_delta);
		if (delta < cone) {
			visible = false;
		}
		
		//Only continue if the satellite isn't behind the Earth
		if(visible){
	
			//Determine the Angle between the GPS SV boresight 
			//and the SV receiver
			VectorN unitrGPS = rGPS.unitVector();
			VectorN minusUnitRGPS = unitrGPS.times(-1.0);
			VectorN minusLosu = losu.times(-1.0);
			double thetaGPS = (180/(Math.PI))*Math.acos(minusUnitRGPS.dotProduct(minusLosu));
			if(thetaGPS > 40 || thetaGPS < -40)
			{
				visible = false;
			}
			
			//Lookup the appropriate gain from the transmitting antenna
			double Gt = interpG.get_value(thetaGPS);
			
			//Determine the angle between the receiver boresight 
			//and the GPS SV
//			Generate a RTN rotation matrix and rotate the LOS vector into it
			//Matrix T = RSW_Frame.ECI2RIC(r,v);
			//VectorN losRTN = T.times(losu);
			
			//azimuth[prn] = ((180)/Math.PI)*Math.acos(dots/mags);
			//azimuth[prn] = ((180)/Math.PI)*Math.atan2(losRTN.get(0),losRTN.get(1)*-1);

//			Compute the Elevaton to the GPS Satellite
			//elevation[prn]=(180/Math.PI)*Math.asin(losRTN.get(2)/losRTN.mag());
			
			
			
			
			
			VectorN unitR = r.unitVector();
			VectorN minusUnitR = unitR.times(-1.0);
			double thetaR = (180/(Math.PI))*Math.acos(minusUnitR.dotProduct(losu)) + pointingBias;
			elevation = thetaR;
			if(thetaR > 40 || thetaR < -40)
			{
				visible = false;
			}
			
			
			//Lookup the approprate gain from the receiving antenna			
			double Gr = interpR.get_value(thetaR);
			
			//Determine the free space propagation loss
			VectorN los = GPS_Utils.lineOfSight(r, rGPS);
			double Ad = 20*Math.log10((C/L1freq) / (4*Math.PI*los.mag()));
			
			//Determine the attenuated power
			// Attenuated power = transmitted power + transmitted gain 
			//                  + free space loss + atmospheric disturbances
			double Ap = PowerSV + Gt + Ad + Ae;
			
			//Determine the received power
			// Received Power = Attenuated Power + Receiver Gain + 
			//  				Losses prior to LNA
			double Rp = Ap + Gr + As;
			
			
			//Determine the Carrier to Noise Ratio
			//Carrier to Noise ratio = Received power - scaling of system noise
			//						 +  288.6 (not sure what that is) + receiver noise figure
			//						 +  Losses in receiver/LNA (front end?)
			Cn0 = Rp - 10*Math.log10(Ts) + 228.6 + Nf + L;
			
			
			
			//If the satellite isn't visible, reset the AcquisitionFlag
			//to force it to reacquire
			//if(visible == false)
			//	AcquisitionFlag[prn] = false;
			
			
			//System.out.println("Cn0 for PRN : " + prn + " is: "+Cn0);
			//If it is visible and strong enough, "Acquire" satellite
			//if(visible && AcquisitionFlag[prn] == false && Cn0 > AcquisitionThreshold)
			//	AcquisitionFlag[prn] = true;

			//If we have "acquired" a satellite and if the Cn0 is greater than
			//the tracking threshold, then it is visible
			//if(AcquisitionFlag[prn])
		//	{
			//	if(Cn0 > TrackingThreshold)
			//	{
			//		visible = true;
			//		trackedCn0 = Cn0;
			//	}
			//}
			//else
			//{
			//	visible = false;
			//	trackedCn0 = 0;
			//}
		
		
		}
		if(visible != true)
			elevation = 0.0;
//		if(visible == true)
//		{
//			FileWriter out;
//			try {
//				out = new FileWriter("C:\\GOESR\\tmp.txt",true);
//				BufferedWriter writer = new BufferedWriter(out);
//				String pp = " " + Cn0 + "\n";
//				writer.write(pp);
//				writer.close();
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
		return visible;
		
	}
	public void readFromFile(String file, double [] antennaDeg, double [] antennaGain) {
		try {
			FileReader fr = new FileReader(file);
			BufferedReader in = new BufferedReader(fr);
			String line;
			int j = 0;
			// loop through the file, one line at a time
			while ((line = in.readLine()) != null) {
				StringTokenizer tok = new StringTokenizer(line, "\t");
				int total = tok.countTokens();
				
				// check for consistent number of columns
				if (total != 2) {
					System.out.println(
					"AntennaGainPattern Error: Number of columns do not match");
					System.exit(-99);
				}
				//antennaDeg[0] = 1.0;
				//Parse out the valuse
				String token = tok.nextToken();
				antennaDeg[j] = Double.parseDouble(token);
				token = tok.nextToken();
				antennaGain[j] = Double.parseDouble(token);
				
				//System.out.println("Degree: "+ antennaDeg[j] + "Gain: " + antennaGain[j]);
				j++;
			}	
			in.close();
			fr.close();
		} catch (IOException e) {
			System.err.println("Error opening:" + file);
			return;
		}
	}


	public boolean visible(VectorN losu, VectorN r, VectorN rGPS, VectorN v, int prn) {
		// TODO Auto-generated method stub
		return false;
	}
	public boolean visible(VectorN losu, VectorN r, VectorN v, VectorN rGPS, VectorN vGPS, int prn, double mjd) {
		// TODO Auto-generated method stub
		return false;
	}

}

