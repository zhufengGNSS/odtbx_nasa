/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2006 United States Government as represented by the
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
 */
package jat.alg.estimators;

import jat.sim.*;
import jat.cm.Constants;
import java.util.HashMap;

import jat.eph.DE405;
import jat.forces.*;
import jat.forces.gravity.*;
import jat.forces.gravity.earth.*;
import jat.gps.GPS_Utils;

import jat.alg.integrators.Derivatives;
import jat.matvec.data.*;
import jat.spacetime.*;
import jat.util.FileUtil;
import java.util.Random;



/**
* This file has to be reconstructed for each state used.  It uses a 
* 4 x 4 gravity model as well as lunar and solar gravity and Solar
* radiation pressure modelds
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/ 



public class JGM4x4SRPEOM9state implements Derivatives {

	private static double re = 6378136.3; // radius of earth in meters
	private static double h_0 = 920000.0; // atmosphere model parameter
	private static double rho_0 = 4.36E-14; // atmosphere model parameter
	private static double gamma_0 = 5.381E-06; // atmosphere model parameter
	//private static double omega_e = 7.2921157746E-05; // earth rotation rate
	public static int n;
	public HashMap hm;
	public static double mass0,mass1,area0,area1;
	public static double Cr0,Cr1;
	public DE405 jpl_ephem;
	public static double mjd0;
	double j2 = Constants.j2;
	double mu = Constants.mu*1e9;
	boolean firsttime;
	UniverseModel universe;
    private Random rn;
	GravityModel earth_grav;
	
	SolarRadiationPressure srp0; 
	SolarRadiationPressure srp1 ;
	public double hc;
	
	
	public JGM4x4SRPEOM9state(HashMap hm){
		this.hm = hm;
		String fs, dir_in;
		fs = FileUtil.file_separator();
		try{
			dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"input"+fs;
		}catch(Exception e){
			dir_in = "";
		}
		//hm = initializer.parse_file(dir_in+"initialConditions.txt");
		mass0 = initializer.parseDouble(hm,"jat.0.mass");
		area0 = initializer.parseDouble(hm,"jat.0.area");
		Cr0 = initializer.parseDouble(hm,"jat.0.Cr");
		//mass1 = initializer.parseDouble(hm,"jat.1.mass");
		//area1 = initializer.parseDouble(hm,"jat.1.area");
		//Cr1 = initializer.parseDouble(hm,"jat.1.Cr");
		mjd0 = initializer.parseDouble(hm,"init.MJD0")+initializer.parseDouble(hm, "init.T0")/86400.0;

		
		hc = -2.87956633585E-10 * GPS_Utils.c;
		
		earth_grav = new GravityModel(4,4,EarthGravityType.JGM3);
		
		srp0 = new SolarRadiationPressure(mass0, area0, Cr0);
		srp1 = new SolarRadiationPressure(mass1, area1, Cr1);
		
		//Set the Gravitational parameter path
		try{
			dir_in = FileUtil.getClassFilePath("jat.eph","DE405")+"DE405data"+fs;
		}catch(Exception e){
			dir_in = "";
		}
		jpl_ephem = new DE405("dir_in");
		

		universe = new UniverseModel(mjd0);
		n = initializer.parseInt(hm,"FILTER.states");
		rn = new Random();
	}
	
	
	/**
	 * Compute the time derivatives
	 * @param t time
	 * @param y double [] containing current state.
	 * @return double[] containing the time derivatives.
	 */
	public double[] derivs(double t, double[] y) {
		
		//The out vector is sized by the number of states squared 
		//plus the number of states.  (n^2 + n)
		
		
		double out[] = new double[n*n + n];
		
		//Obtain thet the correct time
		int ctr = 0;
		Time tt = new Time(mjd0);
		tt.update(t);
		double newttt = TimeUtils.UTCtoTT(t/86400 + mjd0);
		
		
		if(firsttime == false)
		{
			earth_grav.initialize();
			firsttime = true;
			universe.set_use_iers(true);
			srp0.updateArea(area0);
			srp1.updateArea(area1);
			srp0.updateMass(mass0);
			srp1.updateMass(mass1);
			
		}
		
		// Generate some vectors for use later on
		VectorN r0 = new VectorN(y[0], y[1], y[2]);
		//VectorN v0 = new VectorN(y[3], y[4], y[5]);
		
		// store elements of incoming state in more familiar looking variables
		double xx0 = y[0];
		double yy0 = y[1];
		double zz0 = y[2];
		double vx0 = y[3];
		double vy0 = y[4];
		double vz0 = y[5];
	
		// compute derived variables
		double rmag0 = r0.mag();

		
		double rcubed0 = rmag0 * rmag0 * rmag0;
		
		double rsq0 = rmag0 * rmag0;
		
		double re_r0 = re / rmag0;
		
		re_r0 = re_r0 * re_r0;
				
		double zsq_rsq0 = (5.0 * zz0 * zz0 / rsq0) - 1.0;

		
		//Get the acceleration directly from the model
		universe.update(t);
		Matrix M = universe.earthRef.eci2ecef(universe.time.mjd_ut1(),universe.time.mjd_tt());
		VectorN acc0 = earth_grav.gravity(r0,M);
		
		//Use the acceleration model directly for integration
		//of the state.  For the generation of the state transition
		//Matirx, use the normal J2 forulation
		double ax0 = acc0.get(0);
		double ay0 = acc0.get(1);
		double az0 = acc0.get(2);
		
		
		// compute accelerations due to gravity for each satellite
		double ll0 = -1.0 * (mu * xx0 / rcubed0) * (1.0 - 1.5 * re_r0 * j2 * zsq_rsq0);
		double mm0 = -1.0 * (mu * yy0 / rcubed0) * (1.0 - 1.5 * re_r0 * j2 * zsq_rsq0);
		double nn0 = -1.0 * (mu * zz0 / rcubed0) * (1.0 - 1.5 * re_r0 * j2 * (zsq_rsq0 - 2.0));
		
		// compute accelerations due to SRP
		double AU_sqrd = Constants.AU*Constants.AU;
				
		//compute acceleration due to lunar gravity
		double ttt = TimeUtils.TTtoTDB(newttt) + 2400000.5;
        VectorN r_moon = universe.earthRef.moonVector(newttt);
        
        VectorN d0 = r0.minus(r_moon);
        
        double dmag0 = d0.mag();
        
        double dcubed0 = dmag0 * dmag0 * dmag0;
        
        VectorN temp0 = d0.divide(dcubed0);
        
        double smag = r_moon.mag();
        double scubed = smag * smag * smag;
        
        VectorN temp2 = r_moon.divide(scubed);
        VectorN sum0 = temp0.plus(temp2);
  
        VectorN lunarAcceleration0 = sum0.times(Constants.GM_Moon);
        //* TODO Watch this
        //lunarAcceleration0.set(0);
    
        //Compute the acceleration due to the solar gravity
        VectorN r_sun = universe.earthRef.sunVector(newttt);
        d0 = r0.minus(r_sun);
        
        dmag0 = d0.mag();
        
        dcubed0 = dmag0 * dmag0 *dmag0;
        
        temp0 = d0.divide(dcubed0);
        
        smag = r_sun.mag();
        double magSun3 = smag * smag * smag;
        
        temp2 = r_sun.divide(magSun3);
        
        sum0 = temp0.plus(temp2);

        VectorN solarAcceleration0 = sum0.times(Constants.GM_Sun);
        //* TODO watch
        //solarAcceleration0.set(0);
        
        //Determine if the sun is visible
        if(r0.x[0]==650869.3445073176){
        	int donothing = 0;
        }
        double visible0 = srp0.partial_illumination(r0,r_sun);
		
		double pressureConstant = 4.5344321837439e-06;
		double SRPscale0 = pressureConstant*(area0/mass0)*visible0;
        
        //Compute the relevant SRP information
		double Xsun = r_sun.get(0);
		double Ysun = r_sun.get(1);
		double Zsun = r_sun.get(2);
		//double aa = 1.0 *Cr*(Xsun/magSun3)*AU_sqrd*visible;
		//double bb = 1.0 *Cr*(Ysun/magSun3)*AU_sqrd*visible;
		//double cc = 1.0 *Cr*(Zsun/magSun3)*AU_sqrd*visible;
        
		VectorN srpacc0 = srp0.accelSRP(r0,r_sun);

		srpacc0 = srpacc0.times(visible0);
		//* TODO watch this
		//srpacc0.set(0);
		
		// compute state derivatives

        //Velocity for spacecraft 0
		out[0] = y[3];
		out[1] = y[4];
		out[2] = y[5];
		
		//acceleration for spacecraft 0
		
		out[3] = ax0   - solarAcceleration0.get(0) + srpacc0.get(0) - lunarAcceleration0.get(0);
		out[4] = ay0    - solarAcceleration0.get(1) + srpacc0.get(1) - lunarAcceleration0.get(1);
		out[5] = az0     - solarAcceleration0.get(2) + srpacc0.get(2) - lunarAcceleration0.get(2);
	
	
		
		//Clock model
		double w_f = rn.nextGaussian()*.036;
		out[6] = w_f + hc;
		double w_g = rn.nextGaussian()*7.106E-05;
		out[7] = w_g;
		
		//* TODO watch this - set clock terms to zero
		//out[6] = 0;
		//out[7] = 0;
//		if(Math.abs(w_f) > 1e3 || Math.abs(w_g) > 1e3){
//			out[6] = out[6];
//		}
	
		//Solar radiation Pressure states
		out[8] = 0;
		
		
		// compute A matrix
		
		Matrix a = new Matrix(n,n); // creates a (15x15 matrix with all zeros

	
		//Generate some utility variables.  It is a bit messy as we will need
		//two now.
		
		double r50 = rsq0 * rcubed0;

		
		double mur50 = mu / r50;

		
		double mur30 = mu / rcubed0;
		
		
		double sz2r20 = 7.0 * zz0 * zz0 / rsq0;

		double muxyr50 = mu * xx0 * yy0 / r50;
		
		double muxzr50 = mu * xx0 * zz0 / r50;
		
		double muyzr50 = mu * yy0 * zz0 / r50;

		double bracket10 = 3.0 - 7.5 * re_r0 * j2 * (sz2r20 - 1.0);
		
		double bracket30 = 3.0 - 7.5 * re_r0 * j2 * (sz2r20 - 3.0);
		
		double bracket50 = 3.0 - 7.5 * re_r0 * j2 * (sz2r20 - 5.0);
		
		double bracket20 = 1.5 * re_r0 * (5.0 * zz0 * zz0 / rsq0 - 1.0);
		
		//Note:  use this formulation for ll to avoid a singularity
		ll0 = -1.0 * (mu  / rcubed0) * (1.0 - 1.5 * re_r0 * j2 * zsq_rsq0);
		
		double dldx0 = ll0 + mur50 * xx0 * xx0 * bracket10;
		
		double dldy0 = muxyr50 * bracket10;
		
		double dldz0 = muxzr50 * bracket30;
		
		double dldj20 = mur30 * xx0 * bracket20;
		
		double dmdx0 = dldy0;		
		
		//Note:  use a different definition of mm to avoid singularity issues
		
		mm0 = -1.0 * (mu / rcubed0) * (1.0 - 1.5 * re_r0 * j2 * zsq_rsq0);
		
		double dmdy0 = mm0  + mur50 * yy0 * yy0 * bracket10;
		
		double dmdz0 = muyzr50 * bracket30;
		
		double dmdj20 = mur30 * yy0 * bracket20;
		
		double dndx0 = muxzr50 * bracket30;
		
		double dndy0 = dmdz0;
		
		//Note:  Use this definition of nn0.  It is uglier, but it removes a potential singularity
		nn0 = -1.0 * (mu / rcubed0) * (1.0 - 1.5 * re_r0 * j2 * (zsq_rsq0 - 2.0));
		
		//double dndz0 = nn0 / zz0 + mur50 * zz0 * zz0 * bracket50;
		double dndz0 = nn0 + mur50 * zz0 * zz0 * bracket50;
		
		
		double dndj20 = mur30 * zz0 * (1.5 * re_r0 * (5.0 * zz0 * zz0 / rsq0 - 3.0));

		SRPscale0 = 0;
		
		a.A[0][3] = 1.0;
		a.A[1][4] = 1.0;
		a.A[2][5] = 1.0;
		
		a.A[3][0] = dldx0 ;
		a.A[3][1] = dldy0 ;
		a.A[3][2] = dldz0 ;
		a.A[3][3] = 0;
		a.A[3][4] = 0;
		a.A[3][5] = 0;
		//a.A[3][8] = 1.0 *(Xsun/magSun3)*AU_sqrd*SRPscale0;

		
		a.A[4][0] = dmdx0;
		a.A[4][1] = dmdy0;
		a.A[4][2] = dmdz0;
		a.A[4][3] = 0;
		a.A[4][4] = 0;
		a.A[4][5] = 0;
		//a.A[4][8] = 1.0 *(Ysun/magSun3)*AU_sqrd*SRPscale0;

		a.A[5][0] = dndx0;
		a.A[5][1] = dndy0;
		a.A[5][2] = dndz0;
		a.A[5][3] = 0;
		a.A[5][4] = 0;
		a.A[5][5] = 0;
		//a.A[5][8] = 1.0 *(Zsun/magSun3)*AU_sqrd*SRPscale0;
		

		
		
		// compute phi derivatives

		Matrix phi = this.phi(y);

		Matrix dphi = a.times(phi);

		//	   dphi.print("dphi");

		// put phi derivatives into output array
		int k = n;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				out[k] = dphi.A[i][j];
				k = k + 1;
			}
		}

		return out;
	}

	private Matrix phi(double[] in) {
		Matrix out = new Matrix(n, n);
		int k = n;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				out.A[i][j] = in[k];
				k = k + 1;
			}
		}
		return out;
	}

}
