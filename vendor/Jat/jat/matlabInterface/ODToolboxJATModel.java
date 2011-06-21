/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2005 United States Government as represented by the
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
package jat.matlabInterface;

import java.io.Serializable;
import jat.spacecraft.*;
import jat.spacetime.UniverseModel;
import jat.matvec.data.VectorN;
import jat.matvec.data.Matrix;
import jat.alg.integrators.*;
import jat.cm.Constants;
import jat.forces.GravitationalBody;
import jat.forces.gravity.*;
import jat.forces.gravity.earth.*;
import jat.forces.density.earth.*;
import jat.forces.drag.AtmosphericDrag;
import jat.forces.Moon;
import jat.forces.SolarRadiationPressure;
import jat.forces.Sun;

/**
 * This is the primary class for the OD Toolbox interface. Each instantiation
 * encapsulates information about a spacecraft and the universe force models in
 * use.
 * 
 * @author Derek M. Surka
 * 
 */
public class ODToolboxJATModel implements Derivatives, Serializable {

	private static final long serialVersionUID = 1L;

    public UniverseModel spacetime;
    /**
     * Spacecraft model used when propagating a single spacecraft
     */
    public Spacecraft sc;
    /**
     * The starting epoch in Modified Julian Date Universal Coordinated Time (UTC)
     */
    public double mjd_utc_start;
    
    /**
     * Flag for computing 2-body gravity partials
     */
    private boolean use_j2_gravity 		= false;

    /**
     * Flag for computing J2 gravity partials
     */
    private boolean use_2body_gravity 	= false;

    /**
     * Flag for computing drag partials
     */
    private boolean use_drag 			= false;
       
    /**
     * Index of drag atmosphere model for computing drag partials
     */
    private int drag_index;
       
    /**
     * Constructor initializes a single spacecraft simulation given relevant parameters
     * @param cr Coefficient of Reflectivity used for Solar Radiation Pressure
     * @param cd Coefficient of Drag used for Atmospheric Drag calculations
     * @param dragArea Cross sectional area used for drag
     * @param srpArea  Cross sectional area used for Solar Radiation Pressure
     * @param mass Mass of the spacecraft
     * @param utc Starting epoch, Modified Julian Date in Universal Coordinated Time
     */
    public ODToolboxJATModel(double cr, double cd, double dragArea, double srpArea, double mass, double utc){
        sc = new Spacecraft();
        sc.set_dragArea(dragArea);
        sc.set_srpArea(srpArea);
        sc.set_cd(cd);
        sc.set_cr(cr);
        sc.set_mass(mass);
        sc.set_use_params_in_state(false);
        spacetime = new UniverseModel(utc);
    }
     
    /**
     * Initialize the forces present in the universe model during the Simulation.
     * @param force_flag An array of boolean values indicating the forces in order:
     *  [0] true: Solar gravity 
     *  [1] true: Lunar gravity 
     *  [2] true: Atmospheric drag
     *  [3] true: Solar Radiation Pressure
     * @param grav_model String specifiying the gravity model to use
     * @param degree of the gravity model
     * @param order of the gravity model
     * @param atm_model "NRL" for NRLMSISE2000 or "HP" for Harris Priester
     * @param atm_params Array of doubles containing parameters for atmosphere models
     *  [0] n:         HP inclincation parameter ranging from 2(low inclination) to 6(high)
     *  [1] f107Daily: Daily value of 10.7cm solar flux
     *  [2] f107Avg:   81-day aveage of 10.7cm solar flux
     *  [3] ap:		   Geomagnetic flux index
     */
    public void initializeForces(boolean[] force_flag, String grav_model, int degree, int order,
			String atm_model, double[] atm_params) {

		VectorN zero = new VectorN(0, 0, 0);

		if (grav_model.equalsIgnoreCase("2Body")) {
			GravitationalBody earth = new GravitationalBody(Constants.GM_Earth);
			spacetime.addForce(earth);
			use_2body_gravity = true;
			
		} else { // DMS Is there a better way to find the gravity type?
			Boolean foundModel = false;
			EarthGravityType[] index = EarthGravityType.index;
			for (int i = 0; i < EarthGravityType.num_models; i++) {
				if (grav_model.equalsIgnoreCase(index[i].toString())) {
					GravityModel earth_grav = new GravityModel(degree, order,
							index[i]);
					spacetime.addForce(earth_grav);
					spacetime.set_use_iers(true);
					foundModel = true;
					use_j2_gravity = true;
					break;
				}
			}
			if (!foundModel) {
				System.out
						.println("Warning: GravityModel "
								+ grav_model
								+ " type not found. No Earth gravity forces will be used.");

			}
		}
        
		if (force_flag[0]) {
	        spacetime.set_compute_sun(true);
	        Sun sun =
	            new Sun(Constants.GM_Sun,zero,zero);
	        spacetime.addForce(sun);
	    }
	    if(force_flag[1]){
	        spacetime.set_compute_moon(true);
	        Moon moon =
	            new Moon(Constants.GM_Moon,zero,zero);
	        spacetime.addForce(moon);
	    }
	    if (force_flag[2]) {
			double n_param 		= atm_params[0];  
			double f107Daily 	= atm_params[1];
			double f107Average 	= atm_params[2];
			double ap 			= atm_params[3];
			
			spacetime.set_compute_sun(true);
			
			if (atm_model.equalsIgnoreCase("NRL")) {
				NRLMSISE_Drag drag = new NRLMSISE_Drag(sc);
				
				drag.setF107Daily(f107Daily);
				drag.setF107Average(f107Average);
				drag.setAP(ap);
				
				spacetime.addForce(drag);
				spacetime.set_use_iers(true);

				use_drag 	= true;
				drag_index 	= spacetime.getForceSize() - 1;
				
			} else if (atm_model.equalsIgnoreCase("HP")) {
				HarrisPriester atmos = new HarrisPriester(sc, 150);
				atmos.setParameter(n_param);
				atmos.setF107(f107Daily);

				spacetime.addForce(atmos);
				spacetime.set_use_iers(true);

				use_drag 	= true;
				drag_index 	= spacetime.getForceSize() - 1;
				
			} else {
				System.out
						.println("Warning: AtmosphereModel "
								+ atm_model
								+ " type not found. No atmospheric drag forces will be used.");
	        	
	        }
	        	
	    }
	    if(force_flag[3]){
	        spacetime.set_compute_sun(true);
	        SolarRadiationPressure srp = new SolarRadiationPressure(sc);
	        spacetime.addForce(srp);
	    }
    }

    /**
     * Implements the Derivatives interface for use with the RungeKutta integrator.
     * Given a time in seconds since epoch and a state vector, return the derivatives
     * of the state vector.
     *
     * @param t Time since epoch [sec]
     * @param X State vector. [x y z xdot ydot zdot] in meters and m/s
     */
    public double[] derivs(double t, double[] x) {
        //* Update spacecraft
    	sc.updateState(x);    
        spacetime.update(t);  //DMS this is necessary because updates time, earthRef, and IERS
        //* Get non-control derivatives
        double[] xdot = spacetime.derivs(t,sc); //DMS this only updates time (duplicate)
        return xdot;
    }

    /**
     * Get the gravity gradient partials that are computed here for convenience.
     * 2-body gradient is from r2bp.m
     * J2 gradient is from jat.forces.J2Gravity
     * These should be computed in gravity files.
     * @param r Spacecraft position vector, meters
     * @return 3x3 array of gravity partials
     */
    public double[][] gravityPartials(double[] r){
    	double j2 = 0.00108263;
        double mu = Constants.GM_Earth;  // [m^3/s^2]; JGM3
        double re = Constants.rEarth_STKJGM3;  // Radius Earth [m]; JGM3;

		// strip out components of r
		double xx = r[0];
		double yy = r[1];
		double zz = r[2];
		
		// compute some intermediate values
		double rmag 	= Math.sqrt(xx*xx + yy*yy + zz*zz);
		double rsq 		= rmag * rmag;
		double rcubed 	= rsq * rmag;
		double r5 		= rsq * rcubed;
		
		double mur5 = mu / r5;
		double mur3 = mu / rcubed;

		double re_r = re / rmag;
		re_r = re_r * re_r;

		// Initialize partials to 0
		double dldx = 0;
		double dldy = 0;
		double dldz = 0;
		double dmdx = 0;
		double dmdy = 0;
		double dmdz = 0;
		double dndx = 0;
		double dndy = 0;
		double dndz = 0;
        
        if( use_2body_gravity ){
    		dldx = (3.0 * xx * xx / rsq - 1) * mur3;
    		dldy = 3.0 * xx * yy * mur5;
    		dldz = 3.0 * xx * zz * mur5;
    		
    		dmdx = dldy;
    		dmdy = (3.0 * yy * yy / rsq - 1) * mur3;
    		dmdz = 3.0 * yy * zz * mur5;
    		
    		dndx = dldz;
    		dndy = dmdz;
    		dndz = (3.0 * zz * zz / rsq - 1) * mur3;    		
    	}
    	else {     		
    		// compute some intermediate values
       		double zsq_rsq = (5.0 * zz * zz / rsq) - 1.0;		
       		double sz2r2 = 7.0 * zz * zz / rsq;
    		double muxyr5 = mu * xx * yy / r5;
    		double muxzr5 = mu * xx * zz / r5;
    		double muyzr5 = mu * yy * zz / r5;
    		double bracket1 = 3.0 - 7.5 * re_r * j2 * (sz2r2 - 1.0);
    		double bracket3 = 3.0 - 7.5 * re_r * j2 * (sz2r2 - 3.0);
    		double bracket5 = 3.0 - 7.5 * re_r * j2 * (sz2r2 - 5.0);

    		// compute "accelerations" due to gravity (these would be multiplied by r[i]
    		double ll = -1.0 * mur3 * (1.0 - 1.5 * re_r * j2 * zsq_rsq);
    		double mm = -1.0 * mur3 * (1.0 - 1.5 * re_r * j2 * zsq_rsq);
    		double nn = -1.0 * mur3 * (1.0 - 1.5 * re_r * j2 * (zsq_rsq - 2.0));
    		
    		// partials of ll wrt r
    		dldx = ll + mur5 * xx * xx * bracket1;
    		dldy = muxyr5 * bracket1;
    		dldz = muxzr5 * bracket3;
    		
    		// partials of mm wrt r
    		dmdx = dldy;
    		dmdy = mm + mur5 * yy * yy * bracket1;
    		dmdz = muyzr5 * bracket3;

    		// partials of nn wrt r
    		dndx = dldz;
    		dndy = dmdz;
    		dndz = nn + mur5 * zz * zz * bracket5;   		
    	}
        
		// set the gravity gradient matrix
    	double[][] partials = new double[3][3];
     	partials[0][0] = dldx;
     	partials[0][1] = dldy;
     	partials[0][2] = dldz;
     	partials[1][0] = dmdx;
     	partials[1][1] = dmdy;
     	partials[1][2] = dmdz;
     	partials[2][0] = dndx;
     	partials[2][1] = dndy;
     	partials[2][2] = dndz;
     	
    	return partials;
    }
    
    /**
     * Get the atmospheric drag partials wrt v
     * @return 3x3 array of drag partials
     */
    public double[][] dragPartials(){
    	AtmosphericDrag atmModel = (AtmosphericDrag)spacetime.getForce(drag_index); 
    	Matrix dadv = atmModel.partialV();
    	
    	double[][] partials = new double[3][3];
     	partials[0][0] = dadv.get(0,0);
     	partials[0][1] = dadv.get(0,1);
     	partials[0][2] = dadv.get(0,2);
     	partials[1][0] = dadv.get(1,0);
     	partials[1][1] = dadv.get(1,1);
     	partials[1][2] = dadv.get(1,2);
     	partials[2][0] = dadv.get(2,0);
     	partials[2][1] = dadv.get(2,1);
     	partials[2][2] = dadv.get(2,2);
     	
    	return partials;
    }
    
    /**
     * Set the starting epoch.
     * @param mjd_utc UTC Epoch in MJD
     */
    public void set_epoch(double mjd_utc){
        spacetime.set_time(mjd_utc);
    }

    /**
     * Get the 2-body partials flag
     * @return 2-body partials flag
     */
    public boolean get_2body_gravity_flag(){
    	return use_2body_gravity;
    }

    /**
     * Get the J2 partials flag
     * @return J2 partials flag
     */
    public boolean get_j2_gravity_flag(){
    	return use_j2_gravity;
    }

    /**
     * Get the drag partials flag
     * @return drag partials flag
     */
    public boolean get_drag_flag(){
    	return use_drag;
    }

}
