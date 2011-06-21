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

package jat.forces.density.earth;
import jat.matvec.data.*;
import jat.spacecraft.Spacecraft;
import jat.spacetime.BodyRef;
import jat.spacetime.EarthTrueOfDateRef;
import jat.spacetime.ReferenceFrameTranslater;
import jat.spacetime.Time;
import jat.timeRef.*;
import jat.audio.*;
import jat.forces.*;
import jat.forces.drag.AtmosphericDrag;

/**
 * <P>
 * The CIRA_ExponentialDrag class computes the acceleration due to drag on a satellite
 * using an exponential Earth atmosphere model. The min altitude is
 * currently 200 km. To go lower, just need to add more values from 
 * the table.
 * 
 * Reference: Vallado, Table 8-4.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class CIRA_ExponentialDrag extends AtmosphericDrag {

	private static final long serialVersionUID = 1L;
	
	private int brack;
	
	private Matrix dadr;
	
	private final static double[] rho_0 = {
		2.789E-10, 7.248E-11, 2.418E-11, 9.158E-12, 3.725E-12, 1.585E-12, 
		6.967E-13, 1.454E-13, 3.614E-14, 1.170E-14, 5.245E-15, 3.019E-15};
		
	private final static double[] H = {
		37.105, 45.546, 53.628, 53.298, 58.515, 60.828, 63.822, 71.835,
		88.667, 124.64, 181.05, 268.0};
		
	private final static double[] h0 = {
		200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 600.0, 700.0, 
		800.0, 900.0, 1000.0}; 
		
	/**
	 * Constructor
	 * @param cd coefficient of drag
	 * @param area drag cross-sectional area
	 * @param mass mass
	 */	
    public CIRA_ExponentialDrag (double cd, double area, double mass){
    	super(cd, area, mass);
    }
    /**
     * Constructor 
	 * @param sc Spacecraft parameters
	 */
	public CIRA_ExponentialDrag(Spacecraft sc) {
		super(sc);
	}
    /** Compute the atmospheric density using an exponential atmosphere model.
     * @param ref EarthRef object. Not used.
     * @param r ECI position vector in meters.
     * @return Atmospheric density in kg/m^3.
     */
    public double computeDensity(EarthRef ref, VectorN r){
        r.checkVectorDimensions(3);
        
        // Get the J2000 to TOD transformation
        Matrix N = ref.TOD();
        
        // Transform r from J2000 to TOD
        VectorN r_tod = N.times(r);
        double rmag = r_tod.mag();

        // Satellite height
        Geodetic geod = new Geodetic(r_tod);
        double height = geod.getHAE()/1000.0;              //  convert to [km]
        
        // check to see if too low
        if (height < h0[0]) {
        	System.out.println("CIRA_ExponentialDrag: altitude = "+height+" too low. Min altitude = "+h0[0]);
			SoundPlayer.play("C:\\Jat\\jat\\jat\\audio\\sounds\\humanerror.wav");
			System.exit(99);

//        	return 0.0;
        }
        
        // find the right height bracket
        int n = h0.length;
        int bracket = 0;
        if (height >= h0[n-1]) {
        	bracket = n - 1;
        }
        else {
        	for (int i = 0; i < (n-1); i++) {
        		if ((height >= h0[i]) && (height < h0[i+1])){
        			bracket = i;
        		}
        	}
        }
        
        // compute the density
        this.brack = bracket;
        double rho = rho_0[bracket] * Math.exp((h0[bracket] - height)/H[bracket]);
        
//        System.out.println("ced density: "+rho);
        
        return rho;

    }
    
    /** Compute the atmospheric density using an exponential atmosphere model.
     * @param t Time reference object. Not used.
     * @param ref EarthRef object. Not used.
     * @param r ECI position vector in meters.
     * @return Atmospheric density in kg/m^3.
     */
    public double computeDensity(Time t, BodyRef ref, VectorN r) {
        r.checkVectorDimensions(3);
        
        // Translate from J2000 to TOD
        ReferenceFrameTranslater xlater =
          new ReferenceFrameTranslater(ref, new EarthTrueOfDateRef(), t);
        VectorN r_tod = xlater.translatePoint(r);
        double rmag = r_tod.mag();

        // Satellite height
        Geodetic geod = new Geodetic(r_tod);
        double height = geod.getHAE()/1000.0;              //  convert to [km]
        
        // check to see if too low
        if (height < h0[0]) {
        	System.out.println("CIRA_ExponentialDrag: altitude = "+height+" too low. Min altitude = "+h0[0]);
			SoundPlayer.play("C:\\Jat\\jat\\jat\\audio\\sounds\\humanerror.wav");
			System.exit(99);

//        	return 0.0;
        }
        
        // find the right height bracket
        int n = h0.length;
        int bracket = 0;
        if (height >= h0[n-1]) {
        	bracket = n - 1;
        }
        else {
        	for (int i = 0; i < (n-1); i++) {
        		if ((height >= h0[i]) && (height < h0[i+1])){
        			bracket = i;
        		}
        	}
        }
        
        // compute the density
        this.brack = bracket;
        double rho = rho_0[bracket] * Math.exp((h0[bracket] - height)/H[bracket]);
        
//        System.out.println("ced density: "+rho);
        
        return rho;
    }
    
    /** Computes the acceleration due to drag in m/s^2.
     * @param ref EarthRef object.
     * @param beta Satellite ballistic coefficient (Cd*A/m)
     * @param r ECI position vector in meters.
     * @param v ECI velocity vector in meters.
     * @return acceleration due to drag in m/s^2.
     */
    public void compute(EarthRef ref, VectorN r, VectorN v){

        r.checkVectorDimensions(3);
        v.checkVectorDimensions(3);
        double rmag = r.mag();
        double beta = cd * area / mass;

        // compute the atmospheric density
        double rho = computeDensity(ref, r);

        // compute the relative velocity vector and magnitude
        VectorN we = new VectorN(0.0, 0.0, omega_e);
        VectorN wxr = we.crossProduct(r);
        VectorN vr = v.minus(wxr);
        double vrmag = vr.mag();

        // form -1/2 (Cd*A/m) rho
        double coeff = -0.5 * beta * rho;
        double coeff2 = coeff * vrmag;

        // compute the acceleration in ECI frame (km/s^2)
        this.drag = vr.times(coeff2);
        
        // form partial of drag wrt v
        Matrix vrvrt = vr.outerProduct(vr);
        vrvrt = vrvrt.divide(vrmag);
        Matrix vrm = new Matrix(3);
        vrm = vrm.times(vrmag);
        this.dadv = (vrvrt.plus(vrm)).times(coeff);
        
		// form partial of drag wrt cd
		double coeff3 = coeff2 / this.cd;
		this.dadcd = vr.times(coeff3);

        // form partial of drag wrt r, see Montenbruck, p. 249
        double Hh = H[this.brack];
        double coeff4 = -1.0 / (Hh * rmag);
        VectorN drhodr = r.times(coeff4);
        Matrix part1 = vr.outerProduct(drhodr);
        part1 = part1.times(coeff2);
        
        Matrix cross = we.cross();
        Matrix part2 = this.dadv.times(cross);
        this.dadr = part1.minus(part2);        
        
    }
    
    /**
     * Return the partial derivative of acceleration wrt position
     * @return Matrix containing the partial derivative of acceleration wrt position
     */
    public Matrix partialR(){
    	return this.dadr;
    }

        
}
