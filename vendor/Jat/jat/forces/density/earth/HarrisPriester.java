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
import jat.timeRef.*;
import jat.timeRef.EarthRef;
import jat.forces.drag.AtmosphericDrag;
import jat.forces.harrispriester.*;
import jat.spacetime.*;
//import jat.eph.*;

/**
 * <P>
 * The HarrisPriester class computes the Harris-Priester atmosphere model.
 * This code is from Montenbruck. Good for 100 - 1000 km altitude only.
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @author Richard C. Page III
 * @version 1.0
 */

public class HarrisPriester extends AtmosphericDrag{

	private static final long serialVersionUID = 1L;

    private static final double upper_limit = 2000.0;           // Upper height limit [km]
    private static final double lower_limit =    100.0;           // Lower height limit [km]
    private static final double ra_lag      = 0.523599;           // Right ascension lag [rad]

    protected double[][] data;
    protected int data_n=0;
    protected int h_index = 0;
    
    private double n_prm = 2;           // Harris-Priester parameter
    // 2(6) low(high) inclination
    /**
     * F10.7 cm Mean Solar Flux
     */
    protected double f107=157;
    
	/**
	 * Constructor
	 * @param cd coefficient of drag
	 * @param area drag cross-sectional area
	 * @param mass mass
	 */	
    public HarrisPriester (double cd, double area, double mass){
    	super(cd, area, mass);
    	initialize();
    }

    /**
     * Constructor 
	 * @param sc Spacecraft parameters
	 */
	public HarrisPriester(Spacecraft sc) {
		super(sc);
		VectorN h = sc.r().crossProduct(sc.v());
		double i = Math.acos(h.x[2] / h.mag()); // inclination
		//if(i>Math.PI/2) i = Math.PI-i;
		this.n_prm = 2.0+i*8.0/Math.PI;
		initialize();
	}
	
	/**
	 * Constructor
	 * @param sc Spacecraft parameters
	 * @param f Value of F10.7 for the density model
	 */
	public HarrisPriester(Spacecraft sc, double f){
	    super(sc);
	    f107 = f;
	    VectorN h = sc.r().crossProduct(sc.v());
		double i = Math.acos(h.x[2] / h.mag()); // inclination
		this.n_prm = 2.0+i*8.0/Math.PI;
		initialize();
	}
	
	/**
	 * Initialize the data file and store the data
	 *
	 */
	private void initialize(){
	    HarrisPriesterData hpdata = new HarrisPriesterData();
	    hpdata.process(f107);
	    data = hpdata.get_data();
	    data_n = hpdata.get_data_n();
	}
	
	/**
	 * Set the Harris Priester parameter used during density calculation.
	 * Use 2 for low inclination orbits and 6 for polar orbits.
	 * The default (without calling setParameter) is 2. 
	 * @param n
	 */
	public void setParameter(double n){
	    n_prm = n;
	}
	
	/**
	 * Sets the value of F10.7 which determines the density model
	 * @param f
	 */
	public void setF107(double f){
	    f107 = f;
	    initialize();
	}
	/**
	 * Get the minimum and maximum densities at the given altitude
	 * @param alt True Altitude [km]
	 * @return [min_density max_density]
	 */
    public double[] get_min_max(double alt){
        h_index = alt_index(alt);
        double[] out = new double[4];
        if(alt != data[h_index][0]){
            out[0] = data[h_index][1];//+(data[h_index+1][1]-data[h_index][1])*(alt-data[h_index][0])/
            						//(data[h_index+1][0]-data[h_index][0]);
            out[1] = data[h_index][2];//+(data[h_index+1][2]-data[h_index][2])*(alt-data[h_index][0])/
									//(data[h_index+1][0]-data[h_index][0]);
            out[2] = data[h_index+1][1];
            out[3] = data[h_index+1][2];
        } else {
            out[0] = data[h_index][1];
            out[1] = data[h_index][2];
        }
        return out;
    }


	
	/**
	 * Search the data for the particular altitude
	 * @param a The altitude to find in the table.
	 * @return The index of the altitude which is the
	 * last altitude which does not excede the argument.
	 */
	private int alt_index(double a){
	    int index = 0;
	    for(int i=0; i<data_n-1; i++){
	        if(a >= data[i][0] && a < data[i+1][0]){
	            index = i;
	            break;
	        }
	    }
	    return index;
	}
	
//	private double computeMin(double f, double alt){
//	    int i = alt_index(alt);
//	    if(i==hf.length) return 0;
//	    double pmin1=0, pmin2=0, pmin=0;
//	    double k=0;
//        pmin1 = Amin[i]*Math.pow((f-65),Alphamin[i])+
//       		Bmin[i]*(2-Math.exp(-Betamin[i]*(f-65)));
//        pmin2 = Amin[i+1]*Math.pow((f-65),Alphamin[i+1])+
//   			Bmin[i+1]*(2-Math.exp(-Betamin[i+1]*(f-65)));
//        k = (alt-hf[i])/(hf[i+1]-hf[i]);
//        pmin = pmin1*Math.pow((pmin2/pmin1),k);
//	    return pmin;
//	}
//
//	private double computeMax(double f, double alt){
//	    int i = alt_index(alt);
//	    if(i==hf.length) return 0;
//	    double pmax1=0, pmax2=0, pmax=0;;
//	    double k=0;
//        pmax1 = Amax[i]*Math.pow((f-65),Alphamax[i])+
//       		Bmax[i]*(2-Math.exp(-Betamax[i]*(f-65)));
//        pmax2 = Amax[i+1]*Math.pow((f-65),Alphamax[i+1])+
//   			Bmax[i+1]*(2-Math.exp(-Betamax[i+1]*(f-65)));
//        k = (alt-hf[i])/(hf[i+1]-hf[i]);
//        pmax = pmax1*Math.pow((pmax2/pmax1),k);
//	    return pmax;
//	}
	

	

	/**
	 * 
	 * Compute the atmospheric density using the Harris-Priester atmosphere model.
     * @param ref Earth reference object.
     * @param r ECI position vector in meters.
     * @return Atmospheric density in kg/m^3.
     */
	public double computeDensity(Time t, BodyRef ref, VectorN r){
	    double density = 0;
	    r.checkVectorDimensions(3);
        // Translate from J2000 to TOD
        ReferenceFrameTranslater xlater =
          new ReferenceFrameTranslater(ref, new EarthTrueOfDateRef(), t);
        VectorN r_tod = xlater.translatePoint(r);
        double rmag = r_tod.mag();
        //* Variables
        int    i, ih;                              // Height section variables
        double dec_Sun, ra_Sun, c_dec;             // Sun declination, right asc.
        //* Satellite true altitude
      //Matrix eci2ecef = ref.eci2ecef();     //*debug
      //VectorN r_ecef = eci2ecef.times(r);   //*debug
      //Geodetic geod = new Geodetic(r_ecef); //*debug  - no change
        Geodetic geod = new Geodetic(r_tod);
        double alt = geod.getHAE()/1000.0;		// [km]
        
        //alt = 360.013234313;
        
	    //*	Exit with zero density outside height model limits
        if ( alt >= upper_limit || alt <= lower_limit )
           return 0.0;
        //* Sun right ascension, declination
        VectorN r_Sun = ref.get_JPL_Sun_Vector(t);
        ra_Sun  = Math.atan2( r_Sun.x[1], r_Sun.x[0] );
        dec_Sun = Math.atan2( r_Sun.x[2], Math.sqrt( Math.pow(r_Sun.x[0],2)+Math.pow(r_Sun.x[1],2) ) );
        //* Unit vector u towards the apex of the diurnal bulge
        //* in inertial geocentric coordinates
        c_dec = Math.cos(dec_Sun);
        double [] u_ = new double[3];  // Apex of diurnal bulge
        u_[0] = c_dec * Math.cos(ra_Sun + ra_lag);
        u_[1] = c_dec * Math.sin(ra_Sun + ra_lag);
        u_[2] = Math.sin(dec_Sun);
        VectorN u = new VectorN(u_);
        //* Cosine of half angle between satellite position vector and
        //* apex of diurnal bulge
        double c_psi2 = (0.5 + 0.5 * r.dotProduct(u)/r.mag());
        //* Density interpolation
        //double pmin = computeMin(f107,alt);	//[kg/km^3]
        //double pmax = computeMax(f107,alt); //[kg/km^3]
        double[] minmax = this.get_min_max(alt);
        double cmin0 = minmax[0];
        double cmax0 = minmax[1];
        double cmin1 = minmax[2];
        double cmax1 = minmax[3];
        double h0 = data[h_index][0];
        double h1 = data[h_index+1][0];
        //* Density computation
//        double pmin = minmax[0];
//        double pmax = minmax[1];
//        density = pmin+(pmax-pmin)*Math.pow(c_psi2,n_prm/2);
//        System.out.println("hp density: "+density);
        
        
        // Height index search and exponential density interpolation

        double h_min = ( h0 - h1 )/Math.log( cmin1/cmin0 );
        double h_max = ( h0 - h1 )/Math.log( cmax1/cmax0 );

        double d_min = cmin0 * Math.exp( (h0-alt)/h_min );
        double d_max = cmax0 * Math.exp( (h0-alt)/h_max );

        // Density computation

        density = d_min + (d_max-d_min)*Math.pow(c_psi2,n_prm/2);
        return density * 1.0e-9;

        //      * Convert units and return density
//        double rho1 = 0; //-0.0045; //-0.052;
//        return density * 1.0e-9 * (1+rho1);       // [kg/m^3]	    
	}
	

	/**
	 * 
	 * Compute the atmospheric density using the Harris-Priester atmosphere model.
	 * @deprecated
     * @param ref Earth reference object.
     * @param r ECI position vector in meters.
     * @return Atmospheric density in kg/m^3.
     */
	public double computeDensity(EarthRef ref, VectorN r){
	    double density = 0;
	    r.checkVectorDimensions(3);
        //* Get the J2000 to TOD transformation
        Matrix N = ref.trueOfDate();
        //* Transform r from J2000 to TOD
        VectorN r_tod = N.times(r);
        double rmag = r_tod.mag();
        //* Variables
        int    i, ih;                              // Height section variables
        double dec_Sun, ra_Sun, c_dec;             // Sun declination, right asc.
        //* Satellite true altitude
      //Matrix eci2ecef = ref.eci2ecef();     //*debug
      //VectorN r_ecef = eci2ecef.times(r);   //*debug
      //Geodetic geod = new Geodetic(r_ecef); //*debug  - no change
        Geodetic geod = new Geodetic(r_tod);
        double alt = geod.getHAE()/1000.0;		// [km]
        
        //alt = 360.013234313;
        
	    //*	Exit with zero density outside height model limits
        if ( alt >= upper_limit || alt <= lower_limit )
           return 0.0;
        //* Sun right ascension, declination
        VectorN r_Sun = ref.get_JPL_Sun_Vector(new Time(ref.mjd_utc()));
        ra_Sun  = Math.atan2( r_Sun.x[1], r_Sun.x[0] );
        dec_Sun = Math.atan2( r_Sun.x[2], Math.sqrt( Math.pow(r_Sun.x[0],2)+Math.pow(r_Sun.x[1],2) ) );
        //* Unit vector u towards the apex of the diurnal bulge
        //* in inertial geocentric coordinates
        c_dec = Math.cos(dec_Sun);
        double [] u_ = new double[3];  // Apex of diurnal bulge
        u_[0] = c_dec * Math.cos(ra_Sun + ra_lag);
        u_[1] = c_dec * Math.sin(ra_Sun + ra_lag);
        u_[2] = Math.sin(dec_Sun);
        VectorN u = new VectorN(u_);
        //* Cosine of half angle between satellite position vector and
        //* apex of diurnal bulge
        double c_psi2 = (0.5 + 0.5 * r.dotProduct(u)/r.mag());
        //* Density interpolation
        //double pmin = computeMin(f107,alt);	//[kg/km^3]
        //double pmax = computeMax(f107,alt); //[kg/km^3]
        double[] minmax = this.get_min_max(alt);
        double cmin0 = minmax[0];
        double cmax0 = minmax[1];
        double cmin1 = minmax[2];
        double cmax1 = minmax[3];
        double h0 = data[h_index][0];
        double h1 = data[h_index+1][0];
        //* Density computation
//        double pmin = minmax[0];
//        double pmax = minmax[1];
//        density = pmin+(pmax-pmin)*Math.pow(c_psi2,n_prm/2);
//        System.out.println("hp density: "+density);
        
        
        // Height index search and exponential density interpolation

        double h_min = ( h0 - h1 )/Math.log( cmin1/cmin0 );
        double h_max = ( h0 - h1 )/Math.log( cmax1/cmax0 );

        double d_min = cmin0 * Math.exp( (h0-alt)/h_min );
        double d_max = cmax0 * Math.exp( (h0-alt)/h_max );

        // Density computation

        density = d_min + (d_max-d_min)*Math.pow(c_psi2,n_prm/2);
        return density * 1.0e-9;

        //      * Convert units and return density
//        double rho1 = 0; //-0.0045; //-0.052;
//        return density * 1.0e-9 * (1+rho1);       // [kg/m^3]	    
	}

	
	public static void main(String[] args) {
	    VectorN r = new VectorN(-4453783.586,-5038203.756,-426384.456);
	    VectorN v = new VectorN(3831.888,-2887.221,-6.018232);
	    Spacecraft sc = new Spacecraft(r,v,0,2.2,20,1000);
		HarrisPriester atmos = new HarrisPriester(sc);
//		EarthRef ref = new EarthRef(53157.5);
//		double rmag = ref.get_JPL_Sun_Vector().mag();
//		atmos.setF107(Math.pow(149597870.0/rmag,2)*150);
//		double d = atmos.computeDensity(ref,r);
//		double tmp = d;
//		atmos.compute(ref,r,v);
//		System.out.println("* new model *   f10.7: "+atmos.f107);
//		System.out.println("Density: "+d);
//		System.out.println("Accel:   "+atmos.dragAccel());
	}

}
