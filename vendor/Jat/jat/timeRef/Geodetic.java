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

package jat.timeRef;
import jat.matvec.data.*;
import jat.math.*;
import jat.cm.*;
/**
 * <P>
 * The Geodetic Class provides methods of converting between geocentric and geodetic reference systems.
 * Reference: Satellite Orbits by Montenbruck and Gill, Section 5.5. This is basically their C++ code converted to Java.
 * 
 * @deprecated
 * @see jat.spacetime.Geodetic
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class Geodetic {

    /** latitude in radians.
     */
    protected double latitude = 0.0;

    /** longitude in radians.
     */
    protected double longitude = 0.0;

    /** height above the ellipsoid in m.
     */
    protected double hae = 0.0;

    private double R_equ = Constants.R_Earth;

    private double f = Constants.f_Earth;

    /** Creates a new instance of Geodetic from latitude, longitude and altitude.
     * @param lambda longitude in radians.
     * @param phi latitude in radians.
     * @param alt altitude in meters.
     */
    public Geodetic(double lambda, double phi, double alt) {
        this.longitude = lambda;
        this.latitude = phi;
        this.hae = alt;
    }

    /** create a Geodetic from a Geocentric Position Vector.
     * @param r Geocentric position in m.
     */
    public Geodetic(VectorN r){
        double  eps = 1000.0 * MathUtils.MACHEPS;   // Convergence criterion
        //double  eps = MathUtils.MACHEPS;   // Convergence criterion
        double  epsRequ = eps*R_equ;
        double  e2      = f*(2.0-f);        // Square of eccentricity

        double  X = r.x[0];                   // Cartesian coordinates
        double  Y = r.x[1];
        double  Z = r.x[2];
        double  rho2 = X*X + Y*Y;           // Square of distance from z-axis

        // Iteration
        double  SinPhi;
        double  ZdZ = 0.0;
        double Nh = 0.0;
        double N = 0.0;
        double dZ = e2*Z;
        double dZ_new = 0.0;

        if (r.mag() > 0.0){

            while (Math.abs(dZ-dZ_new) > epsRequ ) {
                ZdZ    =  Z + dZ;
                Nh     =  Math.sqrt ( rho2 + ZdZ*ZdZ );
                SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
                N      =  R_equ / Math.sqrt(1.0-e2*SinPhi*SinPhi);
                dZ = dZ_new;
                dZ_new =  N*e2*SinPhi;
            }

            // Longitude, latitude, altitude
            this.longitude = Math.atan2 ( Y, X );
            this.latitude = Math.atan2 ( ZdZ, Math.sqrt(rho2) );
            this.hae = Nh - N;
        }
        else{
            this.longitude = 0.0;
            this.latitude = 0.0;
            this.hae = -Constants.R_Earth;
        }
    }


    /** Return the latitude.
     * @return latitude in radians.
     */
    public double getLatitude(){
        return this.latitude;
    }
    /** return the longitude.
     * @return longitude in radians.
     */
    public double getLongitude(){
        return this.longitude;
    }
    /** return the height above the ellipsoid.
     * @return height above the ellipsoid in m.
     */
    public double getHAE(){
        return this.hae;
    }
    /** set the equatorial radius.
     * @param r_eq equatorial radius in m.
     */
    public void setEquatorialRadius(double r_eq){
        this.R_equ = r_eq;
    }

    /** get the equatorial radius.
     * @return equatorial radius in m.
     */
    public double getEquatorialRadius(){
        return this.R_equ;
    }

    /** set the flattening factor.
     * @param fin the flattening factor.
     */
    public void setFlattening(double fin){
        this.f = fin;
    }

    /** get the flattening factor.
     * @return flattening factor.
     */
    public double getFlattening(){
        return this.f;
    }


    /** computes the geocentric position vector.
     * @return geocentric position in m.
     */
    public VectorN geocentricPosition () {
        double  e2     = f*(2.0-f);        // Square of eccentricity
        double  CosLat = Math.cos(this.latitude);         // (Co)sine of geodetic latitude
        double  SinLat = Math.sin(this.latitude);
        double  CosLon = Math.cos(this.longitude);         // (Co)sine of geodetic latitude
        double  SinLon = Math.sin(this.longitude);
        double  N;
        VectorN  r = new VectorN(3);

        // Position vector

        N = this.R_equ / Math.sqrt(1.0-e2*SinLat*SinLat);
        double Nh = N + this.hae;
        r.x[0] =  Nh*CosLat*CosLon;
        r.x[1] =  Nh*CosLat*SinLon;
        r.x[2] =  ((1.0-e2)*Nh)*SinLat;
        return r;
    }

    /** Compute the transformation from ECEF to SEZ (topocentric horizon) reference frame.
     * @return ECEF to SEZ transformation matrix
     */
    public Matrix ECEF2SEZ (){
        double lambda = this.longitude;
        double phi = this.latitude;

        RotationMatrix M = new RotationMatrix( 3, lambda, 2, (Constants.pi/2.0 - phi));
        Matrix out = new Matrix(M);
        return  out;
    }


}
