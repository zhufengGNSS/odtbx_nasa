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
 
package jat.gps;

import jat.matvec.data.*;
import jat.cm.*;
import jat.alg.integrators.*;

/**
 * <P>
 * The ISS Class provides a model of the ISS for the BlockageSim
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class ISS {
    public TwoBody trajectory;
    private double diam;
    
    /**
     * Constructor
     * @param mu GM earth
     * @param a semi-major axis
     * @param e eccentricity
     * @param i inclination
     * @param raan right ascension of the ascending node
     * @param w argument of perigee
     * @param ta true anomaly
     */
    public ISS(double mu, double a, double e, double i, double raan, double w, double ta){
        this.trajectory = new TwoBody(mu, a, e, i, raan, w, ta);
    }
    
    /**
     * Set the diameter of the ISS
     * @param d diameter in meters
     */
    public void setDiameter(double d){
        this.diam = d;
    }
    
    /**
     * Return the diameter of the ISS
     * @return the diameter of the ISS in meters
     */
    public double diameter(){
        return this.diam;
    }
    
    /** 
     * Determine the zenith mask
     * @param dr distance from the ISS
     * @return the zenith mask angle in radians
     */
    public double zenith_mask(double dr){
        double radius = this.diam/2.0;
        double out = Math.atan2(radius, -dr);
        return out;
    }
    
    /**
     * propagate the ISS trajectory
     * @param t0 initial time in seconds
     * @param tf final time in seconds
     * @param pr Printable
     * @param print_switch boolean, true = print, false = don't print
     */
    public VectorN propagate(double t0, double tf, Printable pr, boolean print_switch){
        this.trajectory.propagate(t0, tf, pr, print_switch);
        return this.trajectory.rv;
    }
    
}