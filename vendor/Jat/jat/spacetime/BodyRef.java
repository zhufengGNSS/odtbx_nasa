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
package jat.spacetime;

import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;

/**
 * Interface used as a standard for central body reference frame objects.  Classes
 * which implement BodyRef should represent the coordinate rotation, transformation,
 * and other reference functions commonly used in astrodynamics.
 * 
 * @author Richard C. Page III
 *
 */
public interface BodyRef extends ReferenceFrame {
     
    /**
     * Spin rate
     * @param t Time object
     * @return scalar spin rate [rad/s]
     */
    public double get_spin_rate(Time t);
    /**
     * Get mean radius
     * @return radius [m]
     */
    public double get_mean_radius();
    /**
     * Get gravitational constant
     * @return grav_constant [m^3/s^2]
     */
    public double get_grav_const();
    /**
     * Transformation - inertial to body frame
     * @param t Time object
     * @return Transformation matrix
     */
    public Matrix inertial_to_body(Time t);
    /**
     * Transformation - body to inertial
     * @param t Time object
     * @return Transformation matrix
     */
    public Matrix body_to_inertial(Time t);

    /**
     * Get the current JPL vector to the Sun [m]
     * @param t - Time object
     * @return Vector 3 [m]
     */
    public VectorN get_JPL_Sun_Vector(Time t);

    /**
     * Creates a translate that can translate between two reference
     * frames at a given time
     * @param other another reference frame
     * @param t the time at which translation will be done
     * @return a translating object or null if the reference frame
     * does not know how to translate to the other reference frame.
     */
    public ReferenceFrameTranslater getTranslater(ReferenceFrame other, Time t);   
    
}
