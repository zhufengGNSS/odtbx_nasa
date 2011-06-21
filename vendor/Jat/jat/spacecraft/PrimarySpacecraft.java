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
package jat.spacecraft;

import jat.matvec.data.VectorN;
import jat.matvec.data.Matrix;
import jat.matvec.data.Quaternion;
/**
 * Interface representing the primary or central spacecraft of a formation.
 * 
 * @author Richard C. Page III
 *
 */
public interface PrimarySpacecraft {

    /**
     * Get absolute position
     * @return position [m]
     */
    public VectorN get_abs_pos();
    /**
     * Get absolute velocity
     * @return velocity [m/s]
     */
    public VectorN get_abs_vel();
    /**
     * Get attitude
     * @return quaternion
     */
    public Quaternion get_attitude();
    
    /**
     * Get the rotation matrix from inertial to radial-intrack-crosstrack
     * @return Transformation Matrix
     */
    //* Get the local Radial/In-track/Cross-track frame
    public Matrix get_inertial2RIC();
    /**
     * Get the vector angular rotation of the radial-intrack-crosstrack frame
     * @return omega vector
     */
    //* Get the rotation rate of the RIC frame
    public VectorN get_omegaRIC();
    
}
