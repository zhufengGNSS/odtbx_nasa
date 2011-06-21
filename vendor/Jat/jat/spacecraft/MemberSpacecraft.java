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
import jat.matvec.data.Quaternion;
import jat.timeRef.RSW_Frame;

/**
 * Interface representing a member spacecraft of a formation.
 * 
 * @author Richard C. Page III
 *
 */
public interface MemberSpacecraft {
    
    /**
     * Get the relative position (inertial)
     * @param ps Primary Spacecraft
     * @return relative position [m]
     */
    public VectorN get_rel_pos(PrimarySpacecraft ps);
    /**
     * Get the relative velocity (inertial)
     * @param ps Primary Spacecraft
     * @return relative position [m/s]
     */
    public VectorN get_rel_vel(PrimarySpacecraft ps);
    
    /**
     * Return the state vector containing the relative state.
     * @param ps Primary Spacecraft
     * @return Relative State
     */
    public double[] get_rel_state(PrimarySpacecraft ps);
    
    /**
     * Get the relative position to another member spacecraft.
     * @param ms Member spacecraft
     * @param ps Primary spacecraft
     * @return position [m]
     */
    //* get relative position of MemberSpacecraft ms
    public VectorN get_rel_pos(MemberSpacecraft ms, PrimarySpacecraft ps);
    /**
     * Get the relative velocity (inertial) to another member spacecraft.
     * @param ms Member spacecraft
     * @param ps Primary spacecraft
     * @return velocity [m/s]
     */
    //* get relative velocity of MemberSpacecraft ms
    public VectorN get_rel_vel(MemberSpacecraft ms, PrimarySpacecraft ps);
    
    
    /**
     * Get the attitude
     * @return quaternion
     */
    public Quaternion get_attitude();
    /**
     * Send a control message with the distance between this and another member spacecraft
     * @param target distance
     * @param s Member Spacecraft
     */
    public void send_control(double distance, MemberSpacecraft s);
    /**
     * Send a control message with a vector from the primary spacecraft
     * @param pos target relative position vector
     * @param s Primary spacecraft
     */
    public void send_control(VectorN pos, PrimarySpacecraft s);
    
}
