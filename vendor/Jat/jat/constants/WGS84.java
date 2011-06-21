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

package jat.constants;

public interface WGS84 {
	
    /** Mean Earth Radius in m from WGS-84.
     */
    public final static double R_Earth = 6378.137e3;      // Radius Earth [m]; WGS-84
    
    /** Flattening factor of earth from WGS-84
     */
    public final static double f_Earth = 1.0/298.257223563; // Flattening; WGS-84
    
    /** Earth gravity constant in m^3/s^2 from WGS84
     */
    public final static double GM_Earth    = 398600.5e+9;    // [m^3/s^2]; WGS-84

    /** Earth's rotation rate in rad/s.
     */
    public final static double omega_Earth = 7.2921151467E-05;  // earth rotation rate

}
