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

public interface IERS_1996 {
    /** Speed of Light in m/s
     */
    public static final double c = 299792458.0;
    /** Obliquity of the Ecliptic at J2000 in arcsec
     */
    public static final double e0 = 84381.412;
    /** Mean Earth Radius in m
     */
    public final static double R_Earth = 6378136.49;    
    /** Earth gravity constant in m^3/s^2
     */    
    public static final double GM_Earth = 3.986004418E14;
    /** Earth's rotation rate in rad/s.
     */    
    public static final double w_Earth = 7.292115E-05;    

    public static final double J2_Earth = 1.0826359E-03;
    /** Flattening factor of earth
     */
    public final static double f_Earth = 1.0/298.25642;
    /** Earth mean equatorial gravity in m/s^2
     */
    public static final double g_Earth = 9.780327;
    /** Sun gravity constant in m^3/s^2
     */        
    public static final double GM_Sun = 1.327124E20;
    /** Moon - Earth mass ratio
     */
    public static final double Moon_Earth_Ratio = 0.0123000345;
    
    public static final double GM_Moon = GM_Earth * Moon_Earth_Ratio;
    
    
}
