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

public interface Units {
    /** Feet to Meters.
     */
    public final static double FEET2METERS = 3.048;      // Feet to Meters (exact)
    /** Meters to Feet
     */
    public final static double METERS2FEET = 1.0/3.048; // Meters to Feet
    /** Miles to Meters
     */
    public final static double MILES2METERS    = 1.609344E+03;    // Miles to Meters (exact)

    /** Meters to Miles.
     */
    public final static double METERS2MILES = 1.0/1.609344E+03;  // Meters to Miles
    /** Nautical Miles to Meters
     */
    public final static double NAUTICALMILES2METERS    = 1.852E+03;    // Nautical Miles to Meters (exact)

    /** Meters to Nautical Miles.
     */
    public final static double METERS2NAUTICALMILES = 1.0/1.852E+03;  // Meters to Nautical Miles
    /** Kilometers to Meters
     */
    public final static double KM2METERS    = 1000.0;    // Kilometers to Meters (exact)

    /** Meters to Kilometers.
     */
    public final static double METERS2KM = 1.0/1000.0;  // Meters to Kilometers
    

}
