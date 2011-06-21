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

package jat.eph;

/**
 * Enumeration of the DE405 result types.  See DE405 for a complete
 * description of each term.
 * @see DE405
 */
public enum DE405_Body {	
	
	// Developer warning: do not add, remove, or reorder these without making
	// changes in DE405!
	
		EARTH,
		MERCURY,
	    VENUS,
	    EM_BARY,  // Earth-Moon Barycenter
	    MARS,
	    JUPITER,
	    SATURN,
	    URANUS,
		NEPTUNE,
		PLUTO,
		GEOCENTRIC_MOON,
		SUN,
		MOON_LIB, // Moon Librations
		GEOCENTRIC_SUN,
		SOLAR_SYSTEM_BARY,
		MOON}