package jat.gps;

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
 * 
 * File Created on Jul 13, 2003
 */
 import jat.matvec.data.*;
 
/**
 * The Visible interface is used to provide a means for GPS SV visibility checking. 
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public interface Visible {
	/**
	 * check to see if the GPS SV is visible
	 * @param losu GPS SV line of sight unit vector
	 * @param r receiver position vector
	 * @param rISS ISS position vector
	 * @return boolean true = visible, false = blocked
	 */
	public boolean visible(VectorN losu, VectorN r, VectorN rISS);
}
