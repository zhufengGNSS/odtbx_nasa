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
import jat.math.*;

/**
 * <P>
 * The ElevationMask Class provides a model of GPS signal blockage due to 
 * a minimum elevation constraint or elevation mask.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class ElevationMask implements Visible {
	
	private double elevationMask;
	
	public ElevationMask() {
		this.elevationMask = 10.0 * MathUtils.DEG2RAD;
	}
		
	/** Constructor
	 * @param em elevation mask angle in degrees
	 */
	public ElevationMask(double em) {
		this.elevationMask = em * MathUtils.DEG2RAD;
	}
	
	 /**
     * Determine if the GPS satellite is visible
     * Used by GPS Measurement Generator.
     * @param losu Line of sight unit vector
     * @param r    current position vector
     * @return boolean true if the GPS satellite is visible
     */
	public boolean visible(VectorN losu, VectorN r, VectorN rISS) {
		boolean visible = true;
		double pi2 = MathUtils.PI / 2.0;
		VectorN bore = r.unitVector();
		double cos_theta = bore.dotProduct(losu);
		double theta = Math.acos(cos_theta);
		if (theta > (pi2 - this.elevationMask)) {
			visible = false;
		}
		return visible;
	}

}
