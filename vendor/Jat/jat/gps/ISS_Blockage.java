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
 * The ISS_Blockage Class provides a model of GPS signal blockage due to 
 * a spherical ISS.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class ISS_Blockage implements Visible {
	
	private double elevationMask;
	private static final double issRadius = 50.0;
	private ElevationMask elev;
	
	public ISS_Blockage() {
		this.elevationMask = 10.0 * MathUtils.DEG2RAD;
		this.elev = new ElevationMask();
	}
	
	/** Constructor
	 * @param em elevation mask angle in degrees
	 */
	public ISS_Blockage(double em) {
		this.elevationMask = em * MathUtils.DEG2RAD;
		this.elev = new ElevationMask(this.elevationMask);
		
	}
	
    /**
     * Determine if the GPS satellite is visible, including ISS blockage
     * Used by GPS Measurement Generator.
     * @param losu Line of sight unit vector
     * @param r    current position vector
     * @param rISS current position vector of the ISS
     * @return boolean true if the GPS satellite is visible
     */
	public boolean visible(VectorN losu, VectorN r, VectorN rISS) {
		
		// check elevation mask
		boolean visible = elev.visible(losu, r, rISS);
		
		// check ISS visibility
		VectorN rrel = rISS.minus(r);
		double dist = rrel.mag();
		double cone = Math.atan2(issRadius,dist);
		VectorN rel = rrel.unitVector();
		double cos_delta = rel.dotProduct(losu);
		double delta = Math.acos(cos_delta);
		if (delta < cone) {
			visible = false;
		}
		
		return visible;
	}	

}
