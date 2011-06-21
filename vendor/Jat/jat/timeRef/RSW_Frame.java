package jat.timeRef;

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
 * File Created on Aug 4, 2003
 */
 
import java.io.Serializable;
import jat.matvec.data.*;
/**
 * <P>
 * The RSW_Frame Class provides methods of converting between the ECI and RSW frame,
 * accounting for the Coriolis effect. The RSW frame is also referred to as:
 * RTN, UVW, LVLH. The first axis is along the radius vector, the third axis is along
 * the angular momentum vector and the second is along track.
 * Reference: Vallado.
 * 
 * @deprecated
 * @see jat.spacetime.RSW_Frame
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class RSW_Frame implements Serializable {

	private static final long serialVersionUID = 1L;
	
	/** Reference ECI position vector */
	private VectorN rref;
	
	/** Reference ECI velocity vector */
	private VectorN vref;
	
	/** ECI to RSW direction cosine matrix */
	private Matrix eci2rsw;
	
	/** The angular velocity of RSW frame wrt ECI frame, expressed in RSW frame */
	private VectorN omega;
	
	/** Construct an RSW reference frame
	 * @param r reference ECI position vector
	 * @param v reference ECI velocity vector
	 */		
	public RSW_Frame(VectorN r, VectorN v){
		this.rref = r;
		this.vref = v;
		
		this.eci2rsw = this.ECI2RSW(r, v);
	}
	
	/** Return the ECI to RSW direction cosine matrix
	 * @return the ECI to RSW direction cosine matrix
	 */	
	public Matrix ECI2RSW(){
		return this.eci2rsw;
	}
	
	public VectorN omega(){
		return this.omega;
	}
	
	/** Transform a position vector from ECI to RSW frame.
	 * @param r_eci position vector in ECI frame
	 * @return VectorN containing RSW position vector.
	 */			
	public VectorN transform(VectorN r_eci) {
		VectorN r_rsw = this.eci2rsw.times(r_eci);
		return r_rsw;
	}
	
	/** Transform a position and velocity vector pair from ECI to RSW frame,
	 * accounting for Coriolis.
	 * @param r_eci position vector in ECI frame
	 * @param v_eci velocity vector in ECI frame
	 * @return VectorN containing RSW position and velocity vector.
	 */			
	public VectorN transform(VectorN r_eci, VectorN v_eci){
		VectorN r_rsw = this.eci2rsw.times(r_eci);
		VectorN v_rsw = this.transformVelocity(r_rsw, v_eci);
		VectorN out = new VectorN(r_rsw, v_rsw); 
		return out;
	}
		
	/** Compute the transformation from ECI to the RSW frame.
	 * @param rref target position vector in ECI frame.
	 * @param vref target velocity vector in ECI frame.
	 * @return the ECI to RSW transformation matrix
	 */
	private Matrix ECI2RSW(VectorN rref, VectorN vref) {
		VectorN h = rref.crossProduct(vref);
		VectorN rhat = rref.unitVector();
		VectorN what = h.unitVector();
		VectorN s = what.crossProduct(rhat);
		VectorN shat = s.unitVector();
		Matrix out = new Matrix(3, 3);
		out.setColumn(0, rhat);
		out.setColumn(1, shat);
		out.setColumn(2, what);
		
		double rmag = rref.mag();
		double num = shat.dotProduct(vref);
		double w = num/rmag;
		this.omega = new VectorN(0.0, 0.0, w);
		
		return out.transpose();
	}
	
	/**
	 * Utility function to transform between inertial and radial-intrack-crosstrack
	 * @param rref Position (inertial)
	 * @param vref Velocity (inertial)
	 * @return Transformatin Matrix
	 */
	public static Matrix ECI2RIC(VectorN rref, VectorN vref){
	    VectorN h = rref.crossProduct(vref);
		VectorN rhat = rref.unitVector();
		VectorN what = h.unitVector();
		VectorN s = what.crossProduct(rhat);
		VectorN shat = s.unitVector();
		Matrix out = new Matrix(3, 3);
		out.setColumn(0, rhat);
		out.setColumn(1, shat);
		out.setColumn(2, what);
				
		return out.transpose();
	}
	
	/** Compute the transformation from ECI to the SWR frame.
	 * For use with a sun tracking nadir pointing trajectory
	 * and calculating the required jaw
	 * @param rref target position vector in ECI frame.
	 * @param vref target velocity vector in ECI frame.
	 * @return the ECI to SWR transformation matrix
	 */
	public static Matrix ECI2SWR(VectorN rref, VectorN vref) {
		
		VectorN h = rref.crossProduct(vref);
		h = h.times(-1.0);
		VectorN rhat = rref.unitVector();
		rhat = rhat.times(-1.0);
		VectorN what = h.unitVector();
		VectorN s = what.crossProduct(rhat);
		VectorN shat = s.unitVector();
		Matrix out = new Matrix(3, 3);
		out.setColumn(0, rhat);
		out.setColumn(1, shat);
		out.setColumn(2, what);
		
		
		
		return out.transpose();
	}
	
	
	/** Transform a velocity vector from ECI to the RSW frame, accounting
	 * for Coriolis.
	 * @param r_rsw position vector in RSW frame.
	 * @param v_eci velocity vector in ECI frame.
	 * @return the velocity in the RSW frame
	 */	
	private VectorN transformVelocity(VectorN r_rsw, VectorN v_eci){
		VectorN v = this.eci2rsw.times(v_eci);
		VectorN cor = this.omega.crossProduct(r_rsw);
		VectorN out = v.minus(cor);
		return out;		
	}	
	

}
