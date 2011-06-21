/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2002 United States Government as represented by the
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

package jat.forces.drag;
import jat.cm.Constants;
import jat.forces.ForceModel;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;
import jat.spacecraft.Spacecraft;
import jat.spacetime.BodyRef;
import jat.spacetime.EarthFixedRef;
import jat.spacetime.EarthRef;
import jat.spacetime.ReferenceFrameTranslater;
import jat.spacetime.Time;

/**
 * <P>
 * The AtmosphericDrag class computes the acceleration due to drag on a satellite
 * using an Earth atmosphere model that conforms to the computeDensity abstract method.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

abstract public class AtmosphericDrag implements ForceModel {

	private static final long serialVersionUID = 1L;

	/**
	 * Drag coefficient
	 */
	protected double cd;
	/**
	 * Cross sectional area [m^2]
	 */
	protected double area;
	/**
	 * Spacecraft mass [kg]
	 */
	protected double mass;
	/**
	 * Atmospheric drag [m/s/s]
	 */
	protected VectorN drag;
	/**
	 * dadv
	 */
	protected Matrix dadv;
	/**
	 * dadcd
	 */
	protected VectorN dadcd;
	/**
	 * Rotation rate of the Earth (taken from Constants.omega_e)
	 */
	public double omega_e = Constants.omega_e;//7.2921157746E-05;  // earth rotation rate - [rad/s]

	/**
	 * Constructor
	 * @param cd coefficient of drag
	 * @param area drag cross-sectional area
	 * @param mass mass
	 */	
	public AtmosphericDrag(double cd, double area, double mass){
		this.cd = cd;
		this.area = area;
		this.mass = mass;
	}

	/**
	 * Constructor using the Spacecraft class for parameters
	 * @param sc Spacecraft parameters
	 */
	public AtmosphericDrag(Spacecraft sc){
		this.cd = sc.cd();
		this.area = sc.dragArea();
		this.mass = sc.mass();
	}

	/** Abstract class requires the subclass to compute the atmospheric density.
	 * @param ref EarthRef object.
	 * @param r Position vector.
	 * @return Atmospheric density in kg/m^3
	 */
	abstract public double computeDensity(Time t, BodyRef ref, VectorN r);

//	/** Abstract class requires the subclass to compute the atmospheric density.
//	* @deprecated
//	* @param ref EarthRef object.
//	* @param r Position vector.
//	* @return Atmospheric density in kg/m^3
//	*/
//	abstract public double computeDensity(EarthRef ref, VectorN r);


	/**
	 * Return the acceleration due to drag
	 * @return VectorN containing the acceleration due to drag
	 */
	public VectorN dragAccel(){
		return this.drag;
	}

	/**
	 * Return the partial derivative of acceleration wrt velocity
	 * @return Matrix containing the partial derivative of acceleration wrt velocity
	 */
	public Matrix partialV(){
		return this.dadv;
	}

	/**
	 * Return the partial derivative of acceleration wrt Cd
	 * @return Matrix containing the partial derivative of acceleration wrt Cd
	 */
	public VectorN partialCd(){
		return this.dadcd;
	}

	/**
	 * Update the mass
	 * @param mass the new mass
	 */
	public void updateMass(double mass){
		this.mass = mass;
	}

	/** Computes the acceleration due to drag in m/s^2.
	 * @param ref EarthRef object.
	 * @param r ECI position vector in meters.
	 * @param v ECI velocity vector in meters.
	 * @return acceleration due to drag in m/s^2.
	 */
	public void compute(Time t, BodyRef ref, VectorN r, VectorN v){

		r.checkVectorDimensions(3);
		v.checkVectorDimensions(3);
		double beta = cd * area / mass;  // [m^2/kg]

		// compute the atmospheric density
		double rho = computeDensity(t, ref, r);	// [kg/m^3]
		// form -1/2 (Cd*A/m) rho
		double coeff = -0.5 * beta * rho;

		// compute the relative velocity vector and magnitude
		EarthRef eciRef = null;
		try {
			eciRef = (EarthRef)ref;
		}
		catch (ClassCastException e) {
			eciRef = new EarthRef(t);
		}
		// Two ways to do this.  Let's do both and make sure they come out
		// the same.
		VectorN v_ecef = null;
		ReferenceFrameTranslater xlater = null;
		boolean USE_METHOD_1 = false;
		if (USE_METHOD_1) {
			// Method 1: Just use eci2ecef methods of EarthRef
			VectorN rv_ecef = eciRef.eci2ecf(r, v, t);
			v_ecef = rv_ecef.get(3, 3);
		}
		else {
			// Method 2: Use the reference frame translater method.
			EarthFixedRef ecefRef = new EarthFixedRef();
			xlater = eciRef.getTranslater(ecefRef, t);
			v_ecef = xlater.translateVelocity(v, r);
		}
		double v_mag = v_ecef.mag();
		VectorN drag_ecef = v_ecef.times(v_mag * coeff);

		// Now translate ECEF drag back to ECI.
		// Again, two ways to do this.
		if (USE_METHOD_1) {
			// Method 1
			Matrix ecef2eci = eciRef.ECI2ECEF().transpose();
			this.drag = ecef2eci.times(drag_ecef);
		}
		else {
			// Method 2
			this.drag = xlater.transformDirectionBack(drag_ecef);
		}


		// form partial of drag wrt v
		Matrix vrvrt = v_ecef.outerProduct(v_ecef);
		vrvrt = vrvrt.divide(v_mag);
		Matrix vrm = new Matrix(3);
		vrm = vrm.times(v_mag);
		this.dadv = (vrvrt.plus(vrm)).times(coeff);

		// form partial of drag wrt cd
		double coeff2 = coeff / this.cd;
		this.dadcd = v_ecef.times(coeff2 * v_mag);

	}

//	/** Computes the acceleration due to drag in m/s^2.
//	* @deprecated
//	* @param ref EarthRef object.
//	* @param r ECI position vector in meters.
//	* @param v ECI velocity vector in meters.
//	* @return acceleration due to drag in m/s^2.
//	*/
//	public void compute(EarthRef ref, VectorN r, VectorN v){

//	r.checkVectorDimensions(3);
//	v.checkVectorDimensions(3);
//	double rmag = r.mag();
//	double beta = cd * area / mass;  // [m^2/kg]

//	// compute the atmospheric density
//	double rho = computeDensity(ref, r);	// [kg/m^3]

//	// compute the relative velocity vector and magnitude
//	//Matrix NP = ref.trueOfDate();
//	//VectorN we = (new VectorN(NP.get(2,0), NP.get(2,1), NP.get(2,2))).times(omega_e);
//	VectorN we = new VectorN(0, 0, omega_e);
//	VectorN wxr = we.crossProduct(r);
//	VectorN vr = v.minus(wxr);
//	double vrmag = vr.mag();

//	// form -1/2 (Cd*A/m) rho
//	double coeff = -0.5 * beta * rho;
//	double coeff2 = coeff * vrmag;

//	// compute the acceleration in ECI frame (km/s^2)
//	this.drag = vr.times(coeff2);

//	// form partial of drag wrt v
//	Matrix vrvrt = vr.outerProduct(vr);
//	vrvrt = vrvrt.divide(vrmag);
//	Matrix vrm = new Matrix(3);
//	vrm = vrm.times(vrmag);
//	this.dadv = (vrvrt.plus(vrm)).times(coeff);

//	// form partial of drag wrt cd
//	double coeff3 = coeff2 / this.cd;
//	this.dadcd = vr.times(coeff3);

//	}

	/** Implemented from the ForceModel interface
	 * @param t Time reference object
	 * @param bRef Earth reference object
	 * @param sc Spacecraft parameters and state
	 * @see jat.forces.ForceModel#acceleration(jat.spacetime.Time, jat.spacetime.BodyRef, jat.spacecraft.Spacecraft)
	 */
	public VectorN acceleration(Time t, BodyRef bRef, Spacecraft sc) {
		this.omega_e = bRef.get_spin_rate(t);
		compute(t, bRef, sc.r(), sc.v());
		return dragAccel();
	}
}



