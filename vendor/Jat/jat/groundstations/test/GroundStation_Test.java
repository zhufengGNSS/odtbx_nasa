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

package jat.groundstations.test;

import jat.groundstations.GroundStation;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;
import junit.framework.TestCase;

import org.junit.Test;

/**
 * JUnit test for {@link jat.groundstations.GroundStation} <br>
 * This test GroundStation exercises the ECEF-ECI transformation of
 * GroundStation. Other tests should be added to test the rest of this class.
 * <br>
 * Significantly: the constructors are not tested yet by this class.
 * @author abrown
 */
public class GroundStation_Test extends TestCase {

	/**
	 * accuracy tolerance for direct comparison of doubles
	 */
	static final double dbltol = 1e-15;

	/**
	 * accuracy tolerance for position comparison (m)
	 */
	static final double postol = 1e-6;

	/**
	 * accuracy tolerance for velocity comparison (m/s)
	 */
	static final double veltol = 1e-5;

	/**
	 * Test object.
	 */
	GroundStation gs = null;

	/**
	 * Set up method that creates a groundstation at 0.0 latitude, 0.0
	 * longitude, at 0.0 height above the surface.
	 */
	protected void setUp() throws Exception {
		super.setUp();
		gs = new GroundStation("TestGroundStation", 0.0, 0.0, 0.0);
	}
	
	/**
	 * Test method for
	 * {@link jat.groundstations.GroundStation#getECEFPosition()}. Ensures
	 * the position for 0,0,0 lat/lon/alt is 6378137, 0, 0 meters.
	 */
	@Test
	public void testGetECEFPosition() {
		double[] recef = gs.getECEFPosition();
		assertTrue(jat.math.MathUtils.abs(recef[0] - 6378137.0) < dbltol);
		assertTrue(jat.math.MathUtils.abs(recef[1] - 0.0) < dbltol);
		assertTrue(jat.math.MathUtils.abs(recef[2] - 0.0) < dbltol);
	}

	/**
	 * Test method for
	 * {@link jat.groundstations.GroundStation#getECEFVelocity()}. Ensures the
	 * default groundstation ECEF velocity is zero.
	 */
	@Test
	public void testGetECEFVelocity() {
		double[] v = gs.getECEFVelocity();
		assertTrue(jat.math.MathUtils.abs(v[0] - 0.0) < dbltol);
		assertTrue(jat.math.MathUtils.abs(v[1] - 0.0) < dbltol);
		assertTrue(jat.math.MathUtils.abs(v[2] - 0.0) < dbltol);
	}

	/**
	 * Test method for
	 * {@link jat.groundstations.GroundStation#getECEFVelocityVector()}. Ensures
	 * the default groundstation ECEF velocity is zero.
	 */
	@Test
	public void testGetECEFVelocityVector() {
		VectorN v = gs.getECEFVelocityVector();
		assertTrue(jat.math.MathUtils.abs(v.get(0) - 0.0) < dbltol);
		assertTrue(jat.math.MathUtils.abs(v.get(1) - 0.0) < dbltol);
		assertTrue(jat.math.MathUtils.abs(v.get(2) - 0.0) < dbltol);
	}

	/**
	 * Test method for
	 * {@link jat.groundstations.GroundStation#getECIPosition()}. <br>
	 * Compares against results found from the Vallado matlab function 
	 * ecef2eci.m for converting from ITRF to ICRF.  This isn't exactly the 
	 * same transformation as implemented in GroundStation, but with the 
	 * proper selection on input parameters the two methods can yield 
	 * equivalent results. <br>
	 * The Vallado MATLAB code can be found at: http://celestrak.com/software/vallado-sw.asp
	 * <br>
	 *  The MATLAB function is called using the same recef and vecef from
	 *  the test GroundStation object, gs, in this class using the test Julian Date. <br>
	 *  <i>[reci,veci,aeci] = ecef2eci  (recef,vecef,aecef,ttt,jdut1,lod,xp,yp,eqeterms,ddpsi,ddeps);</i>
	 *  <br>
	 *  <ul>Where:
	 *  <li>recef = from gs.getECEFPosition()</li>
	 *  <li>vecef = from gs.getECEFVelocity()</li>
	 *  <li>aecef = [0; 0; 0]; %(acceleration) </li>
	 *  <li>ttt = (jd-2451545.0)/36525; % Julian centuries since J2000 epoch </li>
	 *  <li>jdut1 = jd; %(Jan 1, 2006, 0h:0m:0s) </li>
	 *  <li>lod = 0; %excess length of day (sec) </li>
	 *  <li>xp = 0; yp = 0; % polar motion coefficients (arcsec) </li>
	 *  <li>eqeterms = 0; % term for ast calc (0,2) </li>
	 *  <li>ddpsi = 0; % delta psi correction to gcrf (rad) </li>
	 *  <li>ddeps = 0; % delta eps correction to gcrf (rad) </li>
	 *  </ul> <br>
	 *  The ECEF to ECI transformation is calculated directly from the Vallado 
	 *  methods and specified for consistency.
	 */
	@Test
	public void testGetECIPosition() {
		/*
		 * The equivalent ECEF to ECI transformation matrix from the Vallado
		 * methods from sideral motion, precession, and nutation. Mecef2eci =
		 * prec*nut*st; 
		 * Mecef2eci = 
		 * -0.181033211413338 -0.983476812561068 0.000579240965509 
		 * 0.983476980970951 -0.181033218728110 0.000040214478845 
		 * 0.000065311848932 0.000576950312264 0.999999831431336
		 * 
		 * This matrix was formed using the parameters in the javadoc using the
		 * Vallado MATLAB methods: precess.m, nutation.m, and sidereal.m
		 */

		// The above matrix as a vector packed by columns.
		double dblMecef2eci[] = { -0.181033211413338, 0.983476980970951,
				0.000065311848932, -0.983476812561068, -0.181033218728110,
				0.000576950312264, 0.000579240965509, 0.000040214478845,
				0.999999831431336 };
		Matrix Mecef2eci = new Matrix(dblMecef2eci, 3);
		Matrix Meci2ecef = Mecef2eci.transpose();

		// Check against the Vallado-calculated ECI position with the above matrix.
		VectorN reci = gs.getECIPosition(Meci2ecef);
		assertTrue(jat.math.MathUtils.abs(reci.get(0) - -1154654.623944) < postol);
		assertTrue(jat.math.MathUtils.abs(reci.get(1) - 6272750.920979) < postol);
		assertTrue(jat.math.MathUtils.abs(reci.get(2) - 416.567920) < postol);
	}

	/**
	 * Test method for
	 * {@link jat.groundstations.GroundStation#getECIVelocity()}. <br>
	 * Compares against results found from the Vallado matlab function 
	 * ecef2eci.m for converting from ITRF to ICRF.  This isn't exactly the 
	 * same transformation as implemented in GroundStation, but with the 
	 * proper selection on input parameters the two methods can yield 
	 * equivalent results. <br>
	 * The Vallado MATLAB code can be found at: http://celestrak.com/software/vallado-sw.asp
	 * <br>
	 *  The MATLAB function is called using the same recef and vecef from
	 *  the test GroundStation object, gs, in this class using the test JD. <br>
	 *  <i>[reci,veci,aeci] = ecef2eci  (recef,vecef,aecef,ttt,jdut1,lod,xp,yp,eqeterms,ddpsi,ddeps);</i>
	 *  <br>
	 *  <ul>Where:
	 *  <li>recef = from gs.getECEFPosition()</li>
	 *  <li>vecef = from gs.getECEFVelocity()</li>
	 *  <li>aecef = [0; 0; 0]; %(acceleration) </li>
	 *  <li>ttt = (jd-2451545.0)/36525; % Julian centuries since J2000 epoch </li>
	 *  <li>jdut1 = jd; %(Jan 1, 2006, 0h:0m:0s) </li>
	 *  <li>lod = 0; %excess length of day (sec) </li>
	 *  <li>xp = 0; yp = 0; % polar motion coefficients (arcsec) </li>
	 *  <li>eqeterms = 0; % term for ast calc (0,2) </li>
	 *  <li>ddpsi = 0; % delta psi correction to gcrf (rad) </li>
	 *  <li>ddeps = 0; % delta eps correction to gcrf (rad) </li>
	 *  </ul> <br>
	 *  The ECEF to ECI transformation is calculated from the Vallado methods
	 *  and specified for consistency.
	 */
	@Test
	public void testGetECIVelocity() {
		/*
		 * The equivalent ECEF to ECI transformation matrix from the Vallado
		 * methods from sideral motion, precession, and nutation. Mecef2eci =
		 * prec*nut*st; 
		 * Mecef2eci = 
		 * -0.181033211413338 -0.983476812561068 0.000579240965509 
		 * 0.983476980970951 -0.181033218728110 0.000040214478845 
		 * 0.000065311848932 0.000576950312264 0.999999831431336
		 * 
		 * This matrix was formed using the parameters in the javadoc using the
		 * Vallado MATLAB methods: precess.m, nutation.m, and sidereal.m
		 */

		// The above matrix as a vector packed by columns.
		double dblMecef2eci[] = { -0.181033211413338, 0.983476980970951,
				0.000065311848932, -0.983476812561068, -0.181033218728110,
				0.000576950312264, 0.000579240965509, 0.000040214478845,
				0.999999831431336 };
		Matrix Mecef2eci = new Matrix(dblMecef2eci, 3);
		Matrix Meci2ecef = Mecef2eci.transpose();

		// Check against the Vallado-calculated ECI velocity with the above matrix.
		VectorN veci = gs.getECIVelocity(Meci2ecef);
		assertTrue(jat.math.MathUtils.abs(veci.get(0) - -457.416142) < veltol);
		assertTrue(jat.math.MathUtils.abs(veci.get(1) - -84.198748) < veltol);
		assertTrue(jat.math.MathUtils.abs(veci.get(2) - 0.268340) < veltol);
	}

	/**
	 * A single corner-case test of the azEl() method where the spacecraft
	 * and station are coincident.
	 */
	@Test
	public void testAzEl_Coincident() {
		double[] recef = gs.getECEFPosition();
		double[] azel = gs.azEl(recef);
		
		assertTrue(jat.math.MathUtils.abs(azel[0]) < dbltol);
		assertTrue(jat.math.MathUtils.abs(azel[1]) < dbltol);
	}
	
}
