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

package jat.test;

import junit.framework.Test;
import junit.framework.TestSuite;

public final class JatUnitTest {

	// /////////////////////////////////////////////////////////////////////
	// Main operation

	/**
	 * This is the method called to start the test application.
	 * 
	 * @param args
	 */
	public static void main(final String[] args) {
		// junit.textui.TestRunner.run(JatUnitTest.class);
	}

	// /////////////////////////////////////////////////////////////////////
	// Test Operations

	/**
	 * Launch all the tests of the suite.
	 * 
	 * @return Test
	 */
	public static Test suite() {
		TestSuite suite = new TestSuite("JAT Unit Test Suite");
		suite.addTestSuite(jat.eph.unittest.DE405Test.class);
		suite.addTestSuite(jat.forces.unittest.CIRAExponentialDragTest.class);
		suite.addTestSuite(jat.forces.unittest.GravityModelTest.class);
		suite.addTestSuite(jat.forces.unittest.HarrisPriesterTest.class);
		suite.addTestSuite(jat.forces.unittest.NRLMSISEDragTest.class);
		suite.addTestSuite(jat.forces.unittest.SolarRadiationPressureTest.class);
		suite.addTestSuite(jat.groundstations.test.GroundStation_Test.class);
		suite.addTestSuite(jat.ground_tracking.unittest.ChargedParticleModelTest.class);
		suite.addTestSuite(jat.matlabInterface.unittest.JATIntegratorsTest.class);
		suite.addTestSuite(jat.matlabInterface.unittest.BatchCalcsTest.class);
		suite.addTestSuite(jat.spacetime.unittest.BodyCenteredInertialRefTest.class);
		suite.addTestSuite(jat.spacetime.unittest.EarthFixedRefTest.class);
		suite.addTestSuite(jat.spacetime.unittest.LunaFixedRefTest.class);
		suite.addTestSuite(jat.spacetime.unittest.FitIERSTest.class);
		suite.addTestSuite(jat.spacetime.unittest.LeapSecondTableTest.class);
		suite.addTestSuite(jat.spacetime.unittest.TimeUtilsUTCTest.class);
		suite.addTestSuite(jat.spacetime.unittest.TimeUtilsLeapYearTest.class);
		
		return suite;
	}

}
