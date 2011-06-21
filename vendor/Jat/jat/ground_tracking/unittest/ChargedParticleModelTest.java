/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2007 United States Government as represented by the
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
 */
package jat.ground_tracking.unittest;

import jat.ground_tracking.ChargedParticleModel;
import jat.matvec.data.VectorN;
import jat.spacetime.Time;
import junit.framework.TestCase;

import org.junit.Test;

public class ChargedParticleModelTest extends TestCase {

	/**
	 * A simple test for the ChargedParticleModel.
	 * Note, this isn't an absolute test, it is just useful for detecting
	 * a change in behavior.
	 */
	@Test
	public void testReturnDelay() {
		double frequency = 1.6e9;
		Time T = new Time(54626);
		VectorN satPos = new VectorN(2.056020937918350e7,1.411064433075980e7,0.160945891394200e7);

		ChargedParticleModel cpm = new ChargedParticleModel(frequency);
		double delay = cpm.returnDelay(satPos,T);
		
		double diff = Math.abs(delay - 0.0016458877075534145);

		assertTrue(diff < 1e13);
	}

}
