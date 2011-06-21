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

package jat.forces.unittest;

import java.io.IOException;

import jat.forces.SolarRadiationPressure;
import jat.spacecraft.Spacecraft;

public class SolarRadiationPressureTest extends ForceModelTest {

  public static void main(String[] args) {
    junit.textui.TestRunner.run(SolarRadiationPressureTest.class);
  }

  /*
   * Test method for 'jat.forces.SolarRadiationPressure.acceleration(Time, BodyRef, Spacecraft)'
   */
  public void testAccelerationTimeBodyRefSpacecraft() throws IOException {
    Spacecraft sc = new Spacecraft();
    sc.set_area(20);
    sc.set_mass(1000);
    sc.set_cr(0.07);
    SolarRadiationPressure force = new SolarRadiationPressure(sc); 
    
    testForceModelAcceleration(sc, force, "solar_pressure.txt", 
        "solar pressure radiation");
  }
}
