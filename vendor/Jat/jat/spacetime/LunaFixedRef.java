/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2006 United States Government as represented by the
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
 * Emergent Space Technologies
 * File created by Richard C. Page III 
 **/
package jat.spacetime;

import jat.eph.*;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;

/**
 * Represents the Lunar Centered Fixed or Selenographic Reference Frame.
 * 
 * @author Rob Antonucci
 */
public class LunaFixedRef implements ReferenceFrame {

	private static final long serialVersionUID = 1L;
  
    /**
     * Construct a LCF reference frame.
     */
    public LunaFixedRef()
    {
      // Does nothing.
    }
    
    /**
     * Compute the LCI to LCF transformation matrix.
     */
    private Matrix computeLCI2LCF(Time t) {
      
      // First compute the Euler angles
      DE405 jpl_ephem = new DE405();
      VectorN angles = jpl_ephem.get_Moon_Libration(t.mjd_tt());

      // Then compute the transformation matrix
      double alpha = angles.get(0) - Math.PI/2;
      double sina = Math.sin(alpha);
      double cosa = Math.cos(alpha);
      double delta = Math.PI/2 - angles.get(1);
      double sind = Math.sin(delta);
      double cosd = Math.cos(delta);
      double lambda = angles.get(2);
      double sinl = Math.sin(lambda);
      double cosl = Math.cos(lambda);
      Matrix xform = new Matrix(3, 3);
      xform.set(0, 0, -cosl*sina - sinl*sind*cosa);
      xform.set(0, 1, cosl*cosa - sinl*sind*sina);
      xform.set(0, 2, sinl*cosd);
      xform.set(1, 0, sinl*sina-cosl*sind*cosa);
      xform.set(1, 1, -sinl*cosa - cosl*sind*sina);
      xform.set(1, 2, cosl*cosd);
      xform.set(2, 0, cosd*cosa);
      xform.set(2, 1, cosd*sina);
      xform.set(2, 2, sind);
      
      return xform;
    }


    /**
     * Returns a translater to translate into other reference frames.
     * @param other another reference frame
     * @param t time at which translation will be done
     * @return translater object or null if does not know how
     * to translate
     */
    public ReferenceFrameTranslater getTranslater(ReferenceFrame other, Time t)
    {
      ReferenceFrameTranslater xlater = null;
      if (other instanceof LunaFixedRef) {
        // Same reference frame.  No translation needed.
        xlater = new ReferenceFrameTranslater();
      }
      else if (other instanceof BodyCenteredInertialRef) {
        xlater = getTranslater((BodyCenteredInertialRef)other, t);
      }
      else if (other instanceof EarthRef) {
        // EarthRef is just a BodyCenteredInertialRef centered on Earth
        xlater = getTranslater(new BodyCenteredInertialRef(DE405_Body.EARTH), t);
      }
      return xlater;
    }
    
    /**
     * Returns a translater to translate to LCI or ECI or any
     * other something-CI.
     * @param inertialRef an inertial reference frame
     * @param t time at which translation will be done
     * @return translater object
     */
    private ReferenceFrameTranslater 
      getTranslater(BodyCenteredInertialRef inertialRef, Time t)
    {
      // We determine the transformation matrix from LCF to LCI.
      // This can be used for transformation to any body-centered
      // inertial frame.
      Matrix lci2lcf = computeLCI2LCF(t);
      Matrix xform = lci2lcf.transpose();
      
      // Determine the position of the other body relative to the Earth.
      // Then transform it to the ECF reference frame.
      DE405 jpl_ephem = new DE405();
      // note, DE405.get_planet_posvel() returns meters and m/s:
      VectorN state1 = new VectorN(jpl_ephem.get_planet_posvel(DE405_Body.MOON, t.mjd_tt()));
      VectorN state2 = new VectorN(6);
      if (!inertialRef.getBody().equals(DE405_Body.MOON)) {
    	  state2 = jpl_ephem.get_planet_posvel(inertialRef.getBody(), t.mjd_tt());
      }
      // We difference (already in meters)
      VectorN originDiff = state2.get(0, 3).minus(state1.get(0, 3));
      VectorN bodyPos = lci2lcf.times(originDiff);
      VectorN originVel = state2.get(3, 3).minus(state1.get(3, 3));
      VectorN bodyVel = lci2lcf.times(originVel);
      // We use -omega because we want the rotation of inertial frame with respect to
      // the fixed frame.
      double[] rotation = {0, 0, -LunaRef.omega};
      ReferenceFrameTranslater xlater =
        new ReferenceFrameTranslater(xform, bodyPos, bodyVel, new VectorN(rotation));
      
      return xlater;
    }
    
}
