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

package jat.gps_ins;
//import jat.matvec.data.*;
//import jat.gps_ins.*;
import jat.ins.*;

/**
 * <P>
 * The EOM interface provides the mechanism for passing a method
 * that computes the derivatives to an integrator.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public interface Derivs {

    /** Compute the derivatives.
     * @params t    double containing time or the independent variable.
     * @params x    VectorN containing the required data.
     * @return      double [] containing the derivatives.
     */

    public double[] derivs(double t, double[] x, INS_Measurement measl_1, INS_Measurement measl, int sw);

}

