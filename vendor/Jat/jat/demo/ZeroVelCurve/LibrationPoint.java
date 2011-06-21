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

package jat.demo.ZeroVelCurve;

import jat.alg.*;
import jat.matvec.data.*;

/**
 * <P>
 * 
 * @author Tobias Berthold
 * @version 1.0
 */ 

public class LibrationPoint extends NLESolver
{
    double mu;

    /** Runs the example.
     * @param args Arguments.
     */
    
    LibrationPoint(double mu, VectorN xin)
    {
        super(xin);
        this.mu=mu;
    }

    // find libration points
  	public VectorN evaluate(VectorN xin)
	{
	    double [] xout=new double[1];
	    xout[0]=Ux(mu,xin.x[0]);
	    return new VectorN(xout);
	}
    
    static double Ux(double mu,double x)
    {
        double t1=(-1.+x+mu);
        double t2=x+mu;
        double t1_cubed=t1*t1*t1;
        double t2_cubed=t2*t2*t2;
        
        double ret_val=x-mu*t1/Math.abs(t1_cubed)-(1-mu)*t2/Math.abs(t2_cubed);
        return ret_val;        
    }
}
