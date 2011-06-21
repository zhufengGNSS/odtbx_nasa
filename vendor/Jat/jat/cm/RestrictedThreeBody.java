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

package jat.cm;
import jat.alg.integrators.*;

/**
 * <P>
 * The RestrictedThreeBody Class provides the means to generate orbits
 * in the Restricted Three Body Problem (RTBP). 
 * Ref: Vallado
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class RestrictedThreeBody implements Derivatives {
    
    private double mu;
    
    /** Creates a new instance of RestrictedThreeBody
     * @param mu_in Mass ratio.
     */
    public RestrictedThreeBody(double mu_in) {
        this.mu = mu_in;
    }
    
    /** Computes the derivatives for propagating orbits in the RTBP.
     * @param t time.
     * @param state Six element state vector.
     * @return The derivatives of the state at time t.
     */    
    public double[] derivs(double t, double[] state) {
        
        double out[] = new double[6];
        double omu = 1.0 - this.mu;
        
        double xx = state[0];
        double yy = state[1];
        double zz = state[2];
        double xdot = state[3];
        double ydot = state[4];
        double zdot = state[5];
        
        // compute state derivatives
        
        double r1 = Math.sqrt((xx + this.mu)*(xx + this.mu) + yy*yy + zz*zz);
        double r2 = Math.sqrt((xx - omu)*(xx - omu) + yy*yy + zz*zz);
        
        double r13 = r1 * r1 * r1;
        double r23 = r2 * r2 * r2;
        
        out[0] = xdot;
        out[1] = ydot;
        out[2] = zdot;
        out[3] = 2.0*ydot + xx - omu*(xx + this.mu)/r13 - mu*(xx - omu)/r23;
        out[4] = yy - 2.0*xdot - omu*yy/r13 - this.mu*yy/r23;
        out[5] = -omu*zz/r13 - this.mu*zz/r23;
        return out;
    }
    
    /** Computes the value of the Jacobi constant.
     * @param state The current 6 element state vector.
     * @return Jacobi constant.
     */    
    public double jacobiConstant(double state[]) {
        double xx = state[0];
        double yy = state[1];
        double zz = state[2];
        double xdot = state[3];
        double ydot = state[4];
        double zdot = state[5];
        double omu = 1.0 - this.mu;
        double r1 = Math.sqrt((xx + this.mu)*(xx + this.mu) + yy*yy + zz*zz);
        double r2 = Math.sqrt((xx - omu)*(xx - omu) + yy*yy + zz*zz);
        
        double out = xx*xx + yy*yy + 2.0*omu/r1 + 2.0*mu/r2 - xdot*xdot - ydot*ydot - zdot*zdot;
        return out;
    }

    /** Computes the square of the velocity. Used in computing zero-velocity curves.
     * @param state The current 6 element state vector.
     * @param Jacobi Jacobi constant.
     * @return v^2.
     */        
    public double velocitySquared(double state[], double Jacobi) {
        double xx = state[0];
        double yy = state[1];
        double zz = state[2];
        double xdot = state[3];
        double ydot = state[4];
        double zdot = state[5];
        double omu = 1.0 - this.mu;
        double r1 = Math.sqrt((xx + this.mu)*(xx + this.mu) + yy*yy + zz*zz);
        double r2 = Math.sqrt((xx - omu)*(xx - omu) + yy*yy + zz*zz);
        
        double out = xx*xx + yy*yy + 2.0*omu/r1 + 2.0*mu/r2 - Jacobi;
        return out;
    }
}
