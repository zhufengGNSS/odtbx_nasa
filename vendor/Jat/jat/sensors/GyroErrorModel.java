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

package jat.sensors;
import jat.matvec.data.*;
import jat.math.*;

/**
 * <P>
 * The GyroErrorModel Class models a generic optical (RLG or FOG) gyro. It does
 * not include g-sensitive or non-linearity error effects.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class GyroErrorModel {
    public VectorN scaleFactors;   // random constants, just treat as constant for now
    public VectorN misalignments;  // random constants, just treat as constant for now
    public VectorN noiseSigmas;
    public VectorN initialBiases;
    public VectorN biasRWsigmas;
    public GaussianVector noise;
    public GaussianVector randwalk;
    public Matrix sf;

    /** Construct a gyro error model.
     * @param sf Vector3 containing scale factor errors in ppm.
     * @param ma Vector6 containing gyro misalignments in arcsec.
     * @param ns Vector3 containing the gyro measurement noise sigmas in deg/hr.
     * @param ib Vector3 containing the initial gyro drift sigma in deg/hr. This is the random constant part.
     * @param bs Vector3 containing the sigmas for gyro random walk in deg/rt hr.
     */
    public GyroErrorModel(VectorN sf, VectorN ma, VectorN ns, VectorN ib, VectorN bs){
        sf.checkVectorDimensions(3);
        ma.checkVectorDimensions(6);
        ns.checkVectorDimensions(3);
        ib.checkVectorDimensions(3);
        bs.checkVectorDimensions(3);

        // convert from user to internal units

        this.scaleFactors = sf.divide(1.0E+06);
        this.misalignments = ma.times(MathUtils.ARCSEC2RAD);

        double dprth2rps = MathUtils.DEG2RAD / 60.0;
        double dph2rps = dprth2rps / 60.0;
        this.noiseSigmas = ns.times(dph2rps);
        this.initialBiases = ib.times(dph2rps);
        this.biasRWsigmas = bs.times(dprth2rps);
        VectorN zeroMean = new VectorN(3);
        this.noise = new GaussianVector(zeroMean, this.noiseSigmas);
        this.randwalk = new GaussianVector(zeroMean, this.biasRWsigmas);
        this.sf = this.SFmatrix();
    }

    private Matrix SFmatrix(){
        Matrix out = new Matrix(this.scaleFactors);
        out.A[0][1] = -1.0 * this.misalignments.x[0];
        out.A[0][2] = this.misalignments.x[1];
        out.A[1][0] = this.misalignments.x[2];
        out.A[1][2] = -1.0 * this.misalignments.x[3];
        out.A[2][0] = -1.0 * this.misalignments.x[4];
        out.A[2][1] = this.misalignments.x[5];
        return out;
    }

    /** Outputs the gyro's contribution to the INS error quaternion.
     * @param qref Quaternion containing the current reference quaternion.
     * @param wref Vector3 containing the reference angular velocity vector.
     * @param bg Vector3 containing the current gyro bias state.
     * @return Error quaternion due to gyro errors.
     */
    public Quaternion computeErrors(Quaternion qref, VectorN wref, VectorN bg){
        Matrix Q = qref.qMatrix();

        // fill the gyro measurement noise vector

        this.noise.nextSet();
//        this.noise.print("noise");

        // create the scalefactor-misalignment error vector
//        wref.print("wref");
        VectorN sfma = this.sf.times(wref);
//        sfma.print("error due to sfma");
//        bg.print("gyro bias");

        // add up the effects
        VectorN sum = sfma.plus(bg.plus(this.noise));
//        sum.print("sum of all gyro errors");

        // multiply by Q matrix to get gyro contribution to error quaternion
        VectorN temp = Q.times(sum);
//        temp.print("error quaternion due to gyro");
        Quaternion out = new Quaternion(temp);
        return out;
    }
}
