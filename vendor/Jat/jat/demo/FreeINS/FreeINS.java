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

package jat.demo.FreeINS;

import java.io.*;
import jat.alg.integrators.*;
import jat.cm.*;
import jat.sensors.*;
import jat.plot.*;
//import jat.util.*;
import jat.matvec.data.*;

/**
 * <P>
 * The FreeINS Class simulates an unaided INS in orbit. Currently the true orbit
 * is a non-perturbed two-body orbit. The INS error model currently assumes
 * perfect gyros, accelerometers and gravity model.
 *
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class FreeINS implements Derivatives {

    /** Initial position vector, ECI in km.
     */
    public INSErrorState x;
    public GyroErrorModel gyro;
    public AccelerometerErrorModel accel;

    /** Create a 3 plot per page plot.
     */
    public ThreePlots plots = new ThreePlots();



    /** Creates a new instance of FreeINS */
    public FreeINS() {
    }

    /** Generates the plot.
     * @param title String containing a title for the plot.
     */
    public void plotSetup(String title) {
        plots.setTitle(title);

        plots.topPlot.setXLabel("t (s)");
        plots.middlePlot.setXLabel("t (s)");
        plots.bottomPlot.setXLabel("t (s)");

        plots.topPlot.setYLabel("position (m)");
        plots.middlePlot.setYLabel("velocity (m/s)");
        plots.bottomPlot.setYLabel("tilt (arc-s)");

        plots.topPlot.addLegend(0,"radial");
        plots.topPlot.addLegend(1,"intrack");
        plots.topPlot.addLegend(2,"xtrack");
        plots.middlePlot.addLegend(0,"radial");
        plots.middlePlot.addLegend(1,"intrack");
        plots.middlePlot.addLegend(2,"xtrack");
        plots.bottomPlot.addLegend(0,"radial");
        plots.bottomPlot.addLegend(1,"intrack");
        plots.bottomPlot.addLegend(2,"xtrack");
    }
    /** Builds the vector for printing.
     * @return Print vector.
     * @param t Time.
     * @param y State array.
     */
    public VectorN printVector(double t, double[] y) {
        VectorN out = new VectorN(23);

        // strip out the incoming data
        VectorN r = new VectorN(y[0], y[1], y[2]);
        VectorN v = new VectorN(y[3], y[4], y[5]);
        Quaternion q = new Quaternion(y[6], y[7], y[8], y[9]);
        VectorN dr = new VectorN(y[10], y[11], y[12]);
        VectorN dv = new VectorN(y[13], y[14], y[15]);
        Quaternion dq = new Quaternion(y[16], y[17], y[18], y[19]);
        //        dq.print("dq in printVector");

        // construct the INS indicated orbit
        VectorN rstar = r.minus(dr);
        VectorN vstar = v.minus(dv);
        VectorN qr = q.minus(dq);
        Quaternion qref = new Quaternion(qr);
        TwoBody ref_orbit = new TwoBody(rstar, vstar);
        double href = ref_orbit.getHmag();
        double rmag = rstar.mag();
        double rsq = rmag * rmag;
        VectorN omega = new VectorN(0.0, 0.0, (href/rsq));

        // normalize qref

        qref.unitize();

        // get the tilts in arc-sec

        VectorN tilt = dq.tilts(qref);
        //        tilt.print("tilt in rad");
        tilt = tilt.divide(Constants.arcsec2rad);
        //        tilt.print("tilt out of printVector");

        // rotate INS errors from ECI to RSW and km to m

        Quaternion qprime = qref.conjugate();
        VectorN drr = qprime.transform(dr);

        VectorN coriolis = omega.crossProduct(drr);
        VectorN dvr = qprime.transform(dv);
        dvr = dvr.minus(coriolis);
        drr = drr.times(1000.0);
        dvr = dvr.times(1000.0);

        out.set(0, r);
        out.set(3, v);
        out.set(6, q);
        out.set(10, drr);
        out.set(13, dvr);
        out.set(16, dq);
        out.set(20, tilt);

        boolean first = true;
        if (t == 0.0) first = false;

        plots.topPlot.addPoint(0, t, drr.x[0], first);
        plots.topPlot.addPoint(1, t, drr.x[1], first);
        plots.topPlot.addPoint(2, t, drr.x[2], first);
        plots.middlePlot.addPoint(0, t, dvr.x[0], first);
        plots.middlePlot.addPoint(1, t, dvr.x[1], first);
        plots.middlePlot.addPoint(2, t, dvr.x[2], first);
        plots.bottomPlot.addPoint(0, t, tilt.x[0], first);
        plots.bottomPlot.addPoint(1, t, tilt.x[1], first);
        plots.bottomPlot.addPoint(2, t, tilt.x[2], first);

        return out;
    }

    /** Computes the derivatives for integration.
     * @param t Time (not used).
     * @param y State array.
     * @return Array containing the derivatives.
     */
    public double[] derivs(double t, double[] y){
        Constants c = new Constants();
        VectorN out = new VectorN(26);

        // strip out the incoming data
        VectorN r = new VectorN(y[0], y[1], y[2]);
        VectorN v = new VectorN(y[3], y[4], y[5]);
        Quaternion q = new Quaternion(y[6], y[7], y[8], y[9]);
        VectorN dr = new VectorN(y[10], y[11], y[12]);
        VectorN dv = new VectorN(y[13], y[14], y[15]);
        Quaternion dq = new Quaternion(y[16], y[17], y[18], y[19]);
        VectorN bg = new VectorN(y[20], y[21], y[22]);
        VectorN ba = new VectorN(y[23], y[24], y[25]);

        //        q.unitize();
        //        dq.unitize();

        // construct the true orbit
        TwoBody true_orbit = new TwoBody(r, v);
        VectorN g = true_orbit.local_grav();
        Matrix omega_true = new Matrix(4, 4);
        double rmag = r.mag();
        double rsq = rmag * rmag;
        double htrue = true_orbit.getHmag();
        double wtrue = 0.5 * (htrue/rsq);
        omega_true.set(0, 1, wtrue);
        omega_true.set(1, 0, -wtrue);
        omega_true.set(2, 3, wtrue);
        omega_true.set(3, 2, -wtrue);
        VectorN q_deriv = omega_true.times(q);

        // construct the INS indicated orbit

        VectorN rstar = r.minus(dr);
        VectorN vstar = v.minus(dv);

        // compute quantities for derivs evaluated on the INS indicated orbit

        TwoBody ref_orbit = new TwoBody(rstar, vstar);
        Matrix chat = ref_orbit.RSW2ECI();
        Matrix gg = ref_orbit.gravityGradient();
        VectorN dv_grav = gg.times(dr);
        Matrix omega = new Matrix(4,4);
        rmag = rstar.mag();
        rsq = rmag * rmag;
        double href = ref_orbit.getHmag();  // reference orbit angular momentum
        double w = href/rsq;                // angular velocity
        double whalf = 0.5 * w;
        omega.set(0, 1, whalf);
        omega.set(1, 0, -whalf);
        omega.set(2, 3, whalf);
        omega.set(3, 2, -whalf);
        VectorN dq_unforced = omega.times(dq);

        // add in gyro error effects
        VectorN qr = q.minus(dq);
        Quaternion qref = new Quaternion(qr);
        VectorN wref = new VectorN(0.0, 0.0, w);
        VectorN dq_gyro = this.gyro.computeErrors(qref, wref, bg);
        VectorN dq_deriv = dq_unforced.plus(dq_gyro);
        this.gyro.randwalk.nextSet();

        // add in constant acceleration
        double accel = 0.0000001 * 9.81;    // 100 micro-g
        VectorN f = new VectorN(0.0, accel, 0.0);
//        VectorN f = new VectorN(0.0, 0.0, 0.0);

        // add in accelerometer effects
        VectorN del_f = this.accel.computeErrors(f, ba);
//        del_f.print("del_f");

        // construct the measured specific force
        VectorN fhat = f.minus(del_f);

        Matrix abar = new Matrix(4, 3);
        abar.A[0][1] = -fhat.x[2];
        abar.A[0][2] = fhat.x[1];
        abar.A[1][0] = fhat.x[2];
        abar.A[1][2] = -fhat.x[0];
        abar.A[2][0] = -fhat.x[1];
        abar.A[2][1] = fhat.x[0];
        abar.A[3][0] = -fhat.x[0];
        abar.A[3][1] = -fhat.x[1];
        abar.A[3][2] = -fhat.x[2];
        Matrix atrans = abar.transpose();
        Matrix rhat = qref.rMatrix();
        Matrix rtrans = rhat.transpose();
        VectorN dv_t = atrans.times(rtrans.times(dq));
        dv_t = dv_t.times(2.0);

        // sum errors due to accelerometer and tilt
        VectorN df = del_f.plus(dv_t);

        // transform errors from body to ECI
        VectorN dv_err = chat.times(df);

        // compute total dv due to gravity and sensors
        VectorN dvdot = dv_grav.plus(dv_err);

        // derivatives for true orbit
        out.set(0, v);
        out.set(3, g);
        out.set(6, q_deriv);

        //derivatives for the INS errors
        out.set(10, dv);
        out.set(13, dvdot);
        out.set(16, dq_deriv);
        out.set(20, this.gyro.randwalk);
        out.set(23, this.accel.randwalk);
        return out.x;
    }

    /** Main program.
     * @param args Program arguments.
     * @throws IOException If there are problems with the output file.
     */
    public static void main(java.lang.String args[]) throws IOException {

        String directory = "C:\\Temp\\";

        FileOutputStream outf = new FileOutputStream(directory+"ins99.txt");
        PrintWriter pw = new PrintWriter(outf);

        FreeINS ins = new FreeINS();
        ins.plotSetup("Errors Due to Initial 1 meter Radial Position Error");

        // set initial orbit
        TwoBody orbit = new TwoBody(6770.0, 0.0, 51.0, 0.0, 0.0, 0.0);
        VectorN r0 = orbit.getR();
        VectorN v0 = orbit.getV();
        Matrix Cb2i = orbit.RSW2ECI();
        Quaternion q0 = new Quaternion(Cb2i);
        double period = orbit.period();

        // set initial INS errors

        VectorN dr0 = new VectorN(0.0, 0.0, 0.0);
        VectorN dv0 = new VectorN(0.0, 0.0, 0.0);
        VectorN tilts0 = new VectorN(0.0, 0.0, 0.0);

        ins.x = new INSErrorState(r0, v0, q0, dr0, dv0, tilts0);

        // set gyro errors

        VectorN gsf = new VectorN(0.0, 0.0, 0.0); // Vector3 containing scale factor errors in ppm.
        VectorN gma = new VectorN(6);             // Vector6 containing gyro misalignments in arcsec.
//        gma.x[1] = 20.0;
//        gma.x[3] = 20.0;
        VectorN gns = new VectorN(0.0, 0.0, 0.0); // Vector3 containing the gyro measurement noise sigmas in deg/rt-hr.
        VectorN gib = new VectorN(0.0, 0.0, 0.0); // Vector3 containing the initial gyro drift sigma in deg/hr. This is the random constant part.
        VectorN gbs = new VectorN(0.0, 0.0, 0.0); // Vector3 containing the sigmas for gyro random walk in deg/rt hr.

        ins.gyro = new GyroErrorModel(gsf, gma, gns, gib, gbs);
//        ins.gyro = new GyroErrorModel();

        // set accelerometer errors
        VectorN asf = new VectorN(0.0, 350.0, 0.0); // Vector3 containing scale factor errors in ppm.
        VectorN ama = new VectorN(6);             // Vector6 containing accel misalignments in arcsec.
//        ama.x[3] = 20.0;
        VectorN ans = new VectorN(0.0, 0.0, 0.0); // Vector3 containing the accel measurement noise sigmas in micro-g/rt-Hz.
        VectorN aib = new VectorN(0.0, 0.0, 0.0); // Vector3 containing the initial accel bias sigma in mg. This is the random constant part.
        VectorN abs = new VectorN(0.0, 0.0, 0.0); // Vector3 containing the sigmas for accel random walk in ---.

        ins.accel = new AccelerometerErrorModel(asf, ama, ans, aib, abs);
        RungeKutta8 rk8 = new RungeKutta8();
        int neqns = 26;

        double dt = rk8.step_size;
        double t = 0.0;
//        double tf = 2.0 * period;   // go 2 orbit periods
        double tf = 2000.0;   // go 2 orbit periods

        double[] newstate = new double[neqns];
        double[] oldstate = ins.x.stateArray(ins.x.r, ins.x.v, ins.x.q, ins.x.dr, ins.x.dv, ins.x.dq, ins.gyro.initialBiases, ins.accel.initialBiases);

        // check to see if the first step will go past tf

        if ((t + dt) > tf) {
            dt = tf - t;
        }

        VectorN temp = ins.printVector(t, oldstate);
        pw.print(t+"\t");
        temp.print(pw);

        // main integration loop
        //        for (int it = 0; it < 3; it++){
        while (t < tf) {
            newstate = rk8.step(t, oldstate, ins);
            for (int i = 0; i < neqns; i++) {
                oldstate[i] = newstate[i];
            }

            t = t + dt;

            Quaternion q = new Quaternion(oldstate[6], oldstate[7], oldstate[8], oldstate[9]);
            q.unitize();
            //            Quaternion dq = new Quaternion(oldstate[16], oldstate[17], oldstate[18], oldstate[19]);
            //            dq.unitize();

            for (int i = 0; i < 4; i++){
                oldstate[i+6] = q.x[i];
                //                oldstate[i+16] = dq.x[i];
            }

            temp = ins.printVector(t, oldstate);
            pw.print(t+"\t");
            temp.print(pw);
            System.out.println(t);

            if ((t + dt) > tf) {
                dt = tf - t;
            }
        }
        System.out.println("FreeINS run completed");
        ins.plots.setVisible(true);
        pw.close();
        outf.close();
        //        pw2.close();
        //        outf2.close();
    }

}
