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
import jat.cm.*;

/**
 * <P>
 * The INSErrorState Class provides a way to initialize and handling of
 * an INS error state
 *
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class INSErrorState {

    /** Initial position vector, ECI in km.
     */
    public VectorN r;
    /** Initial velocity vector, ECI in km/s.
     */
    public VectorN v;
    /** Initial RSW2ECI quaternion.
     */
    public Quaternion q;
    /** Initial position error, ECI in km.
     */
    public VectorN dr;
    /** Initial velocity error, ECI in km/s.
     */
    public VectorN dv;
    /** Initial tilt errors about RSW axis in radians.
     */
    public Quaternion dq;

    public VectorN tilt;

    public VectorN bg;

    public VectorN ba;


    /** Creates a new instance of FreeINS */
    public INSErrorState() {
    }

    /** Creates a FreeINS from a set of initial conditions.
     * @param r Initial position vector, ECI in km.
     * @param v Initial velocity vector, ECI in km/s.
     * @param q Initial platform to ECI quaternion.
     * @param dr Initial position error, RSW in m.
     * @param dv Initial velocity error, RSW in m/s.
     * @param tilt Initial tilt errors about RSW axes in arc-sec.
     */
    public INSErrorState(VectorN r, VectorN v, Quaternion q, VectorN dr, VectorN dv, VectorN tilt){
        r.checkVectorDimensions(3);
        v.checkVectorDimensions(3);
        dr.checkVectorDimensions(3);
        dv.checkVectorDimensions(3);
        tilt.checkVectorDimensions(3);

        // set initial position and velocity
        this.r = new VectorN(r);
        this.v = new VectorN(v);
        this.q = new Quaternion(q);

        // rotate position and velocity errors from RSW to ECI and convert to km
        TwoBody orbit = new TwoBody(r, v);
        Matrix rsw2eci = orbit.RSW2ECI();
        double h = orbit.getHmag();
        double rmag = r.mag();
        double rsq = rmag * rmag;
        VectorN omega = new VectorN(0.0, 0.0, (h/rsq));
        //        rsw2eci.print();
        VectorN drtmp = rsw2eci.times(dr);
        VectorN coriolis = omega.crossProduct(dr);
        VectorN dvtmp = rsw2eci.times(dv.plus(coriolis));
        this.dr = drtmp.divide(1000.0);
        this.dv = dvtmp.divide(1000.0);

        // convert tilts from arcsec to radians
        this.tilt = tilt.times(MathUtils.ARCSEC2RAD);
        this.dq = this.q.tilt2quat(this.tilt);
        //        tilt.print("tilt out of constructor");
    }

    /** Builds the state array from the pieces.
     * @return Total State array.
     * @param dq Error quaterion.
     * @param r Position vector.
     * @param v Velocity vector.
     * @param dr Position error.
     * @param dv Velocity error.
     * @param q Quaternion error.
     */
    public double[] stateArray(VectorN r, VectorN v, Quaternion q, VectorN dr, VectorN dv, Quaternion dq, VectorN bg, VectorN ba) {
        double out[] = new double [26];
        for (int i = 0; i < 3; i++){
            out[i] = r.x[i];
            out[i+3] = v.x[i];
            out[i+6] = q.x[i];
            out[i+10] = dr.x[i];
            out[i+13] = dv.x[i];
            out[i+16] = dq.x[i];
            out[i+20] = bg.x[i];
            out[i+23] = ba.x[i];
        }
        out[9] = q.x[3];
        out[19] = dq.x[3];
        return out;
    }


}
