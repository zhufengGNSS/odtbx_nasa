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
 *
 */

package jat.gps;
import jat.matvec.data.*;
import jat.matvec.data.arrayTools.*;
//import jat.alg.integrators.*;
//import jat.alg.*;
import jat.math.*;
import java.io.*;

/**
 * <P>
 * The GPS_Receiver Class provides a GPS receiver model for the BlockageSim
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class GPS_Receiver {
            
    private double elmask;
    private double delr;
    private MultipathModel mp;
    
    /**
     * Constructor
     * @param el elevation mask angle in degrees
     * @param dr distance from the ISS in meters
     * @param arcs ISS radar cross-sectional area for multipath
     * @param nrays int containing number of multipath rays per direct signal
     */
    public GPS_Receiver(double el, double dr, double arcs, int nrays){
        this.elmask = el * MathUtils.DEG2RAD;
        this.delr = dr;
        this.mp = new MultipathModel(arcs, nrays);
    }
    
    /**
     * Return the elevation mask angle
     * @return the elevation mask angle in radians
     */
    public double elevationMask(){
        return this.elmask;
    }
    
    /** Compute the RSW to ECI transformation matrix.
     * @param r VectorN containing the ECI position vector
     * @param v VectorN containing the ECI velocity vector
     * @return Matrix containing the transformation.
     */
    public Matrix RSW2ECI(VectorN r, VectorN v) {
        VectorN h = this.getH(r, v);
        VectorN rhat = r.unitVector();
        VectorN what = h.unitVector();
        VectorN s = what.crossProduct(rhat);
        VectorN shat = s.unitVector();
        Matrix out = new Matrix(3,3);
        out.setColumn(0, rhat);
        out.setColumn(1, shat);
        out.setColumn(2, what);
        return out;
    }
    
    /** Compute the ECI to RSW transformation matrix.
     * @param r VectorN containing the ECI position vector
     * @param v VectorN containing the ECI velocity vector
     * @return Matrix containing the transformation.
     */
    public Matrix ECI2RSW(VectorN r, VectorN v){
        Matrix temp = this.RSW2ECI(r, v);
        Matrix out = temp.transpose();
        return out;
    }
    
    /**
     * Return the angular momentum vector
     * @param r VectorN containing the ECI position vector
     * @param v VectorN containing the ECI velocity vector
     * @return VectorN containing the angular momentum.
     */
    public VectorN getH(VectorN r, VectorN v) {
        VectorN out = r.crossProduct(v);
        return out;
    }
    
    /**
     * Return the GPS line of sight vector
     * @param rGPS VectorN containing the ECI position vector of the GPS SV
     * @param r VectorN containing the ECI position vector of the receiver
     * @return VectorN containing the vector from the receiver to the GPS SV.
     */
    public VectorN lineOfSight(VectorN rGPS, VectorN r){
        VectorN out = rGPS.minus(r);
        return out;
    }
    
    /**
     * Return the GPS line of sight boresight angle
     * @param los VectorN containing the GPS line of sight vector
     * @param r VectorN containing the ECI position vector of the receiver
     * @return double containing the angle from the boresight to the line of sight in radians
     */
    public double boresightAngle(VectorN los, VectorN r){
        VectorN bore = r.unitVector();
        VectorN rlos = los.unitVector();
        double cos_theta = bore.dotProduct(rlos);
        double theta = Math.acos(cos_theta);
        return theta;
    }
    
    /**
     * Return the PDOP
     * @param dop Matrix containing the dop matrix
     * @return double containing the value of PDOP
     */
    public double pdop(Matrix dop){
        dop.checkMatrixDimensions(4, 4);
        double d1 = dop.A[0][0];
        double d2 = dop.A[1][1];
        double d3 = dop.A[2][2];
        double out = Math.sqrt(d1 + d2 + d3);
        return out;
    }
    
    /**
     * Get the position of the GPS receiver from the ISS position vector
     * @param rISS ECI position vector of the ISS
     * @return VectorN containing the receiver ECI position vector
     */
    public VectorN getPosition(VectorN rISS){
        // get a unit vector along ISS position vector
        VectorN ISSunit = rISS.unitVector();
        
        // multiply by dr to get the change
        VectorN delta = ISSunit.times(this.delr);
        
        // obtain the STS position
        VectorN out = rISS.plus(delta);
        return out;
    }
    
    /**
     * Process all GPS SVs in View
     * @param t sim time in seconds
     * @param t_mjd MJD time
     * @param con GPS_Constellation
     * @param rISS ISS position vector
     * @param zmask zenith mask
     * @param eci2rsw Matrix containing the ECI to RSW transformation
     * @param pw PrintWriter
     */
    public void processAllInView(double t, double t_mjd, GPS_Constellation con, VectorN rISS, double zmask, Matrix eci2rsw, PrintWriter pw){
        double[][] h_array = new double[con.size()][4];
        VectorN dec = new VectorN(con.size());
        double[] elev = new double[con.size()];
        
        double [] cperr = new double[con.size()];
        
        double [] prerr = new double[con.size()];
        VectorN delx = new VectorN(4);
        
        // compute position based on ISS position
        VectorN r = this.getPosition(rISS);
        
        // iniitialize counters
        int below_horizon = 0;
        int below_elmask = 0;
        int blocked = 0;
        int visible = 0;
        
        // process each GPS satellite
        for (int isv = 0; isv < con.size(); isv++){
            
            // retrieve the SV
            GPS_SV sv = con.getSV(isv);
            
            int prn = sv.prn();
            
            // compute the GPS position vector
            VectorN rGPS = sv.rECI(t_mjd);
            
            // compute the angle between the boresight and LOS vector
            VectorN los =  this.lineOfSight(rGPS, r);
                        
            double theta = this.boresightAngle(los, r);
            dec.x[isv] = theta * MathUtils.RAD2DEG;
            
            double pi2 = MathUtils.PI/2.0;
            
            if (theta > pi2) {
                below_horizon++;
            }
            else {
                if (theta > (pi2 - this.elevationMask())){
                    below_elmask++;
                }
                else{
                    if (theta < zmask){
                        blocked++;
                    }
                    else {
                        // save off the DOP
                        VectorN losu = los.unitVector();
                        VectorN losuRSW = eci2rsw.times(losu);
                        h_array[visible][0] = losuRSW.x[0];
                        h_array[visible][1] = losuRSW.x[1];
                        h_array[visible][2] = losuRSW.x[2];
                        h_array[visible][3] = 1.0;
                        
                        // save off the elevation angle
                        elev[visible] = 90.0 - (theta * MathUtils.RAD2DEG);
                        
                        // compute the multipath error
                        mp.environment(prn, delr, theta, rGPS, rISS, r);
                        cperr[visible] = mp.carrierPhaseError();
                        prerr[visible] = mp.pseudorangeError();
                                               
                        visible++;
                    }
                }
            }
            
        }
        
        // compute the pdop
        double pdop = 1.0E+15;
        if (visible > 3) {
            Matrix h = new Matrix(h_array, visible, 4);
            Matrix htrans = h.transpose();
            Matrix hth = htrans.times(h);
            //            h.print("h");
            Matrix hth_inv = hth.inverse();
            //            hth_inv.print("hth_inv");
            pdop = this.pdop(hth_inv);
            
            // compute multipath errors
            VectorN error = new VectorN(cperr, visible);
            Matrix M = hth_inv.times(htrans);
            delx = M.times(error);
        }
        
//        System.out.println(t+"\t"+visible+"\t"+below_horizon+"\t"+below_elmask+"\t"+blocked+"\t"+pdop+"\t"+err[0]+"\t"+dec.x[0]+"\t"+err[1]+"\t"+dec.x[1]);
        System.out.println(t+"\t"+prerr[0]+"\t"+(this.mp.beta[0]/this.mp.beta0));
        pw.println(t+"\t"+visible+"\t"+below_horizon+"\t"+below_elmask+"\t"+blocked+"\t"+pdop+"\t"+cperr[0]+"\t"+elev[0]+"\t"+cperr[1]+"\t"+elev[1]+"\t"+prerr[0]+"\t"+prerr[1]);
//        delx.print("multipath error");
        //            dec.print("dec");
        //            lp.print(t, rSTS.x);
    }
    
    /**
     * Process the 12 highest elevation GPS SVs
     * @param t sim time in seconds
     * @param t_mjd MJD time
     * @param con GPS_Constellation
     * @param rISS ISS position vector
     * @param zmask zenith mask
     * @param eci2rsw Matrix containing the ECI to RSW transformation
     * @param pw PrintWriter
     */
    public void process12chan(double t, double t_mjd, GPS_Constellation con, VectorN rISS, double zmask, Matrix eci2rsw, PrintWriter pw){
        int numChannels = 12;
        double[][] h_array = new double[numChannels][4];
        VectorN dec = new VectorN(con.size());
        double [] theta = new double[con.size()];
        Matrix losMat = new Matrix(3, con.size());
        
        // compute position based on ISS position
        VectorN r = this.getPosition(rISS);
        
        // iniitialize counters
        int below_horizon = 0;
        int below_elmask = 0;
        int blocked = 0;
        int visible = 0;
        
        // process each GPS satellite
        for (int isv = 0; isv < con.size(); isv++){
            
            // retrieve the SV
            GPS_SV sv = con.getSV(isv);
            
            // compute the GPS position vector
            VectorN rGPS = sv.rECI(t_mjd);
            
            // compute the angle between the boresight and LOS vector
            VectorN los =  this.lineOfSight(rGPS, r);
            
            // store off in LOS Matrix
            losMat.setColumn(isv, los);
            
            theta[isv] = this.boresightAngle(los, r);
            dec.x[isv] = theta[isv] * MathUtils.RAD2DEG;
        }
        
        // pick out the 12 highest elevation SVs
        Sort highEl = new Sort(dec.x);
        double pi2 = MathUtils.PI/2.0;
        
        // process 12 highest elevation SVs
        for (int i = 0; i < numChannels; i++){
            int jsv = highEl.getOrder(i);
            
            if (theta[jsv] > pi2) {
                below_horizon++;
            }
            else {
                if (theta[jsv] > (pi2 - this.elevationMask())){
                    below_elmask++;
                }
                else{
                    if (theta[jsv] < zmask){
                        blocked++;
                    }
                    else {
                        // save off the DOP
                        VectorN los = losMat.getColumnVector(jsv);
                        VectorN losu = los.unitVector();
                        VectorN losuRSW = eci2rsw.times(losu);
                        h_array[visible][0] = losuRSW.x[0];
                        h_array[visible][1] = losuRSW.x[1];
                        h_array[visible][2] = losuRSW.x[2];
                        h_array[visible][3] = 1.0;
                        visible++;
                    }
                }
            }
        }
        
        
        
        
        // compute the pdop
        double pdop = 1.0E+15;
        if (visible > 3) {
            Matrix h = new Matrix(h_array, visible, 4);
            Matrix htrans = h.transpose();
            Matrix hth = htrans.times(h);
            //            h.print("h");
            Matrix hth_inv = hth.inverse();
            //            hth_inv.print("hth_inv");
            pdop = this.pdop(hth_inv);
        }
        
        System.out.println(t+"\t"+visible+"\t"+below_horizon+"\t"+below_elmask+"\t"+blocked+"\t"+pdop);
        pw.println(t+"\t"+visible+"\t"+below_horizon+"\t"+below_elmask+"\t"+blocked+"\t"+pdop);
        //            dec.print("dec");
        //            lp.print(t, rSTS.x);
    }
    
    
}
