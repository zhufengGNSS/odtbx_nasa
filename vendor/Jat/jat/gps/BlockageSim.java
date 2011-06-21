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
//import jat.matvec.data.arrayTools.*;
import jat.alg.integrators.*;
//import jat.alg.*;
import java.io.*;
import jat.math.*;

/**
 * <P>
 * The BlockageSim Class provides a simple simulation of ISS blockage.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class BlockageSim {
    
    private double t0_mjd = 51969.0;  // MJD at beginning of the sim
    private double dr;               // STS 10 meters below ISS
    
    /**
     * Constructor
     * @param d distance chaser is below the ISS in meters
     */
    public BlockageSim(double d){
        this.dr = d;
    }
    
    private double mjd_time(double tsim){
        double tdays = tsim/86400.0;
        double out = this.t0_mjd + tdays;
        return out;
    }
    
    public static void main(String args[]) throws IOException {
        
        String directory = "C:\\Jat\\jat\\output\\";
        
        FileOutputStream outfile = new FileOutputStream(directory+"mp100-5-1000.txt");
        PrintWriter pw = new PrintWriter(outfile);
        
        //        LinePrinter lp = new LinePrinter(directory, "nsats.txt");
        Printable lp = new LinePrinter();
        
        
        BlockageSim x = new BlockageSim(-100.0);
        
        double pi2 = MathUtils.PI/2.0;
        
        // Set Up the GPS Constellation
        String in_directory = "C:\\Jat\\jat\\input\\";
        String filename = "rinex.n";
        String file = in_directory + filename;
        GPS_Constellation constellation = new GPS_Constellation(file);
        
        VectorN dec = new VectorN(constellation.size());
        
        // Set Up the ISS
        double mu = 3.986004418E+14;
        ISS iss = new ISS(mu, 6678000.0, 0.005, 56.0, 0.0, 0.0, 0.0);
        iss.setDiameter(100.0);
        double zmask = iss.zenith_mask(x.dr);
//        System.out.println("zenith mask = "+(zmask * x.math.RAD2DEG));
        
        // Set Up the STS with 10 degree elevation mask
        double arcs = 1000.0;   // ISS radar cross-sectional area for multipath
        int nrays = 5;
        GPS_Receiver sts = new GPS_Receiver(10.0, x.dr, arcs, nrays);
        
        // Initialize time loop
        double t = 0.0;
        double t0 = 0.0;
//        double tf = 86400.0;
        double tf = 3600.0;
        
        while (t < tf){
                        
            // what time is it in MJD?
            double t_mjd = x.mjd_time(t);
            
//            System.out.println(t);
            //            System.out.println("mjd = "+t_mjd);
            
            // propagate ISS orbit
            VectorN rvISS = iss.propagate(t0, t, lp, false);
            VectorN rISS = rvISS.get(0, 3);

            // compute ECI2RSW transformation matrix
            Matrix rsw2eci = iss.trajectory.RSW2ECI();
            Matrix eci2rsw = rsw2eci.transpose();
            
            // process the measurements
            sts.processAllInView(t, t_mjd, constellation, rISS, zmask, eci2rsw, pw);
//            sts.process12chan(t, t_mjd, constellation, rISS, zmask, eci2rsw, pw);
            
            // increment time
            t0 = t;
            t = t + 1.0;
            
        }
        
        pw.close();
        outfile.close();
        //        lp.close();
        System.out.println("Simulation Run Completed.");
        
    }
    
    
}
