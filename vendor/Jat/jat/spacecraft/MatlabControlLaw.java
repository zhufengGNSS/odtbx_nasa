/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2005 United States Government as represented by the
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
package jat.spacecraft;

import jat.matlabInterface.*;
import jat.matvec.data.VectorN;
/**
 * This class represents a call to a Matlab function which runs a control algorithm.  The
 * class implements ControlLaw and can be used in JAT as any other Java based control law
 * except that the propagation must occur from within the Matlab session.
 * 
 * @see jat.sim.Readme.txt
 * @author Richard C. Page III
 *
 */
public class MatlabControlLaw extends ControlLaw {
    /**
     * Matlab function
     */
    MatlabFunc controller;
    /**
     * Constructor
     * @param cmd Matlab command (string)
     */
    public MatlabControlLaw(String cmd){
        controller = new MatlabFunc(cmd);
    }
    /**
     * Given the time in seconds since simulation epoch, the state, and the relative state
     * compute the control thrust generated.
     * @param t seconds since simulation epoch
     * @param x state [x y z xdot ydot zdot ...] in meters and meters per second
     * @param xrel same as x only relative to the primary spacecraft of the formation 
     * (if there is one, otherwise not used)
     */
    public double[] compute_control(double t, double[] x, double[] xrel){
        //Object[] args = new Object[3];
        Object[] args = new Object[3];
        Double in1 = new Double(t);
        Double[] in2 = new Double[10];
        Double[] in3 = new Double[6];
//        double[] arg0 = new double[1];
//        arg0[0] = t;
//        args[0] = arg0;
//        args[1] = x;
//        args[2] = xrel;
        for(int i=0; i<10; i++){
            in2[i] = new Double(x[i]);
        }
        for(int i=0; i<6; i++){
            in3[i] = new Double(xrel[i]);
        }
        args[0] = in1;
        args[1] = in2;
        args[2] = in3;
        //System.out.println("t: "+t+"  "+in[1]+"  "+in[4]+"  "+in[7]+" "+in[11]+" "+in.length);
        double[] out = controller.call(args);
        VectorN ref = new VectorN(xrel);
        //System.out.println("t: "+t+"  "+out[0]+"  "+out[1]+"  "+out[2]+" rel: "+ref.mag());
        return out;
    }
    
}
