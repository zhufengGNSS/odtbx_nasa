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
package jat.matlabInterface;

import com.mathworks.jmi.Matlab;

/**
 * @author Richard C. Page III
 *
 */
public class MatlabFunc {

    private String cmd;
    
    public MatlabFunc(String func){
        cmd = func;
    }
    
    public double[] call(Object[] inputArgs){
        double[] returnVals = null;
		try {
			returnVals = (double[])Matlab.mtFevalConsoleOutput(cmd, inputArgs, 0);
//			returnVals = (double[])Matlab.mtFeval(cmd, inputArgs, 0);
		} catch (Exception e) {
			e.printStackTrace();
		}		
		return returnVals;
    }
    
    public static void main(String[] args) {
        Double test = new Double(1.0);
        double[] out = null;
        MatlabFunc test_func = new MatlabFunc("test_func");
        Double[] input = new Double[2];
        input[0] = test;
        input[1] = new Double(42);
        out = test_func.call(input);
        
    }
}
