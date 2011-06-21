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

package jat.test.propagator;


import jat.traj.*;
import jat.alg.integrators.LinePrinter;
/**
 * @author Richard C. Page III
 *
 */
public class PlotTrajectory {

    /**
     * 
     */
    public PlotTrajectory() {

    }
    
    public static void plot(String sjat, String sstk){
    	Trajectory jat = new Trajectory();
	    Trajectory stk = new Trajectory();
	    //String test = "GEO35_HP_jat.txt";
	    String file = "C:/Code/Jat/jat/test/propagator/output/"+sjat;
	    jat.readFromFile(file,"\t","km");
	    String stkfile = "C:/Code/Jat/jat/test/propagator/output/"+sstk;
	    try{
	    	stk.readFromFile(stkfile," ","km");
	    }catch(NumberFormatException e){
	    	stk = new Trajectory();
	    	stk.readFromFile(stkfile, "\t", "km");
	    }
	    LinePrinter lp = new LinePrinter();
	    RelativeTraj comp = new RelativeTraj(jat,stk,lp,sjat);
	    //comp.process(1e-5);
	    comp.process_RSS(1e-7);
    }
    
	public static void main(String[] args) throws java.io.IOException {
		//String sjat = "GEO35_HP_jat.txt";
		//String sstk = "GEO35_HP_stk.txt";
		//String sjat = "GEO34_jat.txt";
		//String sstk = "GEO34_stk.txt";
		String sjat = "GEO33_jat.txt";
		String sstk = "GEO33_stk.txt";
		//String sjat = "ISS7_HP_jat.txt";
		//String sstk = "ISS7_HP_stk.txt";
		//String sjat = "ISS4_HP_jat.txt";
		//String sstk = "ISS4_HP_stk.txt";
		PlotTrajectory.plot(sjat,sstk);
		
//		String s = "C:/Code/Jat/jat/sim/output/JAT_error_2.txt";
//		Trajectory traj = new Trajectory();
//		traj.readFromFile(s);
		
	}
}
