/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2006 United States Government as represented by the
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
package jat.spacetime;

import java.io.IOException;

import jat.forces.*;
import jat.matvec.data.RotationMatrix;
import jat.matvec.data.VectorN;
import jat.spacecraft.Spacecraft;
import jat.traj.Trajectory;
import jat.util.Celestia;
import jat.util.FileUtil;
import jat.alg.ScalarfromArrayFunction;
import jat.alg.integrators.Derivatives;
import jat.alg.integrators.RungeKutta8;
import jat.alg.opt.DFP;
import jat.alg.opt.GradientSearch;
import jat.audio.SoundPlayer;
import jat.cm.Lambert;
import jat.eph.*;

/**
 * TODO Javadoc
 * @author Richard C. Page III
 */
public class SolarSystemModel implements Derivatives, ScalarfromArrayFunction {

	protected Time time;
	
	protected SolarBarycenterRef origin;
	
	protected DE405 ephem;
	
	public static int SUN = 0;
	public static int MERCURY = 1;
	public static int VENUS = 2;
	public static int EARTH = 3;
	public static int MARS = 4;
	public static int JUPITER = 5;
	public static int SATURN = 6;
	public static int URANUS = 7;
	public static int NEPTUNE = 8;
	public static int MOON = 9;
	
	private int[] bodies = {DE405_Body.SUN.ordinal(), DE405_Body.MERCURY.ordinal(), DE405_Body.VENUS.ordinal(), DE405_Body.EARTH.ordinal(), DE405_Body.MARS.ordinal(), DE405_Body.JUPITER.ordinal(), DE405_Body.SATURN.ordinal(), DE405_Body.URANUS.ordinal(), DE405_Body.NEPTUNE.ordinal(), DE405_Body.PLUTO.ordinal(), DE405_Body.MOON.ordinal()};
	
	/**
	 * Gravitational mass for the planets.  Source: JPL HORIZONS
	 * NOTE: Units [km^3 / s^2] 
	 */
	public static double GM_body[] = {
			1.3271243994e11, 	// sun
			22032.090,	// mercury
			324858.630,	// venus
			398600.440, // earth
			42828.3,	// mars
			126686537,	// jupiter
			37931284.5,	// saturn
			5793947,	// uranus
			6835107,	// neptune
			4902.798	// luna (Earth's moon)			
	};
	
	protected GravitationalBody body[];
	
	/**
     * List of the force models present in the simulation universe.
     */
    protected ForceModelList forces;
    
    /**
     * Epoch - initial ephemeris time.
     */
    protected double epoch;
	
    //* TEMP
    protected Lambert lam;
    protected DE405 ephem2;
    protected double travel_time = 7*29; //days
    protected EarthRef earth;
    protected RotationMatrix rot0,rotf;
    protected VectorN r0;
    
    /**
     * Creates the solar system at the given mean julian date.
     * @param mjd_utc Mean julian date in universal coordinated time.
     */
	public SolarSystemModel(double mjd_utc){
		epoch = mjd_utc;
		time = new Time(mjd_utc);
		origin = new SolarBarycenterRef();
		ephem = new DE405();
		body = new GravitationalBody[10];
		forces = new ForceModelList();
		initializeGravity(this.time);
		initializeSolarPressure(this.time);
		earth = new EarthRef(new Time(this.time.mjd_utc()+travel_time));
		update(time);
		//* TEMP
		lam = new Lambert(this.convert_GM_SI(SolarSystemModel.GM_body[SolarSystemModel.SUN]));
		ephem2 = new DE405();
    	rot0 = earth.EclMatrix(this.time.mjd_tdb());
    	time.update(travel_time*86400);
    	rotf = earth.EclMatrix(this.time.mjd_tdb());
    	time.update(0);
		r0 = this.body[SolarSystemModel.EARTH].getPosition().plus(new VectorN(100000,100000,100000));
	}
	
	public void initializeGravity(Time t){
		for(int i=0; i<10; i++){
			body[i] = new GravitationalBody(convert_GM_SI(GM_body[i]));
			forces.addForce(body[i]);
		}		
	}
	
	public void initializeSolarPressure(Time t){
		//* not doing anything yet
		//force.addForce(??)
	}
	
	public void update(Time t){
		this.time = t;
//		ephem.planetary_ephemeris(t.jd_tdb());
		VectorN r_tmp,v_tmp, tmp;

		tmp = ephem.get_planet_posvel(DE405_Body.SUN, t.mjd_tt());
		r_tmp = tmp.get(0,3);
		v_tmp = tmp.get(3,3);
		//earth.update(time.mjd_ut1(),time.mjd_tt());
		RotationMatrix rot = earth.EclMatrix(t.mjd_tt());
		r_tmp = rot.transform(r_tmp);
		v_tmp = rot.transform(v_tmp);
		body[SUN].updatePosition(r_tmp);
		body[SUN].updateVelocity(v_tmp);
		tmp = ephem.get_planet_posvel(DE405_Body.MOON, t.mjd_tt());
		r_tmp = tmp.get(0,3);
		v_tmp = tmp.get(3,3);
		r_tmp = rot.transform(r_tmp);
		v_tmp = rot.transform(v_tmp);
		body[MOON].updatePosition(r_tmp);
		body[MOON].updateVelocity(v_tmp);
		DE405_Body[] bod = DE405_Body.values();
		for(int i=1; i<8; i++){
			int j = bodies[i];        // DE405 index
			tmp = ephem.get_planet_posvel(bod[j], t.mjd_tt());
			r_tmp = tmp.get(0, 3);
			v_tmp = tmp.get(3, 3);
			r_tmp = rot.transform(r_tmp);
			v_tmp = rot.transform(v_tmp);
			body[i].updatePosition(r_tmp);
			body[i].updateVelocity(v_tmp);
		}
	}

	public void update(double t){
		time.update(t);
		update(time);
	}
	public void update_dt(double dt){
		time.step_seconds(dt);
		update(time);
	}
	
	public VectorN getBodyPosition(int b){
		return body[b].getPosition();
	}
	public VectorN getBodyVelocity(int b){
		return body[b].getVelocity();
	}
	public double get_mjd_utc(){
		return time.mjd_utc();
	}
	
	public double convert_GM_SI(double gm){
		return gm*1000000000;
	}

    /**
     * Though UniverseModel does not implement the Derivatives interface, this
     * method has a similar purpose in providing the derivatives of a spacecraft
     * state.  Given the state and relevant properties for a spacecraft, derivs
     * calculates the acceleration due to the current force list and returns
     * the derivaties.
     * 
     * This method is called by SimModel to provide the derivatives to the 
     * RungeKutta integrator.
     * 
     * @see jat.sim.SimModel
     * 
     * @param t Current simulation time in seconds since the initial epoch.
     * @param sc Spacecraft to opperate upon.
     * @return The derivatives of the spacecraft state.
     */
    public double[] derivs(double t, Spacecraft sc) {
        double[] X = sc.toStateVector();
        double[] out = new double[X.length];
        time.update(t);
        update(time);
        double[] accel = forces.acceleration(time,origin,sc);
        out[0] = X[3];
        out[1] = X[4];
        out[2] = X[5];
        out[3] = accel[0];
        out[4] = accel[1];
        out[5] = accel[2];
        for(int i=6; i<X.length; i++) out[i] = 0;
        return out;
    }
	

	/* USED FOR TESTING
	 * (non-Javadoc)
	 * @see jat.alg.integrators.Derivatives#derivs(double, double[])
	 */
	public double[] derivs(double t, double[] x) {
	    double[] out = new double[x.length];
        update(t);
        Spacecraft sc = new Spacecraft(new VectorN(x[0],x[1],x[2]),
        		                       new VectorN(x[3],x[4],x[5]),1,1,1,1);
        double[] accel = forces.acceleration(time,origin,sc);
        out[0] = x[3];
        out[1] = x[4];
        out[2] = x[5];
        out[3] = accel[0];
        out[4] = accel[1];
        out[5] = accel[2];
        for(int i=6; i<x.length; i++) out[i] = 0;
        return out;
	}
	
    public static void main(String[] args) throws IOException{
    	SolarSystemModel.test();
    	//double mjd_utc = 53602;
    	//SolarSystemModel test = new SolarSystemModel(mjd_utc);
    	//test.optimize();
    }

    /**
     * test method
     * Likely there are units problems in this test.
     */
    public static void test() {
    	double start = System.currentTimeMillis();
    	Celestia cel = new Celestia("C:/Code/Celestia/");
    	Trajectory traj = new Trajectory();
    	double mjd_utc = 53602; //53571;//53745.672118;
    	SolarSystemModel test = new SolarSystemModel(mjd_utc);
    	double step_dt = 86400;
//    	double[] in = new double[6];
//    	for(int i=0; i<100; i++){
//    		in = test.body[SolarSystemModel.MARS].getPosition()
//					.minus(test.body[SolarSystemModel.SUN].getPosition())
//					.times(0.001).x;
//    		//traj.add(test.time.mjd_utc(),in);
//    		test.update_dt(step_dt);
//    	}
//    	in = test.body[SolarSystemModel.MARS].getPosition()
//		.minus(test.body[SolarSystemModel.SUN].getPosition())
//		.times(0.001).x;
    	//traj.add(test.time.mjd_utc(),in);
    	
    	VectorN r = test.body[SolarSystemModel.EARTH].getPosition().plus(new VectorN(100000,100000,100000));
    	VectorN v = test.body[SolarSystemModel.EARTH].getVelocity();
    	//Spacecraft sc = new Spacecraft(r,v, 20, 1, 20, 1000);
    	Lambert lam = new Lambert(test.convert_GM_SI(SolarSystemModel.GM_body[SolarSystemModel.SUN]));
    	DE405 ephem2 = new DE405();
    	double travel_time = 7*29; //days
    	double mjd_tt = test.time.mjd_tt()+travel_time;
    	VectorN tempvec = ephem2.get_planet_posvel(DE405_Body.MARS, mjd_tt);
    	double[] xrf = tempvec.get(0, 3).x;
    	double[] xvf = tempvec.get(3, 3).x;
    	EarthRef earth = new EarthRef(new Time(test.time.mjd_utc()+travel_time));
    	
		RotationMatrix rot = earth.EclMatrix(test.time.mjd_tdb());
		
		VectorN rf = new VectorN(xrf[1]*1000,xrf[2]*1000,xrf[3]*1000);
    	VectorN vf = new VectorN(xvf[1]*1000,xvf[2]*1000,xvf[3]*1000);		
    	rf = rot.transform(rf);
		vf = rot.transform(vf);
		vf = vf.unitVector().times(Math.sqrt(test.convert_GM_SI(SolarSystemModel.GM_body[SolarSystemModel.SUN])/rf.mag()));
    	double dt = travel_time*86400;
    	lam.compute(r,v,rf,vf,dt);
    	//v = v.plus(lam.deltav0.times(1));//0.9999999999999));
    	//v = new VectorN(18768.5485896016, 27159.53500768568, 1437.7670601153093);
    	v = new VectorN(18768.56762641506 , 27159.504157982643 , 1437.7664350751745);
    	double t = 0;
    	RungeKutta8 rk8 = new RungeKutta8(1800);
    	double[] x = {r.x[0],r.x[1],r.x[2],v.x[0],v.x[1],v.x[2]};
    	double[] xp = {x[0]/1000,x[1]/1000,x[2]/1000};
    	traj.add(test.time.mjd_utc(),xp);
    	double tsim = 0.5;
    	boolean flag=true;
    	for(int i=0; i<tsim*dt/1800; i++){
    		x = rk8.step(t,x,test);
    		xp[0] = x[0]/1000;
    		xp[1] = x[1]/1000;
    		xp[2] = x[2]/1000;
    		t = t+1800;
//    		if(t>dt && flag){
//    			x[3] = x[3] + lam.deltavf.x[0]*0.98;
//    			x[4] = x[4] + lam.deltavf.x[1]*0.98;
//    			x[5] = x[5] + lam.deltavf.x[2]*0.98;
//    			flag = false;
//    		}
    		test.time.update(t);
    		traj.add(test.time.mjd_utc(),xp);
    		System.out.println("step: "+t/(3*tsim*dt));
    	}
    	r = new VectorN(x[0],x[1],x[2]);
    	v = new VectorN(x[3],x[4],x[5]);
    	double dt2 = (travel_time + mjd_utc - test.time.mjd_utc())*86400.0;
    	lam.compute(r,v,rf,vf,dt2);
    	x[3] = x[3] + lam.deltav0.x[0];
    	x[4] = x[4] + lam.deltav0.x[1];
    	x[5] = x[5] + lam.deltav0.x[2];
    	
    	for(int i=0; i<1.1*tsim*dt/1800; i++){
    		x = rk8.step(t,x,test);
    		xp[0] = x[0]/1000;
    		xp[1] = x[1]/1000;
    		xp[2] = x[2]/1000;
    		t = t+1800;
    		if(t>dt && flag){
    			x[3] = x[3] + lam.deltavf.x[0];
    			x[4] = x[4] + lam.deltavf.x[1];
    			x[5] = x[5] + lam.deltavf.x[2];
    			flag = false;
    		}
    		test.time.update(t);
    		traj.add(test.time.mjd_utc(),xp);
    		System.out.println("step: "+t/(3*tsim*dt));
    	}
    	
    	double mjd3 = test.time.mjd_utc()+100;
    	VectorN tmp = test.ephem.get_planet_posvel(DE405_Body.MARS, mjd3);
    	rot = earth.EclMatrix(test.time.mjd_tdb());
    	rf = tmp.get(0, 3).times(1000);
    	vf = tmp.get(3, 3).times(1000);
    	rf = rot.transform(rf);
    	vf = rot.transform(vf);
    	r = new VectorN(x[0],x[1],x[2]);
    	v = new VectorN(x[3],x[4],x[5]);
    	double dt3 = (100)*86400.0;
    	lam.compute(r,v,rf,vf,dt3);
    	x[3] = x[3] + lam.deltav0.x[0];
    	x[4] = x[4] + lam.deltav0.x[1];
    	x[5] = x[5] + lam.deltav0.x[2];
    	
    	while(test.time.mjd_utc() < mjd3){
    		x = rk8.step(t,x,test);
    		xp[0] = x[0]/1000;
    		xp[1] = x[1]/1000;
    		xp[2] = x[2]/1000;
    		t = t+1800;
    		test.time.update(t);
    		traj.add(test.time.mjd_utc(),xp);
    		System.out.println("step: "+t/((mjd3+100-mjd_utc)*86400));
    	}
    	
    	x[3] = x[3] + lam.deltavf.x[0];
    	x[4] = x[4] + lam.deltavf.x[1];
    	x[5] = x[5] + lam.deltavf.x[2];
    	
    	while(test.time.mjd_utc() < mjd3+100){
    		x = rk8.step(t,x,test);
    		xp[0] = x[0]/1000;
    		xp[1] = x[1]/1000;
    		xp[2] = x[2]/1000;
    		t = t+1800;
    		test.time.update(t);
    		traj.add(test.time.mjd_utc(),xp);
    		System.out.println("step: "+t/((mjd3+100-mjd_utc)*86400));
    	}
    	
    	cel.set_trajectory(traj);
   	    try {
			cel.write_xyz("jat_MRO");
			cel.write_ssc_heliocentric("jat_MRO", "jat_MRO", TimeUtils.MJDtoJD(mjd_utc));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		double elapsed = (System.currentTimeMillis()-start)*0.001/60;
        System.out.println("Elapsed time [min]: "+elapsed);
    	System.out.println("Finished.");
    	String fs = FileUtil.file_separator();
		String dir = FileUtil.getClassFilePath("jat.audio", "SoundPlayer")+"sounds";
		String file = dir+fs+"itisdone.wav";
		//SoundPlayer.play(file);

    }
    
    public void optimize(){
//    	 TODO Auto-generated method stub
    	double start = System.currentTimeMillis();
    	Celestia cel = new Celestia("C:/Code/Celestia/");
    	Trajectory traj = new Trajectory();
    	double mjd_utc = 53602; //53571;//53745.672118;
    	this.epoch = mjd_utc;
    	Time time = new Time(mjd_utc);
		this.update(time);
//		VectorN r = this.body[SolarSystemModel.EARTH].getPosition().plus(new VectorN(100000,100000,100000));
    	VectorN v = this.body[SolarSystemModel.EARTH].getVelocity();
    	//Spacecraft sc = new Spacecraft(r,v, 20, 1, 20, 1000); 
    	VectorN tempvec = ephem2.get_planet_posvel(DE405_Body.MARS, time.mjd_tt());
    	double[] xrf = tempvec.get(0, 3).x;
    	double[] xvf = tempvec.get(3, 3).x;
    	
		VectorN rf = new VectorN(xrf[1]*1000,xrf[2]*1000,xrf[3]*1000);
    	VectorN vf = new VectorN(xvf[1]*1000,xvf[2]*1000,xvf[3]*1000);		
    	rf = rotf.transform(rf);
		vf = rotf.transform(vf);
		vf = vf.unitVector().times(Math.sqrt(this.convert_GM_SI(SolarSystemModel.GM_body[SolarSystemModel.SUN])/rf.mag()));
    	double dt = travel_time*86400;
    	lam.compute(r0,v,rf,vf,dt);
    	v = v.plus(lam.deltav0);
    	
    	//v = new VectorN(18768.5485896016, 27159.53500768568, 1437.7670601153093);
    	v = new VectorN(18768.56690768879, 27159.50563378906, 1437.766446555035);
    	double[] xinit = v.x;
    	//GradientSearch gs = new GradientSearch(this,xinit);
    	DFP opt = new DFP(this,xinit);
    	//DFP.max_it = 1;
    	opt.max_it = 5;
    	double[] xsol = opt.find_min_DFP();
    	
    	v = new VectorN(xsol);  
    	double t = 0;
    	RungeKutta8 rk8 = new RungeKutta8(1800);
    	double[] x = {r0.x[0],r0.x[1],r0.x[2],v.x[0],v.x[1],v.x[2]};
    	double[] xp = {x[0]/1000,x[1]/1000,x[2]/1000};
    	traj.add(epoch,xp);
    	double tsim = 1.5;
    	boolean flag=false;
    	for(int i=0; i<tsim*dt/1800; i++){
    		x = rk8.step(t,x,this);
    		xp[0] = x[0]/1000;
    		xp[1] = x[1]/1000;
    		xp[2] = x[2]/1000;
    		t = t+1800;
    		if(t>0.9*dt && flag){
    			x[3] = x[3] + lam.deltavf.x[0]*0.98;
    			x[4] = x[4] + lam.deltavf.x[1]*0.98;
    			x[5] = x[5] + lam.deltavf.x[2]*0.98;
    			flag = false;
    		}
    		this.time.update(t);
    		traj.add(this.time.mjd_utc(),xp);
    		//System.out.println("step: "+t/(tsim*dt));
    	}
    	
    	cel.set_trajectory(traj);
    	try {
			cel.write_heliocentric("test_sc2","test_sc2",TimeUtils.MJDtoJD(mjd_utc));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		double elapsed = (System.currentTimeMillis()-start)*0.001/60;
        System.out.println("Elapsed time [min]: "+elapsed);
    	System.out.println("Finished.");
    	String fs = FileUtil.file_separator();
		String dir = FileUtil.getClassFilePath("jat.audio", "SoundPlayer")+"sounds";
		String file = dir+fs+"itisdone.wav";
		SoundPlayer.play(file);

    }
    
	/* (non-Javadoc)
	 * @see jat.alg.ScalarfromArrayFunction#evaluate(double[])
	 */
	public double evaluate(double[] x) {
		
		double t = 0;		
    	RungeKutta8 rk8 = new RungeKutta8(1800);
    	double[] xs = {r0.x[0],r0.x[1],r0.x[2],x[0],x[1],x[2]};
    	double[] xp = {x[0]/1000,x[1]/1000,x[2]/1000};
    	double tsim = 1.1;
    	VectorN m,r;
    	double dt = travel_time*86400;
    	double tmp1=this.body[SolarSystemModel.MARS].getPosition().mag();
    	double tmp2=0;
//    	boolean flag=true;
    	for(int i=0; i<tsim*dt/1800; i++){
    		xs = rk8.step(t,xs,this);
    		xp[0] = xs[0]/1000;
    		xp[1] = xs[1]/1000;
    		xp[2] = xs[2]/1000;    		
    		t = t+1800;    		
    		update(t);
    		m = this.body[SolarSystemModel.MARS].getPosition();
    		r = new VectorN(xs[0],xs[1],xs[2]);
    		tmp2 = m.minus(r).mag();
    		if(tmp2<tmp1) tmp1=tmp2;
//    		if(t>dt && flag){
//    			xs[3] = xs[3] + lam.deltavf.x[0]*0.98;
//    			xs[4] = xs[4] + lam.deltavf.x[1]*0.98;
//    			xs[5] = xs[5] + lam.deltavf.x[2]*0.98;
//    			flag = false;
//    		}
    		//this.time.update(t);
    		//traj.add(this.time.mjd_utc(),xp);
    		//System.out.println("step: "+t/(tsim*dt));
    	}
    	this.update(travel_time*86400);
    	//m = this.body[SolarSystemModel.MARS].getPosition();
    	//VectorN rf = new VectorN(xs[0],xs[1],xs[2]);
    	double val = tmp1; //(m.minus(rf)).mag();
    	//System.out.println("xinit: "+x[0]+" "+x[1]+" "+x[2]);
    	//System.out.println("*** val: "+val);
		return val;
	}
}
