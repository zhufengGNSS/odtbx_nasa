package jat.traj;

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
 * 
 * File Created on Jun 8, 2003
 */

import java.awt.Color;
import java.awt.Graphics;

import jat.matvec.data.*;
import jat.alg.integrators.*;
import jat.plot.*;
import jat.spacetime.*;

/**
* The RelativeTraj.java Class takes two trajectories and computes the relative trajectory.
* The relative trajectory is defined as the difference between them.
* It assumes the input trajectories are in ECI and rotates the differences to RSW.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/
public class RelativeTraj {

	private Trajectory chaser;
	private Trajectory target;
	private double omega;
	private LinePrinter lp;
	SinglePlot traj_plot = new SinglePlot();
	TwoPlots rss_plot = new TwoPlots();
	ThreePlots xyz_plot = new ThreePlots();
	ThreePlots vel_plot = new ThreePlots();
	private boolean verbose = true;
	
	/** Create the relative trajectory from two trajectories
	 * @param sv the satellite trajectory
	 * @param tgt the target or reference trajectory
	 * @param l LinePrinter
	 */
	public RelativeTraj(Trajectory sv, Trajectory tgt, LinePrinter l) {
		this.chaser = sv;
		this.target = tgt;
		this.lp = l;
        traj_plot.setTitle("Relative Trajectory");
        traj_plot.plot.setXLabel("along track");
        traj_plot.plot.setYLabel("radial");
        
        rss_plot.setTitle("RSS Error");
        rss_plot.topPlot.setXLabel("t");
        rss_plot.topPlot.setYLabel("RSS Position");
        rss_plot.bottomPlot.setXLabel("t");
        rss_plot.bottomPlot.setYLabel("RSS Velocity");
        
        xyz_plot.setTitle("Relative Position vs Time");
        xyz_plot.bottomPlot.setXLabel("t");
        xyz_plot.topPlot.setYLabel("radial");
        xyz_plot.middlePlot.setYLabel("along track");
        xyz_plot.bottomPlot.setYLabel("cross track");

        vel_plot.setTitle("Relative Velocity vs Time");
        vel_plot.bottomPlot.setXLabel("t");
        vel_plot.topPlot.setYLabel("radial");
        vel_plot.middlePlot.setYLabel("along track");
        vel_plot.bottomPlot.setYLabel("cross track");
		
	}

	/** Create the relative trajectory from two trajectories
	 * @param sv the satellite trajectory
	 * @param tgt the target or reference trajectory
	 * @param l LinePrinter
	 * @param title the name of the set of figures to be generated e.g. name of spacecraft
	 */
	public RelativeTraj(Trajectory sv, Trajectory tgt, LinePrinter l, String title) {
		this.chaser = sv;
		this.target = tgt;
		this.lp = l;
        traj_plot.setTitle(title+" : Relative Trajectory");
        traj_plot.plot.setXLabel("along track");
        traj_plot.plot.setYLabel("radial");
        
        rss_plot.setTitle(title+" : RSS Error");
        rss_plot.topPlot.setXLabel("t");
        rss_plot.topPlot.setYLabel("RSS Position");
        rss_plot.bottomPlot.setXLabel("t");
        rss_plot.bottomPlot.setYLabel("RSS Velocity");
        
        xyz_plot.setTitle(title+" : Relative Position vs Time");
        xyz_plot.bottomPlot.setXLabel("t");
        xyz_plot.topPlot.setYLabel("radial");
        xyz_plot.middlePlot.setYLabel("along track");
        xyz_plot.bottomPlot.setYLabel("cross track");

        vel_plot.setTitle(title+" : Relative Velocity vs Time");
        vel_plot.bottomPlot.setXLabel("t");
        vel_plot.topPlot.setYLabel("radial");
        vel_plot.middlePlot.setYLabel("along track");
        vel_plot.bottomPlot.setYLabel("cross track");
		
	}


	/** Implements the Printable interface
	 * @param t time
	 * @param y double[] containing data
	 */
	public void print(double t, VectorN x) {

		// handle the first variable for plotting - this is a little mystery but it works
		boolean first = true;
		if (t == 0.0)
			first = false;

		// add data point to the plot
		this.traj_plot.plot.addPoint(0, x.x[1], x.x[0], first);
		this.vel_plot.topPlot.addPoint(0, t, x.x[3], first);
		this.vel_plot.middlePlot.addPoint(0, t, x.x[4], first);
		this.vel_plot.bottomPlot.addPoint(0, t, x.x[5], first);
		this.xyz_plot.topPlot.addPoint(0, t, x.x[0], first);
		this.xyz_plot.middlePlot.addPoint(0, t, x.x[1], first);
		this.xyz_plot.bottomPlot.addPoint(0, t, x.x[2], first);

		// also print to the screen for warm fuzzy feeling or plotting in matlab (boo!)
		if(verbose){
		    this.lp.println(t + "\t" + x.toString());
		    System.out.println(t + "\t" + x.toString());
		}
	}
	/** Implements the Printable interface
	 * @param t time
	 * @param y double[] containing data
	 */
	public void printRSS(double t, double pos, double vel) {

		// handle the first variable for plotting - this is a little mystery but it works
		boolean first = true;
		if (t == 0.0)
			first = false;

		// add data point to the plot
		this.rss_plot.topPlot.addPoint(0, t, pos, first);
		this.rss_plot.bottomPlot.addPoint(0, t, vel, first);

		// also print to the screen for warm fuzzy feeling or plotting in matlab (boo!)
		if(verbose){
		    this.lp.println(t + "\t" + pos);
		    System.out.println(t + "\t" + vel);
		}
	}	
	/** Implements the Printable interface
	 * @param t time
	 * @param y double[] containing data
	 */
	public void printCov(double t, VectorN x, int c) {

		// handle the first variable for plotting - this is a little mystery but it works
		boolean first = true;
		if (t == 0.0)
			first = false;
		Color tmp;
		// add data point to the plot
		//this.traj_plot.plot.addPoint(c, x.x[1], x.x[0], first);
		this.vel_plot.topPlot.addPoint(1, t, x.x[3], first);
		this.vel_plot.middlePlot.addPoint(1, t, x.x[4], first);
		this.vel_plot.bottomPlot.addPoint(1, t, x.x[5], first);		
		this.xyz_plot.topPlot.addPoint(1, t, x.x[0], first);
		this.xyz_plot.middlePlot.addPoint(1, t, x.x[1], first);
		this.xyz_plot.bottomPlot.addPoint(1, t, x.x[2], first);
		//plot the negative
		this.vel_plot.topPlot.addPoint(2, t, -x.x[3], first);
		this.vel_plot.middlePlot.addPoint(2, t, -x.x[4], first);
		this.vel_plot.bottomPlot.addPoint(2, t, -x.x[5], first);		
		this.xyz_plot.topPlot.addPoint(2, t, -x.x[0], first);
		this.xyz_plot.middlePlot.addPoint(2, t, -x.x[1], first);
		this.xyz_plot.bottomPlot.addPoint(2, t, -x.x[2], first);
		
		// also print to the screen for warm fuzzy feeling or plotting in matlab (boo!)
		if(verbose){
		    this.lp.println(t + "\t" + x.toString());
		    System.out.println(t + "\t" + x.toString());
		}
	}

	/** Compute the relative trajectory
	 */
	public void process() {
		this.chaser.reset();
		this.target.reset();
		int nchaser = this.chaser.npts();
		int ntgt = this.target.npts();
		int n = 0;
		if (ntgt > nchaser) {
			n = nchaser;
		} else {
			n = ntgt;
		}
		double units = 1;
		double mjd0 = 0;
		if(this.target.getTimeAt(0)!=0){
			units = 3600.0; //hours
			mjd0 = this.target.getTimeAt(0);
		}
		if(verbose)
		    System.out.println("number of points = "+n+" "+nchaser+" "+ntgt);
		
		for (int i = 0; i < n; i++) {
			double[] chase = chaser.next();
			double[] tgt = target.next();
			double t = tgt[0];
			double dt = chase[0] - t;
			if (dt != 0.0 && verbose)
				System.out.println("times out of sync: " + chase[0] + "\t" + t);

			VectorN r = new VectorN(chase[1], chase[2], chase[3]);
			VectorN v = new VectorN(chase[4], chase[5], chase[6]);
			
			VectorN rtgt = new VectorN(tgt[1], tgt[2], tgt[3]);
			VectorN vtgt = new VectorN(tgt[4], tgt[5], tgt[6]);
			
			RSW_Frame rsw = new RSW_Frame(rtgt, vtgt);

			VectorN dr_eci = r.minus(rtgt);
			VectorN dv_eci = v.minus(vtgt);
			
			VectorN x_rsw = rsw.transform(dr_eci, dv_eci);

			t = (t - mjd0)/units;
			this.print(t, x_rsw);
		}
		System.out.println("done processing");
		traj_plot.setVisible(true);
		vel_plot.setVisible(true);
		xyz_plot.setVisible(true);
		this.lp.close();
	}

	/** Compute the relative trajectory
	 * @param tol Tolerance for time mismatch
	 */
	public void process(double tol) {
		this.chaser.reset();
		this.target.reset();
		int nchaser = this.chaser.npts();
		int ntgt = this.target.npts();
		int n = 0;
		if (ntgt > nchaser) {
			n = nchaser;
		} else {
			n = ntgt;
		}
		double units = 1;
		double mjd0 = 0;
		if(this.target.getTimeAt(0)!=0){
			units = 24.0*3600; //hours
			mjd0 = this.target.getTimeAt(0);
		}
		if(verbose)
		    System.out.println("number of points = "+n+" "+nchaser+" "+ntgt);
		int i=0, j=0;
		while( i<nchaser && j<ntgt ){
			double[] chase = chaser.get(i);//chaser.next();
			double[] tgt = target.get(j);//target.next();
			double t = tgt[0];
			double dt = chase[0] - t;
			if(Math.abs(dt) < tol){
				
				VectorN r = new VectorN(chase[1], chase[2], chase[3]);
				VectorN v = new VectorN(chase[4], chase[5], chase[6]);
				
				VectorN rtgt = new VectorN(tgt[1], tgt[2], tgt[3]);
				VectorN vtgt = new VectorN(tgt[4], tgt[5], tgt[6]);
				
				RSW_Frame rsw = new RSW_Frame(rtgt, vtgt);
				
				VectorN dr_eci = r.minus(rtgt);
				VectorN dv_eci = v.minus(vtgt);
				
				VectorN x_rsw = rsw.transform(dr_eci, dv_eci);
				
				t = (t - mjd0)*units;
				this.print(t, x_rsw);
				i++;
				j++;
			} else if(dt<0){
				i++;
			} else {
				j++;
			}
		}
		System.out.println("done processing");
		traj_plot.setVisible(true);
		vel_plot.setVisible(true);
		xyz_plot.setVisible(true);
		this.lp.close();
	}
	
	/** Compute the relative trajectory
	 * @param tol Tolerance for time mismatch
	 */
	public void process(double tol, Trajectory cov) {
		cov.reset();
		int ncov = cov.npts();
		this.chaser.reset();
		this.target.reset();
		int nchaser = this.chaser.npts();
		int ntgt = this.target.npts();
		int n = 0;
		if (ntgt > nchaser) {
			n = nchaser;
		} else {
			n = ntgt;
		}
		double units = 1;
		double mjd0 = 0;
		if(this.target.getTimeAt(0)!=0){
			units = 24.0*3600; //hours
			mjd0 = this.target.getTimeAt(0);
		}
		if(verbose)
		    System.out.println("number of points = "+n+" "+nchaser+" "+ntgt);
		int i=0, j=0;
		while( i<nchaser && j<ntgt ){
			double[] chase = chaser.get(i);//chaser.next();
			double[] tgt = target.get(j);//target.next();			
			double t = tgt[0];
			double dt = chase[0] - t;
			if(Math.abs(dt) < tol){
				
				VectorN r = new VectorN(chase[1], chase[2], chase[3]);
				VectorN v = new VectorN(chase[4], chase[5], chase[6]);
				
				VectorN rtgt = new VectorN(tgt[1], tgt[2], tgt[3]);
				VectorN vtgt = new VectorN(tgt[4], tgt[5], tgt[6]);
				
				RSW_Frame rsw = new RSW_Frame(rtgt, vtgt);
				
				VectorN dr_eci = r.minus(rtgt);
				VectorN dv_eci = v.minus(vtgt);
				
				VectorN x_rsw = rsw.transform(dr_eci, dv_eci);
				
				t = (t - mjd0)*units;
				this.print(t, x_rsw);
				i++;
				j++;
				if( cov.hasNext()){
					double[] dcov = cov.getNext();
					double tt = dcov[0];
					VectorN rr = new VectorN(dcov[1], dcov[2], dcov[3]);
					VectorN vv = new VectorN(dcov[4], dcov[5], dcov[6]);
						
					//RSW_Frame rsw = new RSW_Frame(r, v);
						
						//VectorN dr_eci = r.minus(rtgt);
						//VectorN dv_eci = v.minus(vtgt);
						
						//VectorN xx_rsw = rsw.transform(rr, vv);
						
						tt = (tt - mjd0)*units;
						this.printCov(tt, new VectorN(rr,vv),1);
				}

			} else if(dt<0){
				i++;
			} else {
				j++;
			}
		}
		
		System.out.println("done processing");
		traj_plot.setVisible(true);
		vel_plot.setVisible(true);
		xyz_plot.setVisible(true);
		this.lp.close();
	}
	
	/** Compute the relative trajectory in ECI (X Y Z)
	 * @param tol Tolerance for time mismatch
	 */
	public void process_ECI(double tol) {
		this.chaser.reset();
		this.target.reset();
		int nchaser = this.chaser.npts();
		int ntgt = this.target.npts();
		int n = 0;
		if (ntgt > nchaser) {
			n = nchaser;
		} else {
			n = ntgt;
		}
		if(verbose)
		    System.out.println("number of points = "+n+" "+nchaser+" "+ntgt);
		int i=0, j=0;
		while( i<nchaser && j<ntgt ){
			double[] chase = chaser.get(i);//chaser.next();
			double[] tgt = target.get(j);//target.next();
			double t = tgt[0];
			double dt = chase[0] - t;
			if(Math.abs(dt) < tol){
				
				VectorN r = new VectorN(chase[1], chase[2], chase[3]);
				VectorN v = new VectorN(chase[4], chase[5], chase[6]);
				
				VectorN rtgt = new VectorN(tgt[1], tgt[2], tgt[3]);
				VectorN vtgt = new VectorN(tgt[4], tgt[5], tgt[6]);
				
				//RSW_Frame rsw = new RSW_Frame(rtgt, vtgt);
				
				VectorN dr_eci = r.minus(rtgt);
				VectorN dv_eci = v.minus(vtgt);
				
				//VectorN x_rsw = rsw.transform(dr_eci, dv_eci);
				
				this.print(t, new VectorN(dr_eci,dv_eci));
				i++;
				j++;
			} else if(dt<0){
				i++;
			} else {
				j++;
			}
		}
		traj_plot.plot.setXLabel("along track");
        traj_plot.plot.setYLabel("radial");
        
        xyz_plot.bottomPlot.setXLabel("t");
        xyz_plot.topPlot.setYLabel("X");
        xyz_plot.middlePlot.setYLabel("Y");
        xyz_plot.bottomPlot.setYLabel("Z");

        vel_plot.bottomPlot.setXLabel("t");
        vel_plot.topPlot.setYLabel("X");
        vel_plot.middlePlot.setYLabel("Y");
        vel_plot.bottomPlot.setYLabel("Z");
		System.out.println("done processing");
		//traj_plot.setVisible(true);
		vel_plot.setVisible(true);
		xyz_plot.setVisible(true);
		this.lp.close();
	}
	
	/** Compute the RSS error relative trajectory in ECI (X Y Z)
	 * @param tol Tolerance for time mismatch
	 */
	public void process_RSS(double tol) {
		this.chaser.reset();
		this.target.reset();
		int nchaser = this.chaser.npts();
		int ntgt = this.target.npts();
		int n = 0;
		if (ntgt > nchaser) {
			n = nchaser;
		} else {
			n = ntgt;
		}
		double units = 1;
		double mjd0 = 0;
		if(this.target.getTimeAt(0)!=0){
			units = 24.0*3600; //hours
			mjd0 = this.target.getTimeAt(0);
		}
		if(verbose)
		    System.out.println("number of points = "+n+" "+nchaser+" "+ntgt);
		int i=0, j=0;
		while( i<nchaser && j<ntgt ){
			double[] chase = chaser.get(i);//chaser.next();
			double[] tgt = target.get(j);//target.next();
			double t = tgt[0];
			double dt = chase[0] - t;
			if(Math.abs(dt) < tol){
				
				VectorN r = new VectorN(chase[1], chase[2], chase[3]);
				VectorN v = new VectorN(chase[4], chase[5], chase[6]);
				
				VectorN rtgt = new VectorN(tgt[1], tgt[2], tgt[3]);
				VectorN vtgt = new VectorN(tgt[4], tgt[5], tgt[6]);
				
				//RSW_Frame rsw = new RSW_Frame(rtgt, vtgt);
				
				VectorN dr_eci = r.minus(rtgt);
				VectorN dv_eci = v.minus(vtgt);
				
				//VectorN x_rsw = rsw.transform(dr_eci, dv_eci);
				t = (t - mjd0)*units;
				this.printRSS(t, dr_eci.mag(),dv_eci.mag());
				i++;
				j++;
			} else if(dt<0){
				i++;
			} else {
				j++;
			}
		}
		
		rss_plot.topPlot.setXLabel("t");
        rss_plot.bottomPlot.setXLabel("t");
        rss_plot.topPlot.setYLabel("Position");
        rss_plot.bottomPlot.setYLabel("Velocity");

        System.out.println("done processing");
		//traj_plot.setVisible(true);
		rss_plot.setVisible(true);
		this.lp.close();
	}
	
	/** Run it
	 */
	public static void main(String[] args) {
		String dir = "C:\\Jat\\jat\\traj\\reference\\";
		String chaserfile = "rvtraj_vbar.jat";
		String tgtfile = "isstraj_vbar.jat";
		String out = "C:\\Jat\\jat\\traj\\reference\\truerel_vbar.txt";
		LinePrinter lp = new LinePrinter(out);
		Trajectory chaser = Trajectory.recover(dir + chaserfile);
		Trajectory tgt = Trajectory.recover(dir + tgtfile);
		RelativeTraj rel = new RelativeTraj(chaser, tgt, lp);
		rel.process();

	}
	
	public double get_max_error(){
	    int nchaser = this.chaser.npts();
		int ntgt = this.target.npts();
		int n = 0;
		if (ntgt > nchaser) {
			n = nchaser;
		} else {
			n = ntgt;
		}
		if(verbose)
		    System.out.println("number of points = "+n+" "+nchaser+" "+ntgt);
		double max = 0;
		for (int i = 0; i < n; i++) {
			double[] chase = chaser.get(i);
			double[] tgt = target.get(i);
			double t = tgt[0];
			double dt = chase[0] - t;
			if (dt != 0.0 && verbose)
				System.out.println("times out of sync: " + chase[0] + "\t" + t);

			VectorN r = new VectorN(chase[1], chase[2], chase[3]);
			VectorN v = new VectorN(chase[4], chase[5], chase[6]);
			
			VectorN rtgt = new VectorN(tgt[1], tgt[2], tgt[3]);
			VectorN vtgt = new VectorN(tgt[4], tgt[5], tgt[6]);
			
			RSW_Frame rsw = new RSW_Frame(rtgt, vtgt);

			VectorN dr_eci = r.minus(rtgt);
			VectorN dv_eci = v.minus(vtgt);
			
			VectorN x_rsw = rsw.transform(dr_eci, dv_eci);
			
			if(dr_eci.mag() > max) max = dr_eci.mag();
			//this.print(t, x_rsw);
		}
		System.out.println("done processing relative trajectory");
		
		return max;
	}
	
	private void get_next_time(int[] i){
		i[0]++; i[1]++;
		double one_time = chaser.getTimeAt(i[0]),two_time = target.getTimeAt(i[1]);
		while(i[0]<chaser.size() && i[1]<target.size() && one_time!=two_time){
			if(one_time < two_time)
				i[0]++;
			else
				i[1]++;
			one_time = chaser.getTimeAt(i[0]);
			two_time = target.getTimeAt(i[1]);
		}
	}
	
	public void setVerbose(boolean b){
	    verbose = b;
	}
}
