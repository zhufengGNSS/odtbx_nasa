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
 * Emergent Space Technologies
 * File created by Richard C. Page III 
 * Some implementation translated from Matlab written by Sun Hur-Diaz.
 **/
package jat.measurements;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import jat.alg.estimators.MeasurementModel;
import jat.alg.integrators.LinePrinter;
import jat.cm.Constants;
import jat.eph.*;
import jat.math.Interpolator;
import jat.math.MathUtils;
import jat.matvec.data.Quaternion;
import jat.matvec.data.RandomNumber;
import jat.matvec.data.VectorN;
import jat.matvec.data.Matrix;
import jat.spacetime.BodyCenteredInertialRef;
import jat.spacetime.EarthFixedRef;
import jat.spacetime.EarthRef;
import jat.spacetime.LunaFixedRef;
import jat.spacetime.LunaRef;
import jat.spacetime.ReferenceFrame;
import jat.spacetime.ReferenceFrameTranslater;
import jat.spacetime.Time;
import jat.spacetime.TimeUtils;
import jat.sim.CEVLunarSim;
import jat.sim.EstimatorSimModel;
import jat.sim.initializer;
import jat.sim.CEVSim;
import jat.util.FileUtil;

public class OpticalMeasurementModel implements MeasurementModel{
	
	public static final int TYPE_YANGLE_STAR = 1;
	public static final int TYPE_YANGLE_LOS = 2;
	public static final int TYPE_RANGE = 3;
	public static final int TYPE_LANDMARK = 4;
    
    /** A list of bodies and the landmarks on them */
    private static ArrayList<CentralBody> gravbody = setupCentralBodies();
	
	public static final int BODY_EARTH = 0;
	public static final int BODY_MOON = 1;
	public static final int BODY_SUN = 2;
	private VectorN H;
	
	private int type;
	private double freq;
	private double t0,tf;
	private VectorN ustar;
	private int cbody;
	private int vbody;
	private Quaternion q;
	private double R;
	private double Rdeg;
	private int i_ydiskbias;
	private int i_yanglebias;
	private int num_landmarks=0;
	private int current_landmark=0;
	
	
	private double mjd0;
	//* global from y_angle
	private double jd0;
	
	//* global from camerr
	private int biasflag = 1;
	
	private DE405 ephem;
	
	public static LinePrinter fobs,fpred;
	
	//private RandomNumber rnd;
	private Random rnd;
	private boolean runMonteCarlo;
	
	public OpticalMeasurementModel(double mjd_epoch,DE405 jpl){
		mjd0=mjd_epoch;
		jd0 = mjd0+2400000.5;
		ephem = jpl;
		//rnd = new RandomNumber();
		rnd = new Random(System.currentTimeMillis());
	}
	
	public OpticalMeasurementModel(HashMap hm, int measNum) {
		mjd0 = initializer.parseDouble(hm,"init.MJD0")+initializer.parseDouble(hm,"init.T0")/86400;
		jd0 = mjd0+2400000.5;
		ephem = new DE405();
		initialize(hm,measNum);
		String dir_in,fs;
		fs = FileUtil.file_separator();
		try{
			dir_in = FileUtil.getClassFilePath("jat.sim","SimModel")+"output"+fs;
		}catch(Exception e){
			dir_in = "";
		}
		int jat_case = CEVSim.JAT_case;
		if(CEVSim.JAT_case == 0) jat_case = CEVLunarSim.JAT_case;
		try{
			runMonteCarlo = initializer.parseBool(hm, "init.runMonteCarlo");
		}catch(NullPointerException ne){runMonteCarlo = false;}
		if(!runMonteCarlo){
			fobs = new LinePrinter(dir_in+"obs_"+jat_case+".txt");
			fpred = new LinePrinter(dir_in+"pred_"+jat_case+".txt");
		}
		rnd = new Random(System.currentTimeMillis());
	}
	
	private void initialize(HashMap hm, int measNum){
		String pref = "MEAS."+measNum+".";
		type = chooseType(initializer.parseString(hm,pref+"type"));
		String tmp;
		switch(type){
		case TYPE_YANGLE_STAR:
			freq = initializer.parseDouble(hm,pref+"frequency");
			t0 = initializer.parseDouble(hm,pref+"t0");
			tf = initializer.parseDouble(hm,pref+"tf");
			ustar = new VectorN(initializer.parseDouble(hm,pref+"ustar.1"),
								initializer.parseDouble(hm,pref+"ustar.2"),
								initializer.parseDouble(hm,pref+"ustar.3"));
			tmp = initializer.parseString(hm,pref+"cbody");
			if(tmp.equalsIgnoreCase("earth"))
				cbody = BODY_EARTH;
			else if(tmp.equalsIgnoreCase("moon"))
				cbody = BODY_MOON;
			else
				cbody = 0;
			R = initializer.parseDouble(hm,pref+"R");
			break;
		case TYPE_YANGLE_LOS:
			freq = initializer.parseDouble(hm,pref+"frequency");
			t0 = initializer.parseDouble(hm,pref+"t0");
			tf = initializer.parseDouble(hm,pref+"tf");
			q = new Quaternion(initializer.parseDouble(hm,pref+"q.1"),
								initializer.parseDouble(hm,pref+"q.2"),
								initializer.parseDouble(hm,pref+"q.3"),
								initializer.parseDouble(hm,pref+"q.4"));
			tmp = initializer.parseString(hm,pref+"cbody");
			if(tmp.equalsIgnoreCase("earth"))
				cbody = BODY_EARTH;
			else if(tmp.equalsIgnoreCase("moon"))
				cbody = BODY_MOON;
			else
				cbody = 0;
			tmp = initializer.parseString(hm,pref+"vbody");
			if(tmp.equalsIgnoreCase("earth"))
				vbody = BODY_EARTH;
			else if(tmp.equalsIgnoreCase("moon"))
				vbody = BODY_MOON;
			else
				vbody = 0;
			R = initializer.parseDouble(hm,pref+"R");
			try{
			i_yanglebias = initializer.parseInt(hm,"FILTER.anglebias");
			}catch(Exception e){ i_yanglebias = 6;}
			break;
		case TYPE_RANGE:
			freq = initializer.parseDouble(hm,pref+"frequency");
			t0 = initializer.parseDouble(hm,pref+"t0");
			tf = initializer.parseDouble(hm,pref+"tf");
			tmp = initializer.parseString(hm,pref+"cbody");
			if(tmp.equalsIgnoreCase("earth"))
				cbody = BODY_EARTH;
			else if(tmp.equalsIgnoreCase("moon"))
				cbody = BODY_MOON;
			else
				cbody = 0;
			tmp = initializer.parseString(hm,pref+"vbody");
			if(tmp.equalsIgnoreCase("earth"))
				vbody = BODY_EARTH;
			else if(tmp.equalsIgnoreCase("moon"))
				vbody = BODY_MOON;
			else
				vbody = 0;
			R = initializer.parseDouble(hm,pref+"R");
			try{
				i_yanglebias = initializer.parseInt(hm,"FILTER.rangebias");
			}catch(Exception e){ i_yanglebias = 6;}
			break;
		case TYPE_LANDMARK:
			freq = initializer.parseDouble(hm,pref+"frequency");
			t0 = initializer.parseDouble(hm,pref+"t0");
			tf = initializer.parseDouble(hm,pref+"tf");
			tmp = initializer.parseString(hm,pref+"cbody");
			if(tmp.equalsIgnoreCase("earth"))
				cbody = BODY_EARTH;
			else if(tmp.equalsIgnoreCase("moon"))
				cbody = BODY_MOON;
			else
				cbody = 0;
			tmp = initializer.parseString(hm,pref+"vbody");
			if(tmp.equalsIgnoreCase("earth"))
				vbody = BODY_EARTH;
			else if(tmp.equalsIgnoreCase("moon"))
				vbody = BODY_MOON;
			else
				vbody = 0;
			R = initializer.parseDouble(hm,pref+"R");
			Rdeg = initializer.parseDouble(hm, pref+"Rdeg");
			num_landmarks = initializer.parseInt(hm, pref+"datacount");
			for(int L=0; L<num_landmarks; L++){
				double diam = initializer.parseDouble(hm, pref+L+".bounding_diameter");
				double lat = initializer.parseDouble(hm, pref+L+".latitude");
				double longd = initializer.parseDouble(hm, pref+L+".longitude");
				double alt = initializer.parseDouble(hm, pref+L+".altitude");
				boolean k = initializer.parseBool(hm, pref+L+".known");
				gravbody.get(vbody).addLandmark(diam, lat, longd, alt,k);
			}			
			break;
		default:
			break;			
		}
	}
	
	public static int chooseType(String s){
		if(s.equalsIgnoreCase("y_angle_star")){
			return OpticalMeasurementModel.TYPE_YANGLE_STAR;
		}else if(s.equalsIgnoreCase("y_angle_los")){
			return OpticalMeasurementModel.TYPE_YANGLE_LOS;
		}else if(s.equalsIgnoreCase("range")){
			return OpticalMeasurementModel.TYPE_RANGE;
		}else if(s.equalsIgnoreCase("landmark")){
			return OpticalMeasurementModel.TYPE_LANDMARK;
		}else{
			return 0;
		}
	}
	
	/**
	 * Code converted from Matlab
	 % [y,dydx]=y_angle(x,t,p,s,inoise)
	 % Computes the predicted cosine angle measurement and its Jacobian given the 
	 % inertial position state, time, and vector observation specification.  The 
	 % use of this function is intended for navigation computation.  Note that
	 % this function can be used for angle measurement between body and star
	 % observation or for a scalar component of the line-of-sight measurement.
	 %
	 % INPUT:
	 %   x      nx1   [x y z xdot ydot zdot other]' inertial position and velocity
	 %                (m,sec)
	 %   t      1x1   simulation time (sec)
	 %   p      1x1   flag indicating type of vector observation
	 %                  1 = Earth
	 %                  2 = Moon
	 %   s      3x1   inertial reference unit vector from the s/c, e.g. star unit 
	 %                vector or a coordinate transformation axis from ECI frame to
	 %                sensor frame.  
	 %   inoise 1x1   flag to indicate apply noise model
	 %
	 % OUTPUT:
	 %   y      1x1   cosine of angle between the unit vector from the s/c to
	 %                the body as specified by p and the inertial reference unit
	 %                vector
	 %   dydx   1x6   Jacobian
	 % 
	 % GLOBAL:
	 %   jd0    1x1   Julian date of initial time when t=0
	 %
	 % Written by S. Hur-Diaz  6/20/2006
	 *
	 */
	private double[] y_angle(VectorN x, double t, int p, VectorN s, int inoise){
		//int lx = x.length;
		
		//% Compute radial distance from Earth
		VectorN r = x.get(0,3);
		//double rmag=r.mag();
		VectorN v = new VectorN(3);
		//% Determine unit vector observation 1 (Either Earth or Moon)
		switch(p){
		case BODY_EARTH: //% Earth obs
			v=r.times(-1.0);
			break;
		case BODY_MOON: //% Moon obs
			//% Get lunar position relative to the Earth
			//* TODO watch units
			VectorN xm=new VectorN(ephem.get_planet_pos(DE405_Body.GEOCENTRIC_MOON, TimeUtils.UTCtoTT(mjd0+t/86400)));
			//getmoon(jd0+t/86400);  
			v=xm.minus(r);
			break;
		default:
			System.err.println("Invalid flag for vector observation 1.");
		break;
		}
		
		VectorN u=v.unitVector();
		Matrix eye = new Matrix(3);
		Matrix dudx=(eye.minus(u.outerProduct(u))).divide(v.mag());  
		
		//% Compute measurement cos(theta)
		double y=s.dotProduct(u);
		
		if (inoise==1){
			double[] arnd_abias = camerr(v.mag(),p); //% 1-sigma noise and bias on angular measurement
			y=Math.cos(Math.acos(y)+arnd_abias[0]*randn()+arnd_abias[1]);
		}else{
			//y= y + x.x[i_yanglebias];
		}
		
		//% Compute Jacobian
		//dydx= -[s'*dudx zeros(1,3)];
		VectorN tmp = new VectorN(s.times(dudx),new VectorN(3)).times(-1.0);
		double[] dydx= tmp.x;
		
		
		//% % Assume the bias is on the pseudo measurement rather than angle itself
		//% if lx==8
		//%     y=y+x(7);
		//%     dydx=[dydx 1 0];
		//% end
		double[] out = new double[x.length+1];
		out[0] = y;
		for(int i=1; i<7; i++) out[i] = dydx[i-1];
		return out;
	}
	
	/**
	 * Returns a normally distributed random variable.
	 */
	private double randn(){
		//* http://www.mathworks.com/access/helpdesk/help/techdoc/matlab.html
		//RandomNumber rnd = new RandomNumber(System.currentTimeMillis());
		//return rnd.normal();
		return rnd.nextGaussian();
		//return 0;
	}
	
	/*
	 * Code converted from Matlab
	 % [arnd,abias]=camerr(r,ibody)
	 % Camera anguler error as provided by Dan Schwab, Boeing
	 % 
	 % INPUT:
	 %   r     1x1   range from the target body (m)
	 %   ibody 1x1   integer indicating Earth (=1) or Moon (=2)
	 % 
	 % OUTPUT:
	 %   arnd  1x1    random error (rad)
	 %   abias 1x1    bias error (rad)
	 %
	 % GLOBAL:
	 %   biasflag   1x1  flag indicating if bias error should be applied
	 %
	 % Written by S. Hur-Diaz   6/20/2006
	 */
	private double[] camerr(double r, int ibody){
		//% Range from Earth in m
		//double[] erange= {1069177850, 213835570, 71278520, 21383560, 13364720};
		double[] erange= {0,0,13364720,21383560,71278520,213835570,1069177850};
		//VectorN erange = new VectorN(erange_tmp);
		
		//% Range from Moon in m
		//double[] mrange = {291342820, 58268560, 19422850, 5826860, 3641790};
		double[] mrange = {1788000,1838000, 3641790,5826860,19422850,58268560,291342820};
		//VectorN mrange = new VectorN(mrange_tmp);
		
		double[] rv = new double[5];
		if (ibody == BODY_EARTH) //% Earth
			rv=erange;
		else if (ibody == BODY_MOON) // % Moon
			rv=mrange;
		else
			System.err.println("Invalid body flag");
		
		
		//double[] angerr_rnd_deg= {0.0022, 0.011, 0.032, 0.105, 0.169};
		double[] angerr_rnd_deg= {0.47,0.44,0.169,0.105,0.032,0.011,0.0022};
		//double[] angerr_bias_deg= {biasflag*0.0046, biasflag*0.023, biasflag*0.070, biasflag*0.235, biasflag*0.375};
		double[] angerr_bias_deg= {biasflag*1.055,biasflag*0.98, biasflag*0.375,biasflag*0.235,biasflag*0.070,biasflag*0.023,biasflag*0.0046};
		
		//% Apollo numbers corresponding to 3km horizon sensing error
		//%angerr_rnd_deg=  [ 0.0002    0.0008    0.0024    0.0080    0.0129]';
		
		//% Interpolate/extrapolate based on actual range
		Interpolator interp1 = new Interpolator(rv,angerr_rnd_deg);
		double arnd=MathUtils.DEG2RAD*(interp1.get_value(r));//interp1(rv,angerr_rnd_deg,r,'linear','extrap')*pi/180;
		interp1 = new Interpolator(rv,angerr_bias_deg);
		double abias=MathUtils.DEG2RAD*(interp1.get_value(r));//interp1(rv,angerr_bias_deg,r,'linear','extrap')*pi/180;
		//*TODO watch this
		//arnd = 0;
		abias = 0;
		double[] out = {arnd, abias};
		return out;
	}
	
	/**
	 * Converted from Matlab
	 * % [y,dydx]=y_disk2(x,t,p,R,inoise)
	 % Similar to y_disk but uses sin formulation which is more accurate than
	 % the tangent formulation.
	 % Compute the sin of the predicted half-disk angle.
	 %
	 % INPUT:
	 %   x       6x1  Earth-centered inertial position and velocity vector (m)
	 %   t       1x1  time (sec)
	 %   p       1x1  integer parameter indicating Earth (=1) or Moon (=2)
	 %   R       1x1  body's mean physical radius (m)
	 %   inoise  1x1  integer noise flag, 1=add noise
	 %
	 % OUTPUT:
	 %   y       1x1  scalar measurement R/r, where r is the distance from the
	 %                s/c to the body
	 %   dydx    1x6  Jacobian of the y wrt x
	 %
	 % GLOBAL:
	 %   jd0     1x1  Julian Date of initial time corresponding to t=0
	 %
	 % Written by S. Hur-Diaz    6/20/2006
	 % 
	 */
	private double[] y_disk2(VectorN state, double t, int p, double R, int inoise){
		//double jd=this.jd0+t/86400;
		//int lx=state.x.length;
		
		VectorN xr = new VectorN(3);
		VectorN pos = state.get(0,3);
		VectorN xm = new VectorN(3);
//		if (p==BODY_EARTH) //% Earth
//			xr=pos;
//		else if (p==BODY_MOON){ //% Moon
			//* TODO watch units
//			xm=ephem.get_Geocentric_Moon_pos(
//					TimeUtils.MJDtoJD(Time.TTtoTDB(Time.UTC2TT(mjd0+t/86400)))).times(1000);
			CentralBody cbod = gravbody.get(cbody);
			CentralBody vbod = gravbody.get(p);
			ReferenceFrameTranslater xlater = 
		        new ReferenceFrameTranslater(cbod.inertialRef, vbod.inertialRef, new Time(mjd0+t/86400.0));
		      //VectorN li = xlater.translatePoint(lf);
		      xr = xlater.translatePoint(pos);
			//xr=(pos.minus(xm));
//		} else
	//		System.err.println("Parameter must be 1 for Earth or 2 for Moon.");
		
		double r=xr.mag();
		double y=R/r;  //% sin(half-angle)		
		
		if (y>1){
			System.out.println("R/r greater than 1, pause");
//			try{
//				System.in.read();
//			}catch(Exception e){}
		}
		
		//% Add noise
		if (inoise==1){
			double[] arnd_abias =camerr(r,p); //% 1-sigma noise and bias on angular measurement
			if(type==TYPE_LANDMARK){
				arnd_abias[0] = arnd_abias[0]*1;
			}
			y=Math.sin(Math.asin(y)+arnd_abias[0]*randn()+arnd_abias[1]);
		}else{
			//y = y + state.x[i_ydiskbias];
		}
		
		//% Compute Jacobian
		VectorN tmp = new VectorN(xr.times(-R/(r*r*r)),new VectorN(3));
		double[] dydx = tmp.x;
		
		//% % Assume the bias is on the pseudo measurement rather than angle itself
		//% if lx==8
		//%     y=y+x(8);
		//%     dydx=[dydx 0 1];
		//% end
		
		
		double[] out = new double[state.length+1];
		out[0] = y;
		for(int i=1; i<7; i++) out[i] = dydx[i-1];
		return out;
	}
	
	/**
	 * Converted from Matlab.
	 * % It either checks to see if there is enough illumination to obtain centroid 
	 % information of the body based on the user-specified limit on the fraction of
	 % illumination or if an earth-fixed landmark position is visible and
	 % illuminated.  The check selection depends on the size of the input
	 % parameter frac_or_xlf.
	 * @param jd
	 * @param cbody
	 * @param state
	 * @param frac_or_xlf
	 * @return
	 */
	private double[] illum(double jd, CentralBody cbody, String eval, VectorN state, double[] frac_or_xlf){
		VectorN r = state.get(0,3);
		//% x is position relative to the central body
		double xnorm=r.mag();
		VectorN xhat=r.unitVector();

        // Figure out the vector from the central body to the sun.
        // We translate the sun's origin to the central body's
        // reference frame.
        Time t = new Time(TimeUtils.JDtoMJD(jd));
        BodyCenteredInertialRef sunRef = 
          new BodyCenteredInertialRef(DE405_Body.SUN);
        ReferenceFrameTranslater xlater =
          new ReferenceFrameTranslater(sunRef, cbody.inertialRef, t);
        VectorN xs = xlater.translatePoint(new VectorN(3));
		VectorN xshat=xs.unitVector();
		
		double R2r=0;
		R2r=cbody.R/xnorm;
		
		int flag = 0;
		double ratio = 0;
		if (frac_or_xlf.length==1){
			double frac=frac_or_xlf[0];
			//% Determine if accurate measurements can be had
			double theta=Math.acos(R2r); //% half-angle of the horizon disk centered at body
			double gam=Math.acos(xshat.dotProduct(xhat));    // % angle between spacecraft and sun centered at body
			ratio=(theta+Math.PI/2-gam)/(2*theta);
			if (ratio>=frac){
				flag=1;  
			} else{
				flag=0;  //% Not enough illumination for accurate measurement
			}
		}else if (frac_or_xlf.length==3){
			VectorN xlf= new VectorN(frac_or_xlf);
			//xlf=xlf(:);
			// Check if the landmark is visible and illuminated
            // We need to translate xlf to an inertial reference frame
            xlater = new ReferenceFrameTranslater(cbody.bodyFixedRef,
                cbody.inertialRef, t);
            VectorN xli = xlater.translatePoint(xlf);
			VectorN xlihat=xli.unitVector();//xli/norm(xli);
			flag=0;          //% initialize
			if (xlihat.dotProduct(xhat) > R2r){
				flag=1;      //% is visible
				if (xlihat.dotProduct(xshat) > 0){
					flag=2;  //% is visible and illuminated
				}
			}
		} else{
			System.err.println("Invalid frac_or_xlf");
		}
		
		double[] out = {flag,ratio};
		return out;
		
	}
	public boolean isIllum(double mjd, VectorN state){
		VectorN xlf = new VectorN(gravbody.get(vbody).getLandmark(current_landmark).lmf);
		double[] in = illum2(TimeUtils.MJDtoJD(mjd),cbody,vbody,state.get(0,3),xlf.x);
		if(in[0]==2) return true;
		else return false;
	}
	/**
	 * Converted from Matlab.
	% More accurate than illum.m
	% It either checks to see if there is enough illumination to obtain centroid 
	% information of the body based on the user-specified limit on the fraction of
	% illumination or if an earth-fixed landmark position is visible and
	% illuminated.  The check selection depends on the size of the input
	% parameter frac_or_xlf.
	% cbody is body that x is relative to
	% vbody is body that you are looking at
	 * @param jd
	 * @param origin
	 * @param target
	 * @param state
	 * @param frac_or_xlf
	 * @return
	 */
	private double[] illum2(double jd, int origin, int target, 
			VectorN state, double[] frac_or_xlf){		
		
		//% xc is position relative to the central body
		VectorN r = state.get(0,3);
		double xnorm=r.mag();
		VectorN xhat=r.unitVector();
		VectorN r_vbody;
//		Figure out the vector from the central body to the sun.
		// We translate the sun's origin to the central body's
		// reference frame.
		Time t = new Time(TimeUtils.JDtoMJD(jd));
		CentralBody originBody,targetBody;
		if(origin == BODY_MOON){
			originBody = gravbody.get(1);
		}else{
			originBody = gravbody.get(0);
		}
		if(target == BODY_MOON){
			targetBody = gravbody.get(1);
		}else{
			targetBody = gravbody.get(0);
		}
		ReferenceFrameTranslater xlater =
			new ReferenceFrameTranslater(originBody.inertialRef, targetBody.inertialRef, t);
		VectorN xsv = xlater.translatePoint(r);
		VectorN xsvhat=xsv.unitVector(); 

		BodyCenteredInertialRef sunRef = 
			new BodyCenteredInertialRef(DE405_Body.GEOCENTRIC_SUN);
		xlater =  new ReferenceFrameTranslater(sunRef, originBody.inertialRef, t);
		VectorN xs = xlater.translatePoint(new VectorN(3));
		VectorN shat=xs.unitVector();

		double R2r = targetBody.R/xsv.mag();
		int flag = 2;
		double ratio = 1;
		//	% Determine if accurate measurements can be had
		if (frac_or_xlf.length==1){
			    double frac=frac_or_xlf[0];
			    //% half-angle of the horizon disk centered at body
			    double theta= Math.acos(R2r); 
			    VectorN z= shat.crossProduct(xsvhat);
			    double znorm=z.mag();
			    if( znorm<1e-20){
			        ratio=1;
			    }else{
			        z=z.divide(znorm);
			        VectorN p=z.crossProduct(shat);
			        double sig=Math.signum(z.dotProduct(p.crossProduct(xsvhat)));
			        double gam=Math.acos(p.dotProduct(xsvhat));
			        ratio=(theta-sig*gam)/(2*theta);
			        ratio=Math.max(0,ratio);
			        ratio=Math.min(1,ratio);
			    }			    
			    if( ratio>=frac){
			        flag=1;  
			    }else{
			        flag=0;  //% Not enough illumination for accurate measurement
			    }
		}else if( frac_or_xlf.length==3){
			    VectorN xlf= new VectorN(frac_or_xlf);
			    //% Check if the landmark is visible and illuminated
			    ReferenceFrame lfr = targetBody.bodyFixedRef;
			    ReferenceFrame lref = targetBody.inertialRef;
			    xlater =  new ReferenceFrameTranslater(lfr, lref , t);
			    VectorN xli = xlater.translatePoint(xlf);
			    VectorN xlihat= xli.divide(xli.mag());
			    flag=0;          //% initialize
			    if( xlihat.dotProduct(xsvhat) > R2r){
			        flag=1;      //% is visible
			        if( xlihat.dotProduct(shat) > 0){
			            flag=2;  //% is visible and illuminated
					}else{
						flag=1;
					}
				}
			}else{
			    System.err.println("Invalid frac_or_xlf");
			}
		int stop;
		if(flag!=2) 
			stop = 0;
		double[] out = {flag,ratio};
		
		//out[0] = 2;
		return out;
		
	}
	public int get_num_landmarks(){
		return this.num_landmarks;
	}
    /**
     * Computes the measurement between the LOS of a landmark on a body 
     * and a star or the pseudorange to the landmark based on the apparent
     * size of the landmark.  Note that the star can be a component vector 
     * of the sensor if the inertial attitude of the spacecraft is known.
     * Note that it assumes the position vector is inertial and relative 
     * to the body that contains the landmark.  
     * @param x state vector containing inertial s/c position and 
     * velocity relative to the body (6 values).  May additionally
     * contain body-fixed coordinates of a landmark.  If not provided,
     * landmark location determined from structure.
     * [x y z xdot ydot zdot xl yl zl]
     * @param t time 
     * @param p index into the gravbody structure which knows landmark
     * positions
     * @param s inertial reference unit vector.  If null is passed in
     * (or a vector of length less than 3) this indicates that the pseudorange
     * should be computed 
     * @param lindex landmark index
     * @param inoise noise flag
     * @return a array of 10 doubles.  The first double is the angle between 
     * unit vector from s/c to landmark and an inertial reference unit vector.
     * The next 9 doubles are the Jacobian matrix (dydx)
     */
    private double[] y_landmark(VectorN x, Time t, int p,
        VectorN s, int lindex, boolean inoise) {
      
      // We'll need a Gaussian random number generator
      //Random randn = new Random(System.currentTimeMillis());    	
      
      double y = 0;
      VectorN dydx = new VectorN(x.length);
      
      // TODO: Do we need t to be relative?
      // jd=jd0+t/86400;
      
      // If the number of states is greater than 6, assume that the 
      // components 7-9 are the body-fixed coordinates of the landmark, 
      // otherwise, it is assumed the landmark is known and specified by 
      // the structure.
      CentralBody body = gravbody.get(p);
      VectorN lf = null;
      if (x.length == 9) {
        lf = x.get(6, 3);
    	  //Landmark mark = new Landmark(20, x.get(6), x.get(7), x.get(8),LunaRef.R_Luna);
    	  //lf = mark.lmf;
      }
      else {
        lf = body.getLandmark(lindex).lmf;
      }
                  
      // Compute the coordinate transformation from body-fixed to inertial
      ReferenceFrameTranslater xlater = 
        new ReferenceFrameTranslater(body.bodyFixedRef, body.inertialRef, t);
      VectorN li = xlater.translatePoint(lf);

      // Compute the vector from the spacecraft to the landmark
      VectorN d = li.minus(x.get(0, 3));
      double dmag = d.mag();
      VectorN dhat = d.unitVector();
      if ((s != null) && (s.length == 3)) {
        // Referece unit vector is provided.  Compute measurment
        // between LOS.
        
        // Compute predicted measurement cos(theta)
        y = dhat.dotProduct(s);

        if (inoise) {
          // Assume similar noise as star and zero bias for now
          // TODO: Make arnd and abias configurable instead of
          // hardcoded.
          double arnd=.0034*Math.PI/180.0;
          //arnd=.07*Math.PI/180.0;
          arnd=Rdeg*Math.PI/180.0;
          double abias=0;
          y=Math.cos(Math.acos(y)+arnd*rnd.nextGaussian()+abias);
        }

        // (1/dmag)*s'*(-eye(3)+dhat*dhat')
        VectorN rdot = s.divide(dmag);
        Matrix tmp = dhat.outerProduct(dhat).minus(Matrix.identity(3,3));
        rdot = rdot.times(tmp);
        dydx.set(0, rdot);
        // vdot is 0
        
        if (x.length == 9) {
          //VectorN ldot = xlater.translatePointBack(rdot);//.times(-1);
          VectorN ldot = xlater.translatePointBack(rdot.times(-1));
          dydx.set(6, ldot);
        }
      }
      else {
        // No reference unit vector is provided.  Compute pseudo range.
        
        // measurement is tan((D/2)/r)
        double D = body.getLandmark(lindex).D;
        y=D*.5/dmag;
          
        if (inoise) {
          // 1-sigma noise and bias on angular measurement
          double[] noiseFactors = camerr(dmag,p);
          double arnd = noiseFactors[0];
          double abias = noiseFactors[1];
          y=Math.tan(Math.atan(y)+arnd*rnd.nextGaussian()+abias);
        }

        double scale = (0.5 * D) / (dmag * dmag * dmag);
        VectorN rpart = d.times(scale);
 
        if (x.length == 9) {
          VectorN ldot = xlater.translatePointBack(rpart).times(-1);
          dydx.set(6, ldot);
        }
      }

      // TODO: Do we really want to shmoosh these into one array
      // just to return it?
      // We combine the angle and the Jacobian into a single array to return it
      double[] out = new double[x.length+1];
      out[0] = y;
      //System.arraycopy(dydx, 0, toReturn, 1, dydx.length);
		for(int i=1; i<x.length+1; i++) out[i] = dydx.x[i-1];
      return out; 
    }
    
	private Matrix nadir_dcm(double jd, VectorN xsc, int cbody, int vbody){

		//xsc=xsc(:);
		//% Need to get x relative to vbody
		//% First get cbody relative to earth
		//xce=feval(cbody.fn,jd);
		VectorN xce = new VectorN(3);
		VectorN xve = new VectorN(3);
		if(cbody==BODY_EARTH){
			xce = new VectorN(3);
		}else if(cbody == BODY_MOON){
			//* TODO watch units
			// (meters)
			xce = new VectorN(ephem.get_planet_pos(DE405_Body.GEOCENTRIC_MOON, TimeUtils.TTtoTDB(TimeUtils.UTCtoTT(TimeUtils.JDtoMJD(jd)))));
		}
		//% Second get vbody relative to earth
		//xve=feval(vbody.fn,jd);
		if(vbody==BODY_EARTH){
			xve = new VectorN(3);
		}else if(vbody == BODY_MOON){
			// (meters)
			xve = new VectorN(ephem.get_planet_pos(DE405_Body.GEOCENTRIC_MOON, TimeUtils.TTtoTDB(TimeUtils.UTCtoTT(TimeUtils.JDtoMJD(jd)))));
		}
		
		//% Finally spacecraft relative to vbody
		VectorN xsv= xsc.get(0,3).plus(xce.minus(xve));

		//% Get unit vector from vbody to spacecraft
		//xsvnorm=norm(xsv(1:3));
		VectorN xsvhat=xsv.unitVector();


		VectorN zm=xsvhat.times(-1.0);
		VectorN z = new VectorN(0.0,0.0,1.0);
		VectorN xm=z.crossProduct(zm);
		xm= xm.unitVector();
		VectorN ym=zm.crossProduct(xm);
		Matrix A_sensor_2_inertial= new Matrix(3);
		A_sensor_2_inertial.A[0][0] = xm.x[0];
		A_sensor_2_inertial.A[1][0] = xm.x[1];
		A_sensor_2_inertial.A[2][0] = xm.x[2];
		A_sensor_2_inertial.A[0][1] = ym.x[0];
		A_sensor_2_inertial.A[1][1] = ym.x[1];
		A_sensor_2_inertial.A[2][1] = ym.x[2];
		A_sensor_2_inertial.A[0][2] = zm.x[0];
		A_sensor_2_inertial.A[1][2] = zm.x[1];
		A_sensor_2_inertial.A[2][2] = zm.x[2];
		
		return A_sensor_2_inertial.transpose();
	}

	public double observedMeasurement(int whichMeas, double t, VectorN x){
		double[] out = new double[1];
		//double t = (mjd-mjd0)*86400;
		double mjd = mjd0+t/86400;
		switch(type){
		case TYPE_YANGLE_STAR:
			out = y_angle(x,t,vbody,ustar,1);
			break;
		case TYPE_YANGLE_LOS:
			//Matrix A = q.quat2DCM();
			Matrix A = nadir_dcm(mjd+2400000.5,x,cbody,vbody);
			VectorN s;
			if(whichMeas==0)
				s = new VectorN(A.getRowVector(0));
			else
				s = new VectorN(A.getRowVector(1));
			out = y_angle(x,t,vbody,s,1);
			break;
		case TYPE_RANGE:
			if(vbody == BODY_EARTH)
				out = y_disk2(x,t,vbody,EarthRef.R_Earth,1);
			else //if(vbody == BODY_MOON)
				out = y_disk2(x,t,vbody,LunaRef.R_Luna,1);
			break;
		case TYPE_LANDMARK:
			Matrix A2 = nadir_dcm(mjd+2400000.5,x,cbody,vbody);
			VectorN s2;
			if(whichMeas==0)
				s2 = new VectorN(A2.getRowVector(0));
			else
				s2 = new VectorN(A2.getRowVector(1));
			out = y_landmark(x,new Time(mjd0+t/86400.0),vbody,s2,current_landmark,true);
			break;
		default:
			return 0;
		}
		
		return out[0];
	}
	
	public double predictedMeasurement(int whichMeas, double t, VectorN x){
		double[] out = new double[1];
		//double t = (mjd-mjd0)*86400;
		double mjd = mjd0+t/86400;
		switch(type){
		case TYPE_YANGLE_STAR:
			out = y_angle(x,t,vbody,ustar,0);
			break;
		case TYPE_YANGLE_LOS:
			VectorN truestate = new VectorN(EstimatorSimModel.truth[0].get_spacecraft().toStateVector());
			Matrix A = nadir_dcm(mjd+2400000.5,truestate,cbody,vbody);//q.quat2DCM();
			VectorN s;
			if(whichMeas==0)
				s = new VectorN(A.getRowVector(0));
			else
				s = new VectorN(A.getRowVector(1));
			out = y_angle(x,t,vbody,s,0);
			break;
		case TYPE_RANGE:
			if(vbody == BODY_EARTH)
				out = y_disk2(x,t,vbody,EarthRef.R_Earth,0);
			else //if(vbody == BODY_MOON)
				out = y_disk2(x,t,vbody,LunaRef.R_Luna,0);
			break;
		case TYPE_LANDMARK:
			Matrix A2 = nadir_dcm(mjd+2400000.5,x,cbody,vbody);
			VectorN s2;
			if(whichMeas==0)
				s2 = new VectorN(A2.getRowVector(0));
			else
				s2 = new VectorN(A2.getRowVector(1));
			out = y_landmark(x,new Time(mjd0+t/86400.0),vbody,s2,current_landmark,false);
			break;
		default:
			return 0;
		}
		H = new VectorN(x.length);
		for(int i=1; i<x.length+1;i++) H.x[i-1] = out[i];
		return out[0];
	}
	
	public int get_type(){
		return type;
	}
	
//	******* Interface Implementation ******* //
	
	
	public VectorN H(VectorN xref) {
		//* Cheating by calculating this in the call to observed measurement
		return this.H;
	}
	
	public double R() {
		return R;
	}
	
	public double zPred(int whichMeas, double t_sim, VectorN state) {
		double out,obs;
		VectorN truth;
		try{
			truth = new VectorN(EstimatorSimModel.truth[0].get_spacecraft().toStateVector());
		}catch(NullPointerException ne){
			truth = new VectorN(CEVLunarSim.truth[0].get_spacecraft().toStateVector());
		}
		double[] in = {1,0};
		if(type==TYPE_LANDMARK){
			VectorN xlf = new VectorN(gravbody.get(vbody).getLandmark(current_landmark).lmf);
			in = illum2(TimeUtils.MJDtoJD(mjd0+t_sim/86400.0),cbody,vbody,truth,xlf.x);
		}
		if(in[0] ==2 || type!=TYPE_LANDMARK){
			obs = observedMeasurement(whichMeas,t_sim,truth);
			double pred = predictedMeasurement(whichMeas, t_sim, state); 
			//* TODO watch this - following is a temporary hack to extract the bias only
			//double pred = predictedMeasurement(whichMeas, t_sim, truth);

			if(obs == 0)
				out = 0.0;
			else
				out = obs-pred;

			if(!runMonteCarlo){
			String typestring = "blank";
			if(type==TYPE_YANGLE_STAR){
				typestring = "yangle_star";
				OpticalMeasurementModel.fobs.println("time: "+t_sim+"  obs: "+Math.acos(obs)*MathUtils.RAD2DEG+"   "+typestring);
				OpticalMeasurementModel.fpred.println("time: "+t_sim+"  obs: "+Math.acos(pred)*MathUtils.RAD2DEG+"   "+typestring);
			}else if(type==TYPE_YANGLE_LOS){
				typestring = "yangle_los";
				OpticalMeasurementModel.fobs.println("time: "+t_sim+"  obs: "+Math.acos(obs)*MathUtils.RAD2DEG+"   "+typestring+" "+whichMeas);
				OpticalMeasurementModel.fpred.println("time: "+t_sim+"  obs: "+Math.acos(pred)*MathUtils.RAD2DEG+"   "+typestring+" "+whichMeas);
			}else if(type== TYPE_RANGE){
				typestring = "range";
				OpticalMeasurementModel.fobs.println("time: "+t_sim+"  obs: "+Math.asin(obs)*MathUtils.RAD2DEG+"   "+typestring);
				OpticalMeasurementModel.fpred.println("time: "+t_sim+"  obs: "+Math.asin(pred)*MathUtils.RAD2DEG+"   "+typestring);
			}else if(type== TYPE_LANDMARK){
//				typestring = "landmark";
//				OpticalMeasurementModel.fobs.println("time: "+t_sim+"  obs: "+obs+"   "+typestring+" L"+current_landmark+"-"+whichMeas);
//				OpticalMeasurementModel.fpred.println("time: "+t_sim+"  obs: "+pred+"   "+typestring+" L"+current_landmark+"-"+whichMeas);
				typestring = "landmark";
				OpticalMeasurementModel.fobs.println(""+t_sim+"  "+obs+"    L: "+current_landmark+" S:"+whichMeas);
				OpticalMeasurementModel.fpred.println(""+t_sim+"  "+pred+"    L: "+current_landmark+" S:"+whichMeas);
			}
			}
		}else{
			out = 0.0;
			this.H = new VectorN(state.length);
			String typestring = "landmark";
			if(!runMonteCarlo){
				OpticalMeasurementModel.fobs.println(""+t_sim+"  "+0+"    L: "+current_landmark+" S:"+whichMeas);
				OpticalMeasurementModel.fpred.println(""+t_sim+"  "+0+"    L: "+current_landmark+" S:"+whichMeas);
			}
		}
		return out;
	}
	
    public static class Landmark {
      
      /** The spatial coordinates for the landmark in the body-fixed
       * reference frame. */
      private final VectorN lmf;
      
      /** The diameter of the landmark. */
      private final double D;
      
      /** Flag indicating whether it is a known or unknown landmark */
      private boolean isKnown = true;
      
      public Landmark(double inD, double lat, double longd, 
          double alt, double planetR) {
    	inD = inD*1000.0;
    	alt = alt*1000.0;
        lmf = latLongAlt2FixedCoords(lat, longd, alt, planetR);
        D = inD;
      }
      public Landmark(double inD, double lat, double longd, 
              double alt, double planetR, boolean known) {
        	inD = inD*1000.0;
        	alt = alt*1000.0;
            lmf = latLongAlt2FixedCoords(lat, longd, alt, planetR);
            D = inD;
            isKnown = known;
      }

      private VectorN latLongAlt2FixedCoords(double lat, double longd, 
          double alt, double planetR) {
        VectorN coords = new VectorN(3);
        double degrees2Radians = Math.PI/180;
        lat = lat * degrees2Radians;
        longd = longd * degrees2Radians;
        double height = planetR + alt;
        coords.set(0, Math.cos(lat) * Math.cos(longd) * height);
        coords.set(1, Math.cos(lat) * Math.sin(longd) * height);
        coords.set(2, Math.sin(lat) * height);
        return coords;
      }
    }
	
	private static class CentralBody {
		public final String name;
		public final double R;
		public final double MU;
        
        /** The body-centered inertial reference frame */
        public final ReferenceFrame inertialRef;
        
        private ArrayList<Landmark> landmarks;
        
        /** The body-centered body-fixed reference frame */
        public final ReferenceFrame bodyFixedRef;
				
		public CentralBody(String inName, double inMu, double inR, 
            ReferenceFrame inInertial, ReferenceFrame inFixed) {
          name = inName;
          MU = inMu;
          R = inR;
          inertialRef = inInertial;
          bodyFixedRef = inFixed;
          landmarks = new ArrayList<Landmark>();
		}
        
        public void addLandmark(double diam, double lat, double longd, 
            double alt) {
          Landmark l = new Landmark(diam, lat, longd, alt, R);
          landmarks.add(l);
        }
        public void addLandmark(double diam, double lat, double longd, 
                double alt,boolean known) {
              Landmark l = new Landmark(diam, lat, longd, alt, R,known);             
              landmarks.add(l);
            }
        public Landmark getLandmark(int i) {
          return landmarks.get(i);
        }
   	}
    
    private static ArrayList<CentralBody> setupCentralBodies() {
      ArrayList<CentralBody> bodies = new ArrayList<CentralBody>();
      
      // Setup Earth
      // TODO: This is a guess at R.  Figure out real R.
//      CentralBody earth = new CentralBody("earth", 398600.436, 6377.8037, 
//          new BodyCenteredInertialRef(DE405.EARTH),
//          new EarthFixedRef());
      CentralBody earth = new CentralBody("earth", EarthRef.GM_Earth, EarthRef.R_Earth, 
              new BodyCenteredInertialRef(DE405_Body.EARTH),
              new EarthFixedRef());
      // Add earth landmark
      // TODO: Read these landmarks from a file somewhere?
      //earth.addLandmark(20, -50.75, 127, 0);
      bodies.add(earth);
      
      // Setup the Moon
//      CentralBody moon = new CentralBody("moon", 4902.801056, 1738,
//          new BodyCenteredInertialRef(DE405.MOON),
//          new LunaFixedRef());
      CentralBody moon = new CentralBody("moon", LunaRef.GM_Luna, LunaRef.R_Luna,
              new BodyCenteredInertialRef(DE405_Body.GEOCENTRIC_MOON),
              new LunaFixedRef());
      // Add moon landmarks
      //moon.addLandmark(20, -50.7, 127.4, 0);
      //moon.addLandmark(20, 0, 127.4, 0);
      //moon.addLandmark(20, 50.7, 127.4, 0);
      bodies.add(moon);
      
      // Setup the Sun
//      CentralBody sun = new CentralBody("sun", 1.32712440018e11, 0,
//          new BodyCenteredInertialRef(DE405.SUN), null);
      CentralBody sun = new CentralBody("sun", Constants.GM_Sun, Constants.R_Sun,
              new BodyCenteredInertialRef(DE405_Body.SUN), null);
      bodies.add(sun);
      return bodies;
    }

	public void currentLandmark(int nmark) {
		this.current_landmark = nmark;		
	}
	
//	public static void main(String args[]){
//		
//		double jd = ;
//		int origin;
//		int target; 
//		VectorN state;
//		double[] frac_or_xlf;
//		OpticalMeasurementModel test = new OpticalMeasurementModel(jd-2400000.5,new DE405());
//		double[] out = test.illum2(jd,origin,target,state,frac_or_xlf);
//	}
}
