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
package jat.attitude;

/**
 * <P>The QuatToDeg converts the quarternions to the corresponding
 *    euler angles for physical interpretation. 
 *
 * @author Noriko Takada
 * @version Last Updated:(8/112/2003)
 * 
 * Reference:	"Spacecraft Dynamics & Control- A Practical Enginnering
 *               Approach" by Marcel J. Sidi
 * 				p. 102-103, p.323
 * 
 * Note: This class uses 3-2-1 Euler angle rotation
 * 			psi (3)
 * 			theta (2)
 * 			phi	(1)
 */
public class DegToQuat 
{
	public double psi = 60;
	public double theta = 10;
	public double phi = 120;
	public double q1;
	public double q2;
	public double q3;
	public double q4;
	
	/**
	 * Constructor 1
	 */
 	public DegToQuat(double psi, double theta, double phi)
 	{
 		this.psi = psi;
 		this.theta = theta;
 		this.phi = phi;
 	}
 	
 	// Default Constructor
 	/**
 	 * Default Constructor
 	 */
 	public DegToQuat()
 	{
 	}
	
	/*
	%Convert degrees to radians
	psi=psi*pi/180;
	theta=theta*pi/180;
	phi=phi*pi/180;

	%P.102 321 rotation Transformation matrix
	a11=cos(theta)*cos(psi);
	a12=cos(theta)*sin(psi);
	a13=-sin(theta);
	a21=-cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
	a22=cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi);
	a23=sin(phi)*cos(theta);
	a31=sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
	a32=-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);
	a33=cos(phi)*cos(theta);

	%P.323- Definitions of Quaternion elements
	trace=a11+a22+a33;
	cos_alpha=(trace-1)/2;
	alpha=acos(cos_alpha);

	e1=(a23-a32)/(2*sin(alpha));
	e2=(a31-a13)/(2*sin(alpha));
	e3=(a12-a21)/(2*sin(alpha));

	q1=e1*sin(alpha/2)
	q2=e2*sin(alpha/2)
	q3=e3*sin(alpha/2)
	q4=cos(alpha/2)
	*/
	
	public double[] calculateQuat()
	{
		double [] out = new double[4];
		
		// Convert degrees to radians
		psi = psi*Math.PI/180;
		theta = theta*Math.PI/180;
		phi = phi*Math.PI/180;
		
		// Sidi p.102, 321 rotation Transformation matrix
		double a11=Math.cos(theta)*Math.cos(psi);
		double a12=Math.cos(theta)*Math.sin(psi);
		double a13=-Math.sin(theta);
		double a21=-Math.cos(phi)*Math.sin(psi)
					+Math.sin(phi)*Math.sin(theta)*Math.cos(psi);
		double a22=Math.cos(phi)*Math.cos(psi)
					+Math.sin(phi)*Math.sin(theta)*Math.sin(psi);
		double a23=Math.sin(phi)*Math.cos(theta);
		double a31=Math.sin(phi)*Math.sin(psi)
					+Math.cos(phi)*Math.sin(theta)*Math.cos(psi);
		double a32=-Math.sin(phi)*Math.cos(psi)
					+Math.cos(phi)*Math.sin(theta)*Math.sin(psi);
		double a33=Math.cos(phi)*Math.cos(theta);
		
		// Sidi p.323 - Definition of Quaternion elements
		double trace = a11 + a22 + a33;
		double cos_alpha  = (trace-1)/2;
		double alpha = Math.acos(cos_alpha);
		
		double e1=(a23-a32)/(2*Math.sin(alpha));
		double e2=(a31-a13)/(2*Math.sin(alpha));
		double e3=(a12-a21)/(2*Math.sin(alpha));
		
		
		q1=e1*Math.sin(alpha/2);
		q2=e2*Math.sin(alpha/2);
		q3=e3*Math.sin(alpha/2);
		q4=Math.cos(alpha/2);
	
		out[0] = q1;
		out[1] = q2;
		out[2] = q3;
		out[3] = q4;
		
		return out;
	}
	
	public static void main(String[] args)
    {
     	double q1 = 0.0;
     	double q2 = 0.0;
     	double q3 = 0.0;
     	double q4 = 1.0;
     	double theta = 10.0;
     	double psi = 60.0;
     	double phi = 120.0;
    
    	q1 = 0.7254;
        q2 = 0.4691;
        q3 = 0.1837;
        q4 = 0.4691;
    	DegToQuat tester = new DegToQuat(psi, theta, phi);
    	
    	double[] quat = new double[4];
    	quat = tester.calculateQuat();
    	q1 = quat[0];
    	q2 = quat[1];
    	q3 = quat[2];
    	q4 = quat[3];

        System.out.println("\nInput: ["+psi+", "+theta+", "+phi+"]");
        System.out.println("\nOutput: q1="+q1+" q2="+q2+" q3="+q3+" q4="+q4);
        
    }
}
