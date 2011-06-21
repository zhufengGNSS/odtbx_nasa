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
 * @version 1.3(3/07/2004)
 * Modification Record: 
 * 		1. Removed: import java.lang.Math
 */


/**
 * Reference: "Spacecraft Dynamics & Control" by Marcel J. Sidi (p.323+)
 *            "Spacecraft Vehicle Dynamics and Control" by Bong Wie (p.318+)
 *            "Cal Poly Aero 451 Project #2, December 2001" 
 * 
 * Theory: Define e = eigenvector of rotation
 *                  = [e1, e2, e3]
 *                theta = the angle of rotation
 *         Then   q1 = e1*sin(theta/2)
 *				  q2 = e2*sin(theta/2)
 *                q3 = e3*sin(theta/2)
 *                q4 = cos(theta/2)
 *
 * How to find quarternions [q1, q2, q3, q4] from the initial condition
 *         1. Calcurate the axis of rotation
 *		   2.
 *         3
 */
 
 	 
public class QuatToDeg
{
 	
 	// Field (variable) declaration
 	public double q1;
 	public double q2;
 	public double q3;
 	public double q4;
 	public double [][] element = new double[3][3];
 	public double [][] rotation_matrix = new double[3][3];
 	public double theta=10;
 	public double  psi, phi;
 	

 	
 	// Constructor
 	public QuatToDeg(double e1, double e2, double e3, double e4)
 	{
 		q1 = e1;
 		q2 = e2;
 		q3 = e3;
 		q4 = e4;
 		
 	}
 	
 	// Default Constructor
 	public QuatToDeg()
 	{
 		q1 = 0.0;
 		q2 = 0.0;
 		q3 = 0.0;
 		q4 = 1.0;
 	}
 	
 	// Compute the rotation matrix between the reference frame
 	// and to the frame of the particular orientation defined by
 	// the given quaternions
 	public double[] calculateAngle()
 	{
 		double [] out = new double[3];
        
 		rotation_matrix[0][0] = 1- 2*(q2*q2 + q3*q3); //c11
 		//System.out.println("\nc11="+rotation_matrix[0][0]);
 		rotation_matrix[0][1] = 2*(q1*q2 + q3*q4);    //c12
 		//System.out.println("\nc12="+rotation_matrix[0][1]);
 		rotation_matrix[0][2] = 2*(q1*q3 - q2*q4);    //c13
 		//System.out.println("\nc13="+rotation_matrix[0][2]);
 		rotation_matrix[1][0] = 2*(q1*q2 - q3*q4);    //c21
   		rotation_matrix[1][1] = 1- 2*(q1*q1 + q3*q3); //c22
 		rotation_matrix[1][2] = 2*(q2*q3 + q1*q4);    //c23
 		rotation_matrix[2][0] = 2*(q1*q3+q2*q4);      //c31
 		rotation_matrix[2][1] = 2*(q2*q3-q1*q4);      //c32
 		rotation_matrix[2][2] = 1-2*(q1*q1 + q2*q2);  //c33
 		
 			
 		for(int j = 0; j<3; j++)
 		{
 			for(int i = 0; i<3; i++)
 			{
 			
 				element[j][i] = rotation_matrix[j][i];
 				rotation_matrix[j][i] = rotation_matrix[i][j];
 			}
 		}
 		//System.out.println("element[0][2](afterfor)="+element[0][2]);
 	 	
 		double c_psi;
 		double s_psi;
 		double c_phi;
 		double s_phi;
 		
 		Quad_check checkPsi;
 		Quad_check checkPhi;
 		
 		
 		theta = (-1)*Math.asin(element[0][2])*180/Math.PI;
 		//System.out.println("Math.asin(1)="+Math.asin(1.0));
 		//System.out.println("c13="+element[0][2]);
 		//System.out.println("Math.asin(c13)="+Math.asin(element[0][2]));
 		//System.out.println("Math.PI="+Math.PI);
 		//System.out.println("theta="+theta);
 		
 		c_psi = element[0][0]/Math.cos(theta*Math.PI/180);
 		s_psi = element[0][1]/Math.cos(theta*Math.PI/180);
 		System.out.println("c_psi="+c_psi);
 		System.out.println("s_psi="+s_psi);
 		
 		// Quadrant Check
 		checkPsi = new Quad_check(c_psi, s_psi);
 		checkPsi.determineAngle();
 		psi = checkPsi.angle;
 		//System.out.println("psi="+psi);
 		
 		s_phi = element[1][2]/Math.cos(theta*Math.PI/180);
 		c_phi = element[2][2]/Math.cos(theta*Math.PI/180);
 		
 		checkPhi = new Quad_check(c_phi, s_phi);
 		checkPhi.determineAngle();
 		phi = checkPhi.angle;
 		
 		
 		out[0] = theta;
        out[1] = psi;
        out[2] = phi;
        
        return out;
 	}// End of calculateAngle
 	
} // End of QuatToDeg
 	
 
class Quad_check
{
 		// Given a value of cosine and sine of an angle,
 		// determine the angle
 		
 		double c_angle; 
 		double s_angle; 
 		double angle;
 		
 		public Quad_check(double cos, double sin)
 		{   //Lesson: do not use the same variable names
 			// for the parameters; that did not work.
 			c_angle = cos;
 			s_angle = sin;
 		}
 		
 		
 		public void determineAngle()
 		{
 			if(c_angle >0) // 1st or 4th
 			{
 			 	if(s_angle > 0) // 1st
 			 		angle = Math.acos(c_angle)*180/Math.PI;
 			 	else if(s_angle < 0) // 4th
 			 		angle = Math.asin(s_angle)*180/Math.PI;
 			 	else if(s_angle == 0)
 			 		angle = 0;
 			 	else
 			 		System.out.println("error!");
 			 		
 			}
 		
 			else if (c_angle < 0) //2nd or 3rd
 			{
 				if(s_angle > 0) // 2nd
 					angle = Math.acos(c_angle)*180/Math.PI;
 				else if(s_angle < 0) // 3rd
 					angle =(Math.acos(c_angle)+Math.PI)*180/Math.PI;
 				else if(s_angle == 0)
 					angle = 180;
 			}
 			//System.out.println("angle="+angle);
 		}// End of determineAngle()
}
 					 
 			 
 
 