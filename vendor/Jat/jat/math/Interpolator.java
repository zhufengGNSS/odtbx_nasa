/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2002 United States Government as represented by the
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

package jat.math;

import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;

import java.io.*;

/**
 * Linearly interpolates values in double arrays
 *
 * @author Tobias Berthold <P>
 * Date        :   11-19-2002
 */
public class Interpolator
{
    private double[] data;
	private double[] vari;      // var is actually a keyword, so can't use
	private int length;
	private  int lengthx;
	private int lengthy;
	private int lower, upper;
	private int lowerx, lowery, upperx, uppery;
	private Matrix dataMat;
	private VectorN Xaxis;
	private VectorN Yaxis;

    /**
     * Create an interpolator where data holds the given values of the function at variable 
	 * @param variable independent variable
	 * @param data dependent variable
	 */
	public Interpolator(double[] variable, double[] data)
    {
        this.data=data;
        this.vari=variable;
        length=data.length-1;
    }
	
	 /**
     * Create inerpolation data for a 2-D bi-linear interpolation 
	 * This will hold the incoming data matrix as well as the 
	 * X and Y axis parameters
	 * @param variable independent variable
	 * @param data dependent variable
	 */
	public Interpolator(String filename)
    {
		try {
			this.Xaxis = new VectorN(readXAxis(filename));
			this.Yaxis = new VectorN(readYAxis(filename));
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		this.dataMat = readData(filename, this.Xaxis.length, this.Yaxis.length);
		
        lengthy=this.dataMat.getColumnDimension()-1;
        lengthx=this.dataMat.getRowDimension()-1;
    }

    /**
     * Return the interpolated value of the function at x
	 * @param x independent variable
	 * @return interpolated value at x
	 */
	public double get_value(double x)
    {
        if(x<vari[0]) return 0.;    // too low: really should catch exception.. later
        if(x==vari[0]) return data[0];    // on lower boundary
        if(x==vari[length]) return data[length];    // on upper boundary
        if(x>vari[length]) return 0.;    // too high: really should catch exception..

        // Find bracketing values
        int i;
        for(i=0;i<length;i++)
            if(x>vari[i] && x<vari[i+1])
            {
                //System.out.println("x "+x+" i "+i);
                lower =i;
                upper =i+1;
                break;
            }

        double delta_y=data[upper]-data[lower];
        double delta_x=vari[upper]-vari[lower];
        // If two consecutive independent var. values are equal:
        if (delta_x<1.e-20)
            return data[lower];
        double slope=delta_y/delta_x;
        return data[lower]+slope*(x-vari[lower]);
    }
    /**
     * Return the interpolated value of the function at x,y
	 * @param x independent variable
	 * @param y independent variable
	 * @return interpolated value at x,y
	 */
	public double get_value(double x, double y)
    {
		double out;
		if(x<Xaxis.get(0) || y < Yaxis.get(0))     // too low: really should catch exception.. later
		{
			System.out.println("Interpolator:  Interpolation value too low");
			System.exit(0);
		}
        if(x>Xaxis.get(Xaxis.getLength()-1)) // too high: really should catch exception..
        {
			System.out.println("Interpolator:  Interpolation value too High!");
			System.out.println("X is : " + x + " Y is: " + y);
			System.exit(0);	
        }
        if(y>Yaxis.get(Yaxis.getLength()-1))
        {
        	y = y -2;
        }
        
		//Determine the X bracketing Value
		int i;
		for(i=0;i<lengthx;i++)
		{
			if(x>Xaxis.get(i) && x<Xaxis.get(i+1))
			{
				//System.out.println("x "+x+" i "+i);
				lowerx =i;
				upperx =i+1;
				break;
			}
		}
		
//		Determine the Y bracketing Value
		for(i=0;i<lengthy;i++)
		{
			if(y>Yaxis.get(i) && y<Yaxis.get(i+1))
			{
				//System.out.println("x "+x+" i "+i);
				lowery =i;
				uppery =i+1;
				break;
			}
		}
		
		//Isolate the square on the grid where our point lies
		//Number them successivly counter clockwise from the
		//Bottom left
		double y1 = dataMat.get(lowerx,lowery);
		double y2 = dataMat.get(upperx,lowery);
		double y3 = dataMat.get(upperx,uppery);
		double y4 = dataMat.get(lowerx,uppery);
		
		//Define some utility variables
		double t = (x - lowerx)/(upperx-lowerx);
		double u = (y - lowery)/(uppery-lowery);
		
		//Compute the interpolated value
		out = (1-t)*(1-u)*y1 + t*(1-u)*y2 + t*u*y3 + (1-t)*u*y4;
		
		
		return out;
    }
	
    /**
     * Return the X Axis values from a file that is formatted
     * with the first row containing the X-axis values, first
     * column containing the Y axis values, and the remainder
     * of the file containing the data
	 * @param String file input file name and path
	 * @return VectorN Xaxis  vector with xaxis values
	 */
	public static double[] readXAxis(String Filename) throws IOException
	{
	      double[] tempArray = new double[1000]; 
	      int length = 0;
	      try {
	          // Create the tokenizer to read from a file
	          FileReader rd = new FileReader(Filename);
	          StreamTokenizer st = new StreamTokenizer(rd);
	      
	          // Prepare the tokenizer for Java-style tokenizing rules
	          st.parseNumbers();
	          st.wordChars('_', '_');
	          st.eolIsSignificant(true);
	      
	          // If whitespace is not to be discarded, make this call
	          st.ordinaryChars(0, ' ');
	      
	          // These calls caused comments to be discarded
	          st.slashSlashComments(true);
	          st.slashStarComments(true);
	      
	          // Parse the file
	          int token = st.nextToken();
	          while (token != StreamTokenizer.TT_EOL) {
	              
	        	  token = st.nextToken();
	              switch (token) {
	              case StreamTokenizer.TT_NUMBER:
	                  // A number was found; the value is in nval
	                  double num = st.nval;
	                  length ++;
	                  tempArray[length] = num;
	                  break;
	              default:
	                  // A regular character was found; the value is the token itself
	                  char ch = (char)st.ttype;
	                  break;
	              }
	          }
	          rd.close();
	      } catch (IOException e) {
	      }

		double[]out = new double[length];
		//Fill the output array
		for(int i = 1;i<length; i++)
			out[i] = tempArray[i];
		
		return out;
	}
	
	 /**
     * Return the Y Axis values from a file that is formatted
     * with the first row containing the Y-axis values, first
     * column containing the Y axis values, and the remainder
     * of the file containing the data
	 * @param String file input file name and path
	 * @return VectorN Xaxis  vector with yaxis values
	 * 
	 */
	public static double[] readYAxis(String Filename) throws IOException
	{
	      //Define a temp array to hold values until we determine how
		  //many array entries we will require
		  double[] tempArray = new double[1000]; 
	      int length = 0;
	      try {
	          // Create the tokenizer to read from a file
	          FileReader rd = new FileReader(Filename);
	          StreamTokenizer st = new StreamTokenizer(rd);
	      
	          // Prepare the tokenizer for Java-style tokenizing rules
	          st.parseNumbers();
	          st.wordChars('_', '_');
	          st.eolIsSignificant(true);
	      
	          // If whitespace is not to be discarded, make this call
	          st.ordinaryChars(0, ' ');
	      
	          // These calls caused comments to be discarded
	          st.slashSlashComments(true);
	          st.slashStarComments(true);
	      
	          
	          // Skip the first line because we read it in already
	          int token = st.nextToken();

	         	          
	          //Read through the file and only accept the element that 
	          //immediatly follows an EOF 
	          while (token != StreamTokenizer.TT_EOF) {
	             
                  while (token != StreamTokenizer.TT_EOL) {
                	  token = st.nextToken();
                	  if(token == StreamTokenizer.TT_EOF)
                    	  break;
                  }
                  token = st.nextToken();
                  double num = st.nval;
                  length ++;
                  tempArray[length] = num;
                  if(token == StreamTokenizer.TT_EOF)
                	  break;
                  
	          }
	          rd.close();
	      } catch (IOException e) {
	      }

		double[]out = new double[length-1];
		
		//Fill the output array
		for(int i = 1;i<length; i++)
			out[i-1] = tempArray[i];
		
		return out;
	}
	
	/**
     * Return the data values from a file that is formatted
     * with the first row containing the Y-axis values, first
     * column containing the Y axis values, and the remainder
     * of the file containing the data
	 * @param String file input file name and path
	 * @return Matrix data  matrix with read data
	 * 
	 */
	public static Matrix readData(String Filename,double lenx, double leny)
	{
		Matrix out = new Matrix((int)lenx,(int)leny);
        //Define a temp array to hold values until we determine how
	    //many array entries we will require 
		  int x = 0;
		  int y = -1;
	      int length = 0;
	      try {
	          // Create the tokenizer to read from a file
	          FileReader rd = new FileReader(Filename);
	          StreamTokenizer st = new StreamTokenizer(rd);
	      
	          // Prepare the tokenizer for Java-style tokenizing rules
	          st.parseNumbers();
	          st.wordChars('_', '_');
	          st.eolIsSignificant(true);
	      
	          // If whitespace is not to be discarded, make this call
	          st.ordinaryChars(0, ' ');
	      
	          // These calls caused comments to be discarded
	          st.slashSlashComments(true);
	          st.slashStarComments(true);
	      
	          
	          // Skip the first line because we read it in already
	          int token = st.nextToken();
              while (token != StreamTokenizer.TT_EOL) {
              	  token = st.nextToken();
                }
	          token = st.nextToken();
              
              
	          //Read through the file and toss out the first element
	          //record the rest in the output matrix
	          while (token != StreamTokenizer.TT_EOF) {
	        	  token = st.nextToken();  //this removes the first column
	        	  x = 0;  
	        	  y  ++;
	        	  token = st.nextToken();
	        	  while (token != StreamTokenizer.TT_EOL){
	        		  token = st.nextToken();
	       
	        		  double num = st.nval;
	        		  out.set(x,y,num);
	        		  x++;
	        		  if(token == StreamTokenizer.TT_EOF)
	        			  break;
	        		  token = st.nextToken();
	        	  }
	        	  token = st.nextToken();
	          }
	          rd.close();
	      } catch (IOException e) {
	      }
	      
		return out;
	}
	
	public static void main(String[] args) throws InterruptedException, IOException {
		String filename = "C:\\GOESR\\LM22deg.txt";
		//double[] Xout = readXAxis(filename);
		//double[] Yout = readYAxis(filename);
		//Matrix data = readData(filename, Xout.length, Yout.length);
		//System.out.println("Got it");
		//Now make the interpolator
		Interpolator interp = new Interpolator(filename);
		double outVal = interp.get_value(2.5,2.5);
		System.out.println("The Interpolated Value is " + outVal);
	
	
	}
	

}
