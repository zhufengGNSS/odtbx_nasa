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
 */

package jat.util;

import jat.matvec.data.Matrix;
import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * Prints a matrix to the screen, nicely formatted so that columns line up.
 * @author Tobias Berthold
 */
public class PrintMatrix
{
	Matrix M;

	enum mode
	{
		FIXED, VARIABLE
	};

	enum decexp
	{
		DECIMAL, EXPONENTIAL
	};

	DecimalFormat df;
	NumberFormat itf;
	private int MinFracDig = 10;
	private int MinIntDig = 8;
	private int[] fractionDigits;
	private int[] integerDigits;
	public String[] titles ;

	public PrintMatrix(Matrix M)
	{
		this.M = M;
		//default values
		int cols = M.getColumnDimension();
		integerDigits=new int[cols];
		fractionDigits=new int[cols];
		titles=new String[cols];		
		for(int i=0;i<cols;i++)
		{
			integerDigits[i]=MinIntDig;
			fractionDigits[i]=MinFracDig;
			titles[i]="Title";
		}
	}

	public void setMinIntDig(int minIntDig)
	{
		MinIntDig = minIntDig;
	}

	public void setMinFracDig(int minFracDig)
	{
		MinFracDig = minFracDig;
		int cols = M.getColumnDimension();
		for(int i=0;i<cols;i++)
		{
			fractionDigits[i]=MinFracDig;
		}
	}

	private void set_print_formats()
	{
		df = (DecimalFormat) NumberFormat.getInstance();
		df.applyPattern("  ###.########; -###.#######");
		// df.applyPattern("####0.000000; -###0.000000");
		// df.setMinimumFractionDigits(MinFracDig);
		df.setMinimumIntegerDigits(5);
		itf = NumberFormat.getInstance();
		itf.setMinimumIntegerDigits(4);
	}

	public void print()
	{ // i:rows j:columns
		// Column Titles
		for (int j = 0; j < M.n; j++)
		{
			int integerwidth= integerDigits[j];
			int fractionwidth= fractionDigits[j];
			int width=integerwidth+fractionwidth+3;
			String s = titles[j] + "                              ";
			String t = s.substring(0,width);
			System.out.print(t);
		}
		System.out.println("");
		// Data
		set_print_formats();
		// System.out.println(M.m + " X " + M.n + " Matrix:");
		for (int i = 0; i < M.m; i++)
		{
			for (int j = 0; j < M.n; j++)
			{
				df.setMinimumIntegerDigits(integerDigits[j]);
				df.setMinimumFractionDigits(fractionDigits[j]);
				// System.out.print("\t");
				System.out.print(df.format(M.A[i][j]));
			}
			System.out.println("");
		}
	}

	public static void main(String[] args)
	{
		// TODO Auto-generated method stub
	}


}
