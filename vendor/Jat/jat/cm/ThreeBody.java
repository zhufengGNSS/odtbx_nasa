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

package jat.cm;

import jat.alg.integrators.*;
import jat.matvec.data.*;

/**
 * <P>
 * The ThreeBody class provides the ability to propagate orbits for a general
 * three body problem
 *
 * @author Tobias Berthold
 * @version 1.0
 */

public class ThreeBody implements Derivatives //, Printable
{

	private double G; // Gravitational constant
	private double m1, m2, m3; // masses

	/**
	 * Method ThreeBody.
	 * @param G gravitational constant
	 * @param m1 mass 1
	 * @param m2 mass 2
	 * @param m3 mass 3
	 */
	public ThreeBody(double G, double m1, double m2, double m3)
	{
		this.G = G;
		this.m1 = m1;
		this.m2 = m2;
		this.m3 = m3;
	}

	/** Compute the derivatives.
	 * @param t    double containing time or the independent variable.
	 * @param x    VectorN containing the required data.
	 * @return      double [] containing the derivatives.
	 */
	// vector r1 is x[0], x[1], x[2]
	// vector v1 is x[3], x[4], x[5]
	// vector r2 is x[6], x[7], x[8]
	// vector v2 is x[9], x[10], x[11]
	// vector r3 is x[12], x[13], x[14]
	// vector v3 is x[15], x[16], x[17]
	public double[] derivs(double t, double[] x)
	{
		//double m1, m2, m3;
		double r1, r2, r3;
		double r12, r23, r13;
		double r12cubed, r23cubed, r13cubed;

		double dxdt[] = new double[18];

		r12 = Math.sqrt((x[6] - x[0]) * (x[6] - x[0]) + (x[7] - x[1]) * (x[7] - x[1]) + (x[8] - x[2]) * (x[8] - x[2]));
		r13 = Math.sqrt((x[12] - x[0]) * (x[12] - x[0]) + (x[13] - x[1]) * (x[13] - x[1]) + (x[14] - x[2]) * (x[14] - x[2]));
		r23 = Math.sqrt((x[12] - x[6]) * (x[12] - x[6]) + (x[13] - x[7]) * (x[13] - x[7]) + (x[14] - x[8]) * (x[14] - x[8]));
		r12cubed = r12 * r12 * r12;
		r13cubed = r13 * r13 * r13;
		r23cubed = r23 * r23 * r23;

		// Derivatives
		dxdt[0] = x[3];
		dxdt[1] = x[4];
		dxdt[2] = x[5];
		dxdt[3] = m2 / r12cubed * (x[6] - x[0]) + m3 / r13cubed * (x[12] - x[0]);
		dxdt[4] = m2 / r12cubed * (x[7] - x[1]) + m3 / r13cubed * (x[13] - x[1]);
		dxdt[5] = m2 / r12cubed * (x[8] - x[2]) + m3 / r13cubed * (x[14] - x[2]);
		dxdt[6] = x[9];
		dxdt[7] = x[10];
		dxdt[8] = x[11];
		dxdt[9] = -m1 / r12cubed * (x[6] - x[0]) + m3 / r23cubed * (x[12] - x[6]);
		dxdt[10] = -m1 / r12cubed * (x[7] - x[1]) + m3 / r23cubed * (x[13] - x[7]);
		dxdt[11] = -m1 / r12cubed * (x[8] - x[2]) + m3 / r23cubed * (x[14] - x[8]);
		dxdt[12] = x[15];
		dxdt[13] = x[16];
		dxdt[14] = x[17];
		dxdt[15] = -m1 / r13cubed * (x[12] - x[0]) - m2 / r23cubed * (x[12] - x[6]);
		dxdt[16] = -m1 / r13cubed * (x[13] - x[1]) - m2 / r23cubed * (x[13] - x[7]);
		dxdt[17] = -m1 / r13cubed * (x[14] - x[2]) - m2 / r23cubed * (x[14] - x[8]);
		return dxdt;
	}

	/** Computes the center of mass
	 * @param x current state
	 * @return center of mass
	 */
	public VectorN center_of_mass(double x[])
	{
		// center of mass
		double M = m1 + m2 + m3;
		double cmx = (m1 * x[0] + m2 * x[6] + m3 * x[12]) / M;
		double cmy = (m1 * x[1] + m2 * x[7] + m3 * x[13]) / M;
		double cmz = (m1 * x[2] + m2 * x[8] + m3 * x[14]) / M;
		VectorN V=new VectorN(cmx,cmy,cmz);
		return V;
	}

	/** Computes the value of the energy of the system
	 * @param x current state
	 * @return The current value of energy.
	 */
	public double Energy(double x[])
	{
		double r12, r23, r13;
		double v1squared, v2squared, v3squared;
		double E = 0., T, U; // total, kinetic, potential energy

		// Kinetic energy
		v1squared = x[3] * x[3] + x[4] * x[4] + x[5] * x[5];
		v2squared = x[9] * x[9] + x[10] * x[10] + x[11] * x[11];
		v3squared = x[15] * x[15] + x[16] * x[16] + x[17] * x[17];
		T = .5 * (m1 * v1squared + m2 * v2squared + m3 * v3squared);

		// Potential energy
		r12 = Math.sqrt((x[6] - x[0]) * (x[6] - x[0]) + (x[7] - x[1]) * (x[7] - x[1]) + (x[8] - x[2]) * (x[8] - x[2]));
		r13 = Math.sqrt((x[12] - x[0]) * (x[12] - x[0]) + (x[13] - x[1]) * (x[13] - x[1]) + (x[14] - x[2]) * (x[14] - x[2]));
		r23 = Math.sqrt((x[12] - x[6]) * (x[12] - x[6]) + (x[13] - x[7]) * (x[13] - x[7]) + (x[14] - x[8]) * (x[14] - x[8]));
		U = m1 * m2 / r12 + m2 * m3 / r23 + m1 * m3 / r13;

		// Total energy
		E = T - U;

		return E;
	}
}
