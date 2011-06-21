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

/**
 * Title:        Non-linear Equation Solver
 * Description:  Solves non-linear equations
 * Copyright:    Copyright (c) 2002
 * Company:
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

package jat.alg;

import jat.matvec.data.*;
import jat.matvec.data.matrixDecompositions.*;

abstract public class NLESolver
{
    public boolean debug=true;
	public int n;             // number of unknowns
	public VectorN x;        // vector of unknowns
	public VectorN scales;   // vector of typical scaling values
	public Matrix d;          // scaling matrix

	// algorithm parameters
	public final double machineps = 1.1102230246251565E-16;  // machine precision
	public double eps = 1.0E-5; // central diff eps
    //	public final double eps = Math.pow(machineps,(1.0/3.0)); // central diff eps
	public final double term = 0.5E-08;                 // termination criterion
	public final double alpha = 1.0E-08;          // line search algorithm param
    //	public final double kmin = 1.1102231E-16;   // minimum correction factor for newton step
	public final double kmin = 1.0E-30;   // minimum correction factor for newton step
	public final int maxit = 30;           // max number of tries in line search

	public NLESolver (VectorN xin)
	// default constructor
	{
		this.x = new VectorN(xin);
		this.n = xin.length;
		this.d = new Matrix(this.n);      // n x n identity matrix
		this.scales = new VectorN(this.n,1.0); // nx1 vector containing all 1's
	}

	public NLESolver (VectorN xin, VectorN xtyp) // used when scaling
	// constructor, inputs are initial x values and typical x's for scaling
	{
		this.x = new VectorN(xin);
		this.n = xin.length;
		this.d = new Matrix(this.n); // n x n identity matrix
		this.scales = new VectorN(this.n, 1.0);

		// fill in the scaling vector and matrix
		for (int i = 0; i < n; i++)
		{
			this.scales.x[i] = 1.0/xtyp.x[i];
			this.d.A[i][i] = this.scales.x[i];
		}
	}

	abstract public VectorN evaluate(VectorN xin);
	// evaluate function to be supplied by the user

    public void set_initial_guess(VectorN x_guess)
    {
        this.x=x_guess;
    }

	protected double figureOfMerit(VectorN xin)
	// computes the rss of the deviations
	{
		VectorN df = this.d.times(xin);
		double out = 0.5 * df.dotProduct(df);
		return out;
	}

	protected Matrix computeJacobian(VectorN xin)
	// computes the Jacobian using central differences
	{
		VectorN nu = new VectorN(xin);
		VectorN nxf = new VectorN(nu);
		VectorN nxb = new VectorN(nu);
		double [][] cx = new double[n][n];

		for (int i = 0; i < n; i++)
		{
			double dnx = Math.abs(eps*nu.x[i]);
			if (Math.abs(nu.x[i]) <= eps)
			{
				dnx = eps*eps;
			}

			nxf.x[i] = nu.x[i] + dnx;
			dnx = nxf.x[i] - nu.x[i];

			VectorN cf = evaluate(nxf);
			nxf.x[i] = nu.x[i];

			nxb.x[i] = nu.x[i] - dnx;
			VectorN cb = evaluate(nxb);
			nxb.x[i] = nu.x[i];

			for (int j = 0; j < n; j++)
			{
				cx[j][i] = 0.5 * (cf.x[j] - cb.x[j])/dnx;
			}
		}
		Matrix Jac = new Matrix(cx, n, n);
		return Jac;
	}

	protected Matrix computeJacobian2(VectorN xin)
	// computes the Jacobian using 4th order central differences
	{
		VectorN nu = new VectorN(xin);
		VectorN nxf = new VectorN(nu);
		VectorN nxb = new VectorN(nu);
		VectorN nxfh = new VectorN(nu);
		VectorN nxbh = new VectorN(nu);

		double [][] cx = new double[n][n];

		for (int i = 0; i < n; i++)
		{
			double dnx = Math.abs(eps*nu.x[i]);
			if (Math.abs(nu.x[i]) <= eps)
			{
				dnx = eps*eps;
			}

			nxf.x[i] = nu.x[i] + dnx;
			dnx = nxf.x[i] - nu.x[i];

			VectorN cf = evaluate(nxf);
			nxf.x[i] = nu.x[i];

			nxb.x[i] = nu.x[i] - dnx;
			VectorN cb = evaluate(nxb);
			nxb.x[i] = nu.x[i];

			nxfh.x[i] = nu.x[i] + 2.0*dnx;
			VectorN cfh = evaluate(nxfh);
			nxfh.x[i] = nu.x[i];

			nxbh.x[i] = nu.x[i] - 2.0*dnx;
			VectorN cbh = evaluate(nxbh);
			nxbh.x[i] = nu.x[i];


			for (int j = 0; j < n; j++)
			{
				cx[j][i] = (8.0*(cf.x[j] - cb.x[j])-(cfh.x[j] - cbh.x[j]))/(12.0*dnx);
			}
		}

		Matrix Jac = new Matrix(cx, n, n);
		return Jac;
	}

	protected VectorN computeGc(Matrix j, VectorN fv)
	// computes the gradient needed by the linesearch
	{
		VectorN gc = new VectorN(n);
		Matrix dsq = d.times(d);      // square the scaling matrix
		Matrix jt = j.transpose();
		Matrix temp = jt.times(dsq);
		gc = temp.times(fv);
		return gc;
	}

	protected Matrix broydenUpdate(VectorN xc, VectorN xplus, VectorN fc, VectorN fplus, Matrix jacprev)
	// computes update to the Jacobian
	{
		Matrix jac = new Matrix(jacprev);
		VectorN s = xplus.plus(xc);
		VectorN den = this.d.times(s);
		double denmag = den.mag();
		double denom = denmag*denmag;
		double temp = 0.0;
		double tempi = 0.0;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				temp = temp + jac.A[i][j]*s.x[j];
			}
			tempi = fplus.x[i] - fc.x[i] - temp;
			if (Math.abs(tempi) >= 1.12E-16)
			{
				tempi = tempi / denom;
				for (int j = 0; j < n; j++)
				{
					jac.A[i][j] = jac.A[i][j]+tempi*s.x[j]*this.scales.x[j]*this.scales.x[j];
				}
			}
		}
		return jac;
	}

	protected VectorN linesearch(VectorN xin, VectorN p, VectorN g)
	// find the best value of k (correction scale factor) by backtracking if necessary
	{
		double init_slope = g.dotProduct(p);

		if (init_slope >= 0.0)
		{
			System.out.println("linesearch: roundoff problem");
		}

		double k = 1.0;
		double ktemp = 1.0;
		double kprev = 1.0;
		VectorN xc = new VectorN(xin);
		VectorN xplus = new VectorN(n);
		VectorN ftemp = new VectorN(n);
		VectorN ab = new VectorN(2);
		VectorN r = new VectorN(2);

		double fc = 0.0;
		double fplus = 0.0;
		double fpprev = 0.0;


		ftemp = evaluate(xin);
		fc = figureOfMerit(ftemp);
		double forig = fc;

		for (int it = 0; it < maxit; it++)
		{
			for (int i = 0; i < n; i++)
			{
				xplus.x[i] = xc.x[i] + k*p.x[i];
			}
			ftemp = evaluate(xplus);
			fplus = figureOfMerit(ftemp);

            //			System.out.println("fplus = "+fplus);

			double slope_crit = fc + this.alpha*k*init_slope;

            //			System.out.println("slope_crit = "+slope_crit);

			// check for acceptable step
			if (fplus <= slope_crit)
			{
            //				System.out.println("acceptable step, returning, k = "+k);
				return xplus;
			}
			else
			{
				if (k < kmin)  // if k goes too low, get out
				{
					System.out.println("linesearch: k went too low, no update");
                    //					System.exit(-1);
                    if (fplus < forig)
                    {
						return xplus;
					}
					else
					{
						return xin;
					}
				}

				// need to reduce k
				if (k == 1.0)   // first backtrack, use quadratic fit
				{
					ktemp = -1.0*init_slope/(2.0*(fplus - fc - init_slope));
                    //		            System.out.println("init_slope = "+init_slope);
                    //		            System.out.println("fplus = "+fplus);
                    //		            System.out.println("fc = "+fc);
                    //					System.out.println("first backtrack, ktemp = "+ktemp);
				}
				else  // subsequent backtrack, use cubic fit
				{
					Matrix q = new Matrix(2,2);
					double k2 = k*k;
					double kp2 = kprev*kprev;
					q.A[0][0] = 1.0/k2;
					q.A[0][1] = -1.0/kp2;
					q.A[1][0] = -1.0*kprev/k2;
					q.A[1][1] = k/kp2;
					r.x[0] = fplus - fc - k*init_slope;
					r.x[1] = fpprev - fc - kprev*init_slope;

					ab = q.times(r);
					ab = ab.times(1.0/(k - kprev));

					double disc = ab.x[1]*ab.x[1] - 3.0*ab.x[0]*init_slope;

					if (ab.x[0] == 0.0) // cubic is really a quadratic
					{
						ktemp = -1.0*init_slope/(2.0*ab.x[1]);
					}
					else  // solve the cubic equation for the new k
					{
						ktemp = (-1.0*ab.x[1]+Math.sqrt(disc))/(3.0*ab.x[0]);
					}

//					System.out.println("subsequent backtrack, ktemp = "+ktemp);

					if (ktemp > (0.5*k)) ktemp = 0.5 * k;
				}

				// store off k and fplus before modifying

				kprev = k;
				fpprev = fplus;
//				k = ktemp;

				if (ktemp <= 0.1*k)
				{
					k = 0.1*k;
				}
				else
				{
					k = ktemp;
				}

			}
		}
		System.out.println("linesearch maxit exceeded, k = "+k);
		return xin;
	}

	public VectorN solveIt(int max_iter, double usr_tol)
	{
		int iter = 0;

		VectorN xc = new VectorN(this.x);
        if (debug)		xc.print("initial unknowns");
		VectorN fvc = evaluate(xc);
	    if(debug)
            fvc.print("initial errors");
		double fc = figureOfMerit(fvc);
		double norm = Math.sqrt(2.0*fc);
		if(debug)
            System.out.println("initial norm = "+norm);
		VectorN gc = new VectorN(n);

		VectorN xplus = new VectorN(xc);
		VectorN fvplus = new VectorN(fvc);
		double fplus = fc;
		Householder house = new Householder(n, n);
//		Matrix jc = computeJacobian(xc);

//        VectorN pert = new VectorN(this.n);
//        int i;
//        for(i=0;i<n;i++)
//            pert.x[i] = 1.0E-15;
        // Generalize later
        /*
        pert.x[1] = 1.0E-15;
        pert.x[2] = 0.0;
        pert.x[3] = 1.0E-15;
        pert.x[4] = 1.0E-15;
        pert.x[5] = 0.0;
        pert.x[6] = 1.0E-15;
        pert.x[7] = 0.0;
        pert.x[8] = 0.0;
        */
		int prt = 0;
		int count = 0;
		int stuck = 0;

		while (fc > usr_tol && iter<max_iter)
		{
			iter++;
			prt++;
			count++;

			// wiggle the central difference pert step size


			// compute Jacobian
			Matrix jc = new Matrix(n);

//			if (eps > 1.0E-15)
//			{
//				eps = eps / 10.0;
//				jc = computeJacobian2(xc);
//			}
//			else
//			{
//				eps = 1.0E-07;
//				jc = computeJacobian2(xc);
//			}

//            count = 3;
//			if (count < 4)
//			{
//				switch (count)
//				{
//					case 1:
//					   eps = 1.0E-14;
//					   break;
////					case 2:
////					   eps = 1.0E-13;
////					   break;
////					case 3:
////					   eps = 1.0E-12;
////					   break;
//					case 2:
//					   eps = 1.0E-11;
//					   break;
////					case 5:
////					   eps = 1.0E-10;
////					   count = 0;
////					   break;
////					case 6:
////					   eps = 1.0E-9;
////					   count = 0;
////					   break;
//					case 3:
//					   eps = 1.0E-8;
//					   count = 0;
//					   break;
////					case 8:
////					   eps = 1.0E-7;
////					   count = 0;
////					   break;
//					default:
//					   eps = 1.0E-8;
//					   break;
//				   }
//
//				jc = computeJacobian2(xc);
//			}
//			else
//			{
				jc = computeJacobian2(xc);
				count = 0;
//			}


			if (debug) {
//				jc.print("Jacobian");
//				double jcdet = jc.det();
//				System.out.println("jcdet = "+jcdet);
			}
//            jc = broydenUpdate(xc, xplus, fvc, fvplus, jc);

			// find newton step direction
			VectorN p = new VectorN(house.compute(jc, fvc.x));
//            Matrix jc_inv=jc.inverse();
//            VectorN p=jc_inv.times(fvc);
			p = p.times(-1.0);

			// compute the gradient of f
			gc = computeGc(jc, fvc);

			// do the linesearch to find the new value of x
			xplus = linesearch(xc, p, gc);


//			VectorN xcorr = xplus.subtract(xc);

			// evaluate the results
			fvplus = evaluate(xplus);
			fplus = figureOfMerit(fvplus);
			norm = Math.sqrt(2.0*fplus);
//			xcorr.print("corrections");

            // check for stuck

//            if ((fc - fplus) <= machineps)
//            {
//				stuck++;
//			}
//			else
//			{
//				stuck = 0;
//			}
//
//			if (stuck > 3)
//			{
//				xplus = xplus.plus(pert);
//				stuck = 0;
//			}

			if(debug)
                //System.out.println(" iter = "+iter+" norm = "+norm+" eps = "+eps+" stuck = "+stuck);
				System.out.println(" iter = "+iter+" f = "+fplus+" x0 = [ "+xplus.x[0]+" "+xplus.x[1]+" "+xplus.x[2]+" "+xplus.x[3]);

			// store for the next iteration
			xc = xplus.copy();
			fc = fplus;
			fvc = fvplus.copy();

			// print the unknowns every 10 iterations

			if (prt > 5)
			{
				prt = 0;
				if(debug)
                {
				    xc.print("      nu.x");
				    fvc.print("errors");
				}
			}

		}

    	if(debug)
	    	xc.print("      nu.x");
		return xc;
	}

	public VectorN solveIt(){
		return solveIt(Integer.MAX_VALUE,term);
	}
}

