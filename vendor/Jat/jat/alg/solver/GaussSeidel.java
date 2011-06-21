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
 * <P>
 *  Solve the problem Ax=b through use of Gauss-Seidel method
 *  as described in "Numerical Analysis" by Richard Burden and J. Douglas Faires
 *
 *
 * @author Tobias Berthold
 * Date        :   2-3-2003
 * @version 1.0
 */


package jat.alg.solver;

import jat.math.*;
import jat.matvec.data.*;

public class GaussSeidel
{
    Matrix A;
    VectorN x,x0,xold,b;
    int n,maxIter;
    double tol;
    
    /**
	 * @param A Matrix
	 * @param b right hand side
	 * @param x0 initial guess
	 * @param maxIter maximum iterations
	 * @param tol error tolerance
	 */
	// Constructor
    public GaussSeidel (Matrix A, VectorN b, VectorN x0, int maxIter, double tol)
    {
        this.A=A;
        this.x0=x0;
        this.b=b;
        this.maxIter=maxIter;
        this.tol=tol;
        this.n=x0.length;
    }
    
    /**
	 * @return solution vector
	 */
	public VectorN iterate ()
    {
        int k=0;
        double err=1.;
        
        // Possible checks: No zeros on diagonal
        
        x=new VectorN (x0);
        xold=new VectorN (x0);
               
        while(k<maxIter && tol<=err )
        {
            k++;
            int i, j;
            // Get new x
            for(i=0; i<n; i++)
            {
                double sum = 0;
                // Solve for current vector component
                for(j=0; j<n; j++)
                {
                    if(j!=i)
                        sum += A.A[i][j]*x.x[j];
                }
                x.x[i] = (-sum + b.x[i])/A.A[i][i];
            }
            x.print ("x");
            
            // Check error tolerance
            err=MathUtils.abs (x.mag ()-xold.mag ());
            System.out.println ("Error: "+err);
            System.out.println (" "+x.mag ());
            System.out.println (" "+xold.mag ());
            xold=x.copy ();
            
            /*
            xOld = xOld.minus(x);
            if( xOld.mag() / x.mag() <= tol) return x;
             */
            //xOld = x;
        }
        return x;
    }
}
