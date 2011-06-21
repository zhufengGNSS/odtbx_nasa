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

package jat.test.GaussSeidel_solve_test;

import jat.alg.solver.*;
import jat.matvec.data.*;

// Book Burden Faires
// Chapter 7.3 example 1

public class GaussSeidel_test
{    
    public static void main(String[] args)
    {
        int n=4;
        int maxIter=10;
        double tol=1.e-10;
        
        double[] in_A=
        {
            10,-1,2,0,
            -1,11,-1,3,
            2,-1,10,-1,
            0,3,-1,8
        };
        double[] in_b={6,25,-11,15};
        double[] in_x0={0,0,0,0};   // Initial guess
        
        Matrix A=new Matrix(in_A,4);
        VectorN x0=new VectorN(n);
        VectorN b=new VectorN(in_b);
        
        A.print("A matrix");
        b.print("b");
        
        GaussSeidel GS=new GaussSeidel(A,b,x0,maxIter,tol);
        VectorN x=GS.iterate();
        
        System.out.println("Solution for x");
        x.print("x");

        /*
                Vector y (n);
                matVecMult (A,x,y);
                y -= b;
                cout << "The error norm is: " << maxNorm (y) << endl;
                cout << "The error is: " << endl;
                cout << y << endl;
                return 0;
        }
         */
    }
}
