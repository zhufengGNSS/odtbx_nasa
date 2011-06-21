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

package jat.matvec.data.arrayTools;

   /** Quick Sort algoritm.
   <P>
   Allows to sort a column quickly,
   Using a generic version of C.A.R Hoare's Quick Sort algorithm.
   <P>
   */

public class Sort {

/* ------------------------
   Class variables
 * ------------------------ */

   /** Array for internal storage of the matrix to sort.
   */

  private double[] A;

   /** Array for internal storage of the order.
   */
  private int[] order;


/* ------------------------
   Constructors
 * ------------------------ */

   /** Construct an ascending order.
   @param a    Array to sort.
   */

  public Sort(double[] a) {
    A = a;
    order = new int[a.length];
    for (int i = 0; i < a.length; i++)
        order[i] = i;

    sort(A);
  }

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Get the ascending order of one line.
   @param i    Line number.
   @return     Ascending order of the line.
   */

  public int getOrder(int i) {
    return order[i];
  }

/* ------------------------
   Private Methods
 * ------------------------ */

     /** This is a generic version of C.A.R Hoare's Quick Sort
    * algorithm.  This will handle arrays that are already
    * sorted, and arrays with duplicate keys.<BR>
    *
    * If you think of a one dimensional array as going from
    * the lowest index on the left to the highest index on the right
    * then the parameters to this function are lowest index or
    * left and highest index or right.  The first time you call
    * this function it will be with the parameters 0, a.length - 1.
    *
    * @param a       A double array.
    * @param lo0       Int.
    * @param hi0       Int.
    */

   private void QuickSort(double a[],int lo0, int hi0) {

      int lo = lo0;
      int hi = hi0;
      double mid;

      if ( hi0 > lo0)
      {

         /* Arbitrarily establishing partition element as the midpoint of
          * the array.
          */
         mid = a[ ( lo0 + hi0 ) / 2 ];

         // loop through the array until indices cross
         while( lo <= hi )
         {
            /* find the first element that is greater than or equal to
             * the partition element starting from the left Index.
             */
	     while( ( lo < hi0 ) && ( a[lo] < mid ))
		 ++lo;

            /* find an element that is smaller than or equal to
             * the partition element starting from the right Index.
             */
	     while( ( hi > lo0 ) && ( a[hi] > mid ))
		 --hi;

            // if the indexes have not crossed, swap
            if( lo <= hi )
            {
               swap(a, lo, hi);
               ++lo;
               --hi;
            }
         }

         /* If the right index has not reached the left side of array
          * must now sort the left partition.
          */
         if( lo0 < hi )
            QuickSort( a, lo0, hi );

         /* If the left index has not reached the right side of array
          * must now sort the right partition.
          */
         if( lo < hi0 )
            QuickSort( a, lo, hi0 );

      }
   }

   /** Swap two positions.
   @param a    Array.
   @param i    Line number.
   @param j    Line number.
   */

   private void swap(double a[], int i, int j) {
      double T;
      T = a[i];
      a[i] = a[j];
      a[j] = T;
      int t;
      t = order[i];
      order[i] = order[j];
      order[j] = t;
   }

   private void sort(double[] a) {
    QuickSort(a, 0, a.length-1);
   }
}