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

package jat.matvec.data;

/**
 * <P>
 * The UniformVector Class provides the means to create vectors of numbers which are samples of uniform
 * distributions. The algorithms come from Chapter 7 of Numerical Recipes.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class UniformVector extends VectorN {

	private static final long serialVersionUID = 1L;
    
    /** Seed for random number generator.
     */    
    protected long idum;

    /** VectorN containing the min for each component.
     */    
    public VectorN min;

    /** VectorN containing the max for each component.
     */    
    public VectorN max;
    
    private RandomNumber rn;    
    
    /** Creates a new instance of UniformVector (Vector3 with zero means). */
    public UniformVector() {
        super(3);
        this.min = new VectorN(3);
        this.max = new VectorN(1.0, 1.0, 1.0);
        this.idum = -1;
        this.rn = new RandomNumber(this.idum);
        fillVector();
    }
    
    /** Create a UniformVector given a VectorN containing the mins and another
     * containing the max values.
     * @param mn VectorN containing the min values.
     * @param mx VectorN containing the max values.
     */
    public UniformVector(VectorN mn, VectorN mx){
        super(mn.length);
        mx.checkVectorDimensions(mn);
        this.min = new VectorN(mn);
        this.max = new VectorN(mx);
        this.idum = -1;
        this.rn = new RandomNumber(this.idum);        
        fillVector();
    }
    
    /** Create a UniformVector given a VectorN containing the mins and another
     * containing the max values.
     * @param mn VectorN containing the min values.
     * @param mx VectorN containing the max values.
     * @param seed Seed for random number generator.
     */
    public UniformVector(VectorN mn, VectorN mx, long seed){
        super(mn.length);
        mx.checkVectorDimensions(mn);
        this.min = new VectorN(mn);
        this.max = new VectorN(mx);
        if (seed > 0) seed = -seed;
        this.idum = seed;
        this.rn = new RandomNumber(this.idum);        
        fillVector();
    }
    
    /** Fills the vector components with random numbers.
     */    
    public void fillVector(){
        int n = min.length;
        for (int i = 0; i < n; i++) {
            this.x[i] = min.x[i] + (max.x[i] - min.x[i])*rn.uniformDeviate();
        }
    }

    /** Fills the vector components with the next set of random numbers.
     */
    public void nextSet(){
        this.fillVector();
    }
    
}
