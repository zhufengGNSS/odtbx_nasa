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
import java.io.*;

/**
 * <P>
 * The GaussianVector Class provides the means to create vectors of numbers which are samples of Gaussian
 * distributions. The algorithms come from Chapter 7 of Numerical Recipes.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class GaussianVector extends VectorN {

	private static final long serialVersionUID = 1L;
    
    /** Seed for random number generator.
     */
    protected long idum;
    /** VectorN containing the mean for each component.
     */
    public VectorN mean;
    /** VectorN containing the sigman for each component.
     */
    public VectorN sigma;
    
    private RandomNumber rn;
    /** Creates a new instance of GaussianVector (Vector3 with zero means). */
    public GaussianVector() {
        super(3);
        this.mean = new VectorN(3);
        this.sigma = new VectorN(1.0, 1.0, 1.0);
        this.idum = -1;
        this.rn = new RandomNumber(this.idum);
        fillVector();
    }
        
    /** Create a GaussianVector given a VectorN containing the means and sigmas.
     * @param mn VectorN containing the mean values.
     * @param sig VectorN containing the sigmas values.
     */
    public GaussianVector(VectorN mn, VectorN sig){
        super(mn.length);
        sig.checkVectorDimensions(mn.length);
        this.mean = new VectorN(mn);
        this.sigma = new VectorN(sig);
        this.idum = -1;
        this.rn = new RandomNumber(this.idum);
        fillVector();
    }
    
    /** Create a GaussianVector given a VectorN containing the means, sigmas and
     * a seed.
     * @param mn VectorN containing the mean values.
     * @param sig VectorN containing the sigmas values.
     * @param seed Seed for the random number generator.
     */
    public GaussianVector(VectorN mn, VectorN sig, long seed){
        super(mn.length);
        sig.checkVectorDimensions(mn.length);
        this.mean = new VectorN(mn);
        this.sigma = new VectorN(sig);
        if (seed > 0) seed = -seed;
        this.idum = seed;
        this.rn = new RandomNumber(this.idum);
        fillVector();
    }
    
    /** Fills the vector components with random numbers.
     */
    public void fillVector(){
        int n = mean.length;
        for (int i = 0; i < n; i++) {
            this.x[i] = rn.normal()*sigma.x[i] + mean.x[i];
        }
    }
    
    /** Fills the vector components with the next set of random numbers.
     */
    public void nextSet(){
        this.fillVector();
    }
    
    /** Method to test the class.
     * @param args Program arguments.
     * @throws IOException If there is a problem with the output file.
     */    
    public static void main(java.lang.String args[]) throws IOException {
        
        FileOutputStream outf = new FileOutputStream("gaussian.txt");
        PrintWriter pw = new PrintWriter(outf);
        PrintWriter sys = new PrintWriter(System.out, true);
        VectorN means = new VectorN(1.0, 2.0, 3.0);
        VectorN sigmas = new VectorN(1.0, 2.0, 3.0);
        GaussianVector x = new GaussianVector(means, sigmas);
        
        for (int i = 0; i < 5000; i++) {
            x.fillVector();
            x.print(pw);
            x.print(sys);
        }
        pw.close();
        outf.close();
    }
}
