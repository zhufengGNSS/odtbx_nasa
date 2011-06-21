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

package jat.forces.gravity;
import jat.forces.ForceModel;
import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;
import jat.spacecraft.Spacecraft;
import jat.spacetime.BodyRef;
import jat.spacetime.ReferenceFrame;
import jat.spacetime.ReferenceFrameTranslater;
import jat.spacetime.Time;

/** <P>
 * The SphericalHarmonicGravity class computes the acceleration due to gravity on a satellite
 * using a spherical harmonic gravity field. This class will be the base class for various
 * spherical harmonic models such as JGM-3, etc. The subclasses will supply the coefficients.
 * Reference: Montenbruck, section 3.2. This is basically Montenbruck's code translated to Java.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

abstract public class SphericalHarmonicGravity implements ForceModel {

	private static final long serialVersionUID = 1L;

    /** GM in m^3/s^2.
     */
    protected double GM;
    /** Radius of planet in m.
     */
    protected double R_ref;
    /** Acceleration scalefactor
     */
    protected double accel_sf;
    
    /** Unnormalized harmonic coefficients.
     */
    protected double[][] CS;
    /** Maximum degree (20).
     */
    protected int n_max;
    /** Maximum order (20).
     */
    protected int m_max;
    /** Harmonic function V.
     */
    protected double[][] V;            // Harmonic functions
    /** Harmonic function W.
     */
    protected double[][] W;            // work array (0..n_max+1,0..n_max+1)

    /** Degree of gravity model desired.
     */
    protected int n_desired;

    /** Order of gravity model desired.
     */
    protected int m_desired;
    
    /** The reference frame in which the gravity is computed. */
    private ReferenceFrame bodyFixedRef;

    /** The user must supply a method to initialize the gravity model. The initialize() method must set the values of GM, R_ref, n_max, m_max and CS to be used.
     */
    abstract public void initialize();

    /** Initialize the value of GM and R_ref.
     * @param gm the planetary gravity constant.
     * @param rref the planetary R_ref value.
     */
    protected void initializeGM(double gm, double rref){
        this.GM = gm;
        this.R_ref = rref;
        this.accel_sf = GM/(R_ref*R_ref);
    }


    /** Initialize the spherical harmonic coefficients.
     * @param nmax The maximum order of the model.
     * @param mmax The maximum degree of the model.
     * @param cs The unnormalized harmonic coefficients.
     */
    protected void initializeCS(int nmax, int mmax, double[][] cs){
        this.n_max = nmax;
        this.m_max = mmax;
        this.CS = new double[nmax+1][mmax+1];
        for (int i = 0; i < this.n_max+1; i++){
            for (int j = 0; j < this.m_max+1; j++){
                this.CS[i][j] = cs[i][j];
            }
        }
    }

    /** Constructor.
     * @param n Desired degree.
     * @param m Desired order.
     */
    public SphericalHarmonicGravity(int n, int m, 
        ReferenceFrame bodyFixedFrame) {
      
        this.n_desired = n;
        this.m_desired = m;
        this.bodyFixedRef = bodyFixedFrame;
        V =  new double[this.n_desired+2][this.n_desired+2];
        W = new double[this.n_desired+2][this.n_desired+2];
    }

    /** Evaluates the two harmonic functions V and W.
     * @param r_bf fixed body position vector.
     */
    private void computeVW(VectorN r_bf) {

        // Auxiliary quantities
        double r_sqr =  r_bf.dotProduct(r_bf);               // Square of distance
        double rho   =  R_ref*R_ref / r_sqr;

        double x0 = R_ref * r_bf.x[0] / r_sqr;          // Normalized
        double y0 = R_ref * r_bf.x[1] / r_sqr;          // coordinates
        double z0 = R_ref * r_bf.x[2] / r_sqr;

        // Calculate zonal terms V(n,0); set W(n,0)=0.0
        double temp = R_ref / Math.sqrt(r_sqr);
        V[0][0] = temp;
        W[0][0] = 0.0;

        temp = z0 * temp;
        V[1][0] = temp;
        W[1][0] = 0.0;

        for (int n=2; n<=n_desired+1; n++) {
            temp = ( (2*n-1) * z0 * V[n-1][0] - (n-1) * rho * V[n-2][0] ) / n;
            V[n][0] = temp;
            W[n][0] =0.0;
        }

        // Calculate tesseral and sectorial terms
        for (int m=1; m<=m_desired+1; m++) {

            // Calculate V(m,m) .. V(n_max+1,m)

            V[m][m] = (2*m-1) * ( x0*V[m-1][m-1] - y0*W[m-1][m-1] );
            W[m][m] = (2*m-1) * ( x0*W[m-1][m-1] + y0*V[m-1][m-1] );

            if (m<=n_desired) {
                V[m+1][m] = (2*m+1) * z0 * V[m][m];
                W[m+1][m] = (2*m+1) * z0 * W[m][m];
            }

            for (int n=m+2; n<=n_desired+1; n++) {
                V[n][m] = ( (2*n-1)*z0*V[n-1][m] - (n+m-1)*rho*V[n-2][m] ) / (n-m);
                W[n][m] = ( (2*n-1)*z0*W[n-1][m] - (n+m-1)*rho*W[n-2][m] ) / (n-m);
            }
        }
    }

    /** Computes the acceleration due to gravity in m/s^2.
     * @param r ECI position vector.
     * @param E ECI to ECEF transformation matrix.
     * @return ECI acceleration in m/s^2.
     */
    public VectorN gravity(VectorN r, Matrix E){
        // Local variables

        // Evaluate harmonic functions
        // Rotate from ECI to ECEF
        VectorN r_bf = E.times(r); 
        VectorN a_bf = bodyFixedGravity(r_bf);

        // Inertial acceleration
        VectorN out = E.transpose().times(a_bf);
        
        //out.print("acceleration");
        return out;
    }
    
    /** Computes the acceleration due to gravity in m/s^2.
     * @param r_bf body fixed position vector.
     * @return body-fixed acceleration in m/s^2.
     */
    public VectorN bodyFixedGravity(VectorN r_bf){
        // Local variables

        int     n,m;                           // Loop counters
        double  Fac;                           // Auxiliary quantities
        double  C,S;                           // Gravitational coefficients

        // Evaluate harmonic functions
        computeVW(r_bf);

        // Calculate accelerations ax,ay,az
        double ax = 0.0;
        double ay = 0.0;
        double az = 0.0;
        
        for (m=0; m<=m_desired; m++)
            for (n=m; n<=n_desired ; n++)
                if (m==0) {
                    C = CS[n][0];   // = C_n,0
                    ax -=       C * V[n+1][1];
                    ay -=       C * W[n+1][1];
                    az -= (n+1)*C * V[n+1][0];
                }
                else {
                    C = CS[n][m];   // = C_n,m
                    S = CS[m-1][n]; // = S_n,m
                    Fac = 0.5 * (n-m+1) * (n-m+2);
                    ax += 0.5*(-C*V[n+1][m+1] - S*W[n+1][m+1]) + Fac*(C*V[n+1][m-1] + S*W[n+1][m-1]);
                    ay += 0.5*(-C*W[n+1][m+1] + S*V[n+1][m+1]) + Fac*(-C*W[n+1][m-1] + S*V[n+1][m-1]);
                    az += (n-m+1)*(-C*V[n+1][m] - S*W[n+1][m]);
                }

        // Body-fixed acceleration
        VectorN a_bf =  new VectorN(ax, ay, az);
        a_bf = a_bf.times(this.accel_sf);
//        a_bf.print("body fixed gravity");
        
        return a_bf;
    }
    
    /** Computes the acceleration due to non-spherical gravity in m/s^2.
     * For use in Encke integration.
     * @param r_bf body fixed position vector.
     * @return body-fixed acceleration in m/s^2.
     */
    public VectorN bodyFixedPerturbation(VectorN r_bf){
        // Local variables

        int     n,m;                           // Loop counters
        double  Fac;                           // Auxiliary quantities
        double  C,S;                           // Gravitational coefficients

        // Evaluate harmonic functions
        computeVW(r_bf);

        // Calculate accelerations ax,ay,az
        double ax = 0.0;
        double ay = 0.0;
        double az = 0.0;
        
        VectorN a_twobody = new VectorN(3);

        for (m=0; m<=m_desired; m++)
            for (n=m; n<=n_desired ; n++)
                if (m==0) {
                    C = CS[n][0];   // = C_n,0
                    ax -=       C * V[n+1][1];
                    ay -=       C * W[n+1][1];
                    az -= (n+1)*C * V[n+1][0];
                    if (n == 1){
                    	a_twobody.set(0, ax);
                    	a_twobody.set(1, ay);
                    	a_twobody.set(2, az);
                    }
                }
                else {
                    C = CS[n][m];   // = C_n,m
                    S = CS[m-1][n]; // = S_n,m
                    Fac = 0.5 * (n-m+1) * (n-m+2);
                    ax += 0.5*(-C*V[n+1][m+1] - S*W[n+1][m+1]) + Fac*(C*V[n+1][m-1] + S*W[n+1][m-1]);
                    ay += 0.5*(-C*W[n+1][m+1] + S*V[n+1][m+1]) + Fac*(-C*W[n+1][m-1] + S*V[n+1][m-1]);
                    az += (n-m+1)*(-C*V[n+1][m] - S*W[n+1][m]);
                }

        // Total Body-fixed acceleration
        VectorN a_bf =  new VectorN(ax, ay, az);
        a_bf = a_bf.times(this.accel_sf);
//        a_bf.print("body fixed gravity");
        
        // compute the two body acceleration
        VectorN tb = a_twobody.times(this.accel_sf);
//        tb.print("twobody grav");

        // perturbation = total accel - two body accel
        return a_bf.minus(tb);
    }    
    
    /** Implemented from the ForceModel interface
     * @param t Time reference object
     * @param bRef Earth reference object
     * @param sc Spacecraft parameters
     */
    public VectorN acceleration(Time t, BodyRef bRef, Spacecraft sc) {
      
      ReferenceFrameTranslater xlater = 
        new ReferenceFrameTranslater(bRef, bodyFixedRef, t);
      VectorN r_bf = xlater.translatePoint(sc.r());
      VectorN a_bf = bodyFixedGravity(r_bf);
      return xlater.transformDirectionBack(a_bf);
    }
    
    /** Computes the partial derivative of gravity with respect to position.
     * @return ECI gravity gradient matrix.
     * @param r ECI position vector.
     * @param E ECI to ECEF transformation matrix.
     */
    public Matrix gravityGradient(VectorN r, Matrix E){

        // Evaluate harmonic functions
        VectorN r_bf = E.times(r);
        computeVW(r_bf);

        double xx = 0.0;
        double xy = 0.0;
        double xz = 0.0;
        double yy = 0.0;
        double yz = 0.0;
        double zz = 0.0;
        double Fac = 0.0;
        double f1 = 0.0;
        double f2 = 0.0;
        int m;
        int n;
        double C;
        double S;
        Matrix out = new Matrix(3,3);

        for (m=0; m<=m_desired; m++) {
            for (n=m; n<=n_desired ; n++) {
                Fac = (n-m+2)*(n-m+1);
                C = CS[n][m];
                S = CS[m-1][n];
                zz += Fac*(C*V[n+2][m] + S*W[n+2][m]);
                if (m==0) {
                    C = CS[n][0];   // = C_n,0
                    Fac = (n+2)*(n+1);
                    xx += 0.5 * (C*V[n+2][2] - Fac*C*V[n+2][0]);
                    xy += 0.5 * C * W[n+2][2];
                    Fac = n + 1;
                    xz += Fac * C * V[n+2][1];
                    yz += Fac * C * W[n+2][1];
                }
                if (m > 0){
                    C = CS[n][m];
                    S = CS[m-1][n];
                    f1 = 0.5*(n-m+1);
                    f2 = (n-m+3)*(n-m+2)*f1;
                    xz += f1*(C*V[n+2][m+1]+S*W[n+2][m+1])-f2*(C*V[n+2][m-1]+S*W[n+2][m-1]);
                    yz += f1*(C*W[n+2][m+1]+S*V[n+2][m+1])+f2*(C*W[n+2][m-1]-S*V[n+2][m-1]);
                    if (m == 1){
                        Fac = (n+1)*n;
                        xx += 0.25*(C*V[n+2][3]+S*W[n+2][3]-Fac*(3.0*C*V[n+2][1]+S*W[n+2][1]));
                        xy += 0.25*(C*W[n+2][3]-S*V[n+2][3]-Fac*(C*W[n+2][1]+S*V[n+2][1]));
                    }
                    if (m > 1) {
                        f1 = 2.0*(n-m+2)*(n-m+1);
                        f2 = (n-m+4)*(n-m+3)*f1*0.5;
                        xx += 0.25*(C*V[n+2][m+2]+S*W[n+2][m+2]-f1*(C*V[n+2][m]+S*W[n+2][m])+f2*(C*V[n+2][m-2]+S*W[n+2][m-2]));
                        xy += 0.25*(C*W[n+2][m+2]-S*V[n+2][m+2]+f2*(-C*W[n+2][m-2]+S*V[n+2][m-2]));
                    }

                }
            }
            yy = -xx - zz;
            out.set(0, 0, xx);
            out.set(0, 1, xy);
            out.set(0, 2, xz);
            out.set(1, 0, xy);
            out.set(1, 1, yy);
            out.set(1, 2, yz);
            out.set(2, 0, xz);
            out.set(2, 1, yz);
            out.set(2, 2, zz);
        }

        out = out.times(this.accel_sf);

        // Rotate to ECI
        Matrix Etrans = E.transpose();
        out = Etrans.times(out.times(E));
        return out;
    }
    
    public double getGM(){
    	return this.GM;
    }

    /** Print out the gravity model parameters.
     */
    public void printParameters(){
        System.out.println("GM = "+this.GM);
        System.out.println("R_ref = "+this.R_ref);
        System.out.println("n_max = "+this.n_max+" m_max = "+this.m_max);
        System.out.println("n_desired = "+this.n_desired+" m_desired = "+this.m_desired);
        Matrix cs = new Matrix(this.CS);
        cs.print("CS");
    }
}
