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
 * The Quaternion Class provides the fundamental operations of quaternions.
 * Quaternions are used to transform a Vector3 from Reference Frame A
 * to Reference Frame B.
 *
 * The following quaternion definition is used:
 *      q.x[0] = e.x[0] * sin(phi/2)
 *      q.x[1] = e.x[1] * sin(phi/2)
 *      q.x[2] = e.x[2] * sin(phi/2)
 *      q.x[3] = cos(phi/2)
 *
 * Reference: Farrell, Jay A. and Matthew Barth, The Global Positioning System &
 * Inertial Navigation, pp. 39 - 42.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class Quaternion extends VectorN {

	private static final long serialVersionUID = 1L;
    
/* ------------------------
   Constructors
 * ------------------------ */
    
    /** Creates new Quaternion */
    public Quaternion() {
        super(4);
    }
    
    /** Construct a quaternion from a double array.
     * @param x    Array containing the 4 values.
     */
    
    public Quaternion(double[] x) {
        super(x);
    }
    
    /** Construct a quaternion from a Vector4.
     * @param r    Vector4 containing the 4 values.
     */
    
    public Quaternion(VectorN r) {
        super(r);
        r.checkVectorDimensions(4);
    }
    
    /** Construct a quaternion from 4 doubles.
     * @param q1    1st component of the quaternion.
     * @param q2    2nd component of the quaternion.
     * @param q3    3rd component of the quaternion.
     * @param q4    4th component of the quaternion.
     */
    
    public Quaternion(double q1, double q2, double q3, double q4) {
        super(4);
        x[0] = q1;
        x[1] = q2;
        x[2] = q3;
        x[3] = q4;
    }
    
    /** Construct a quaternion from a Vector3 and an angle.
     * @param e    Vector3 axis of rotation (unit vector), expressed in frame A.
     * @param phi  Rotation angle in radians.
     */
    
    public Quaternion(VectorN e, double phi) {
        super(4);
        double cos_phi = Math.cos(phi/2.0);
        double sin_phi = Math.sin(phi/2.0);
        x[0] = e.x[0] * sin_phi;
        x[1] = e.x[1] * sin_phi;
        x[2] = e.x[2] * sin_phi;
        x[3] = cos_phi;
    }
    
    /** Construct a quaternion from a 3x3 Direction Cosine Matrix.
     * @param c    Direction Cosine Matrix from A to B.
     */
    
    public Quaternion(Matrix c) {
        super(4);
        c.checkMatrixDimensions(3, 3);
        double[][] r = c.getArrayCopy();
        x[3] = 0.5 * Math.sqrt(1.0 + r[0][0] + r[1][1] + r[2][2]);
        double den = 4.0 * x[3];
        x[0] = (r[2][1] - r[1][2])/den;
        x[1] = (r[0][2] - r[2][0])/den;
        x[2] = (r[1][0] - r[0][1])/den;
    }
/* ------------------------
   Public Methods
 * ------------------------ */
    
    /** Compute the conjugate of this quaternion.
     * @return     Conjugate of this quaternion.
     */
    
    public Quaternion conjugate() {
        double q1 = -1.0 * this.x[0];
        double q2 = -1.0 * this.x[1];
        double q3 = - 1.0 * this.x[2];
        double q4 = this.x[3];
        Quaternion out = new Quaternion(q1, q2, q3, q4);
        return out;
    }
    
    /** Compute the 4x4 Matrix version of this quaternion.
     * @return     4x4 Matrix version of this quaternion.
     */
    
    public Matrix yMatrix() {
        double[][] temp = new double[4][4];
        double q1 = this.x[0];
        double q2 = this.x[1];
        double q3 = this.x[2];
        double q4 = this.x[3];
        
        temp[0][0] = q4;
        temp[0][1] = -1.0 * q3;
        temp[0][2] = q2;
        temp[0][3] = q1;
        
        temp[1][0] = q3;
        temp[1][1] = q4;
        temp[1][2] = -1.0 * q1;
        temp[1][3] = q2;
        
        temp[2][0] = -1.0 * q2;
        temp[2][1] = q1;
        temp[2][2] = q4;
        temp[2][3] = q3;
        
        temp[3][0] = -1.0 * q1;
        temp[3][1] = -1.0 * q2;
        temp[3][2] = -1.0 * q3;
        temp[3][3] = q4;
        Matrix out = new Matrix(temp);
        return out;
    }

    /** Compute the 4x4 Matrix version of this quaternion. R = [2Q|-q]
     * @return     4x4 Matrix version of this quaternion.
     */
    
    public Matrix rMatrix() {
        double[][] temp = new double[4][4];
        double q1 = this.x[0];
        double q2 = this.x[1];
        double q3 = this.x[2];
        double q4 = this.x[3];
        
        temp[0][0] = q4;
        temp[0][1] = -1.0 * q3;
        temp[0][2] = q2;
        temp[0][3] = -1.0 * q1;
        
        temp[1][0] = q3;
        temp[1][1] = q4;
        temp[1][2] = -1.0 * q1;
        temp[1][3] = -1.0 * q2;
        
        temp[2][0] = -1.0 * q2;
        temp[2][1] = q1;
        temp[2][2] = q4;
        temp[2][3] = - 1.0 * q3;
        
        temp[3][0] = -1.0 * q1;
        temp[3][1] = -1.0 * q2;
        temp[3][2] = -1.0 * q3;
        temp[3][3] = -1.0 * q4;
        Matrix out = new Matrix(temp);
        return out;
    }    
    
    /** Compute the 4x3 Matrix version of this quaternion.
     * @return     4x3 Matrix version of this quaternion.
     */
    
    public Matrix qMatrix() {
        double[][] temp = new double[4][3];
        double q1 = this.x[0];
        double q2 = this.x[1];
        double q3 = this.x[2];
        double q4 = this.x[3];
        
        temp[0][0] = q4;
        temp[0][1] = -1.0 * q3;
        temp[0][2] = q2;
        
        temp[1][0] = q3;
        temp[1][1] = q4;
        temp[1][2] = -1.0 * q1;
        
        temp[2][0] = -1.0 * q2;
        temp[2][1] = q1;
        temp[2][2] = q4;
        
        temp[3][0] = -1.0 * q1;
        temp[3][1] = -1.0 * q2;
        temp[3][2] = -1.0 * q3;
        Matrix out = new Matrix(temp);
        out = out.divide(2.0);
        return out;
    }
    
    
    /** Compute the product of two quaternions.
     * @param q2   Quaternion to multiply
     * @return     q1 x q2.
     */
    public Quaternion times(Quaternion q2) {
        Matrix temp = this.yMatrix();
        VectorN q3 = temp.times(q2);
        Quaternion out = new Quaternion(q3);
        return out;
    }
    
    /** Compute the Direction Cosine Matrix from this quaternion.
     * @return     Direction Cosine Matrix from A to B.
     */
    
    public Matrix quat2DCM() {
        double[][] r = new double[3][3];
        double q1 = this.x[0];
        double q2 = this.x[1];
        double q3 = this.x[2];
        double q4 = this.x[3];
        
        double p1 = q1 * q1;
        double p2 = q2 * q2;
        double p3 = q3 * q3;
        double p4 = q4 * q4;
        double p5 = p2 + p3;
        double p6 = p5 + p4 + p1;
        if (p6 != 0.0) p6 = 2.0 / p6;
        
        r[0][0] = 1.0 - p6 * p5;
        r[1][1] = 1.0 - p6 * (p1 + p3);
        r[2][2] = 1.0 - p6 * (p1 + p2);
        
        p1 = p6 * q1;
        p2 = p6 * q2;
        p5 = p6 * q3 * q4;
        p6 = p1 * q2;
        
        r[0][1] = p6 - p5;
        r[1][0] = p6 + p5;
        
        p5 = p2 * q4;
        p6 = p1 * q3;
        
        r[0][2] = p6 + p5;
        r[2][0] = p6 - p5;
        
        p5 = p1 * q4;
        p6 = p2 * q3;
        
        r[1][2] = p6 - p5;
        r[2][1] = p6 + p5;
        
        Matrix out = new Matrix(r);
        return out;
    }
    
    /** Transform a Vector3 from Frame A to Frame B using this quaternion.
     * @param xin  Vector3 expressed in Frame A.
     * @return     Vector3 expressed in Frame B.
     */
    
    public VectorN transform(VectorN xin) {
        xin.checkVectorDimensions(3);
        Matrix c = this.quat2DCM();
        VectorN out = c.times(xin);
        return out;
    }
    
    /** Compute the perturbation quaternion from 3 tilts, performed on the true
     *  quaternion.
     * @param tilts  Vector3 containing 3 tilt angles in radians.
     * @return      perturbation quaternion.
     */

    public Quaternion tilt2quat(VectorN tilts) {    
        tilts.checkVectorDimensions(3);
        if (tilts.mag() == 0.0) {
            Quaternion out = new Quaternion();
            return out;
        }
        else {
            Matrix Q = this.qMatrix();
            Q = Q.times(-1.0);
            VectorN temp = Q.times(tilts);
            Quaternion out = new Quaternion(temp);
            return out;
        }
    }
        
//    public Quaternion tilt2quat(VectorN tilts) {
//        tilts.checkVectorDimensions(3);
        
//        if (tilts.mag() == 0.0) {
//            Quaternion out = new Quaternion();
//            return out;
//        }
//        else {
//            Matrix ctrue = this.quat2DCM();
 //           double ctdet = ctrue.det();
//            System.out.println("det of ctrue = "+ctdet);
//            RotationMatrix b = new RotationMatrix(tilts.x[0], tilts.x[1], tilts.x[2]);
//            b.print("b matrix");
//            double bdet = b.det();
//            System.out.println("det of b = "+bdet);
//            Matrix btrans = b.transpose();
 //           double btdet = btrans.det();
//            System.out.println("det of btrans = "+btdet);
//            btrans = btrans.divide(btdet);
//            btdet = btrans.det();
//            System.out.println("det of btrans = "+btdet);            
//            Matrix chat = btrans.times(ctrue);
//            ctrue.print("c true");
//            chat.print("chat");
            
//            this.print("q true");
            
//            Quaternion qhat = new Quaternion(chat);
//            qhat.print("qhat");
//            VectorN temp = this.minus(qhat);
//            Quaternion out = new Quaternion(temp);
//            return out;
//        }
//    }
    
    /** Compute the 3 tilts from a perturbation quaternion. Implements equation 29
     *  from Friedland, Bernard, "Analysis of Strapdown Navigation Using Quaternions".
     * @param q     Base quaternion (not the perturbation quaternion)
     * @return      VectorN containing the tilts in radians.
     */
    
    public VectorN tilts(Quaternion q) {
        // form -4Q^T
        Matrix Q = q.qMatrix();
        Matrix Qt = Q.transpose();
        Qt = Qt.times(-4.0);
//        System.out.println("tilts function");
//        double dqmag = this.mag();
//        double qmag = q.mag();
//        double Qtdet = Qt.det();
//        System.out.println("dqmag = "+dqmag+" qmag = "+qmag);
        
        double[] temp = Qt.times(this.x);
        VectorN out = new VectorN(temp);
        return out;
    }
    
    /** Compute the big omega matrix from an angular rate vector.
     * Used to propagate quaternions: qdot = omega * q
     * @param w angular rate vector
     * @return angular rate matrix
     */
	public static Matrix omega(VectorN w) {
		w.checkVectorDimensions(3);
        Matrix out = new Matrix(4,4);
        double w1 = 0.5 * w.x[0];
        double w2 = 0.5 * w.x[1];
        double w3 = 0.5 * w.x[2];
		out.set(0,1, w3);
		out.set(0,2, -w2);
		out.set(0,3, w1);
		out.set(1,0, -w3);
		out.set(1,2, w1);
		out.set(1,3, w2);
		out.set(2,0, w2);
		out.set(2,1, -w1);
		out.set(2,3, w3);
		out.set(3,0, -w1);
		out.set(3,1, -w2);
		out.set(3,2, -w3);
		return out;
	}    
        
}
