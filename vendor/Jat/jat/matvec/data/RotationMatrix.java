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
 * The RotationMatrix Class provides the fundamental operations of direction
 * cosine matrices.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class RotationMatrix extends Matrix {

	private static final long serialVersionUID = 1L;
    
/* ------------------------
   Constructors
 * ------------------------ */
    
    /** Construct a 3x3 identity matrix.
     */
    
    public RotationMatrix() {
        super(3);
    }
    
    /** Construct a copy of a Rotation Matrix.
     * @param in   Rotation matrix to make a copy of.
     */
    
    public RotationMatrix(RotationMatrix in) {
        super(in);
        in.checkMatrixDimensions(3,3);
    }
    
    /** Construct a copy of a Rotation Matrix.
     * @param in   Matrix to make a copy of.
     */
    
    public RotationMatrix(Matrix in) {
        super(in);
        in.checkMatrixDimensions(3,3);
    }

    /** Construct a small angle Rotation matrix from 3 small angles.
     * @param in   Matrix to make a copy of.
     */
    
    public RotationMatrix(double a, double b, double c) {
        super(3);
        this.A[0][1] = c;
        this.A[0][2] = -b;
        this.A[1][0] = -c;
        this.A[1][2] = a;
        this.A[2][0] = b;
        this.A[2][1] = -a;
    }    
    
    /** Construct a single axis rotation matrix.
     * @param n    Axis to rotate about.
     * @param phi  Angle to rotate through in radians.
     */
    
    public RotationMatrix(int n, double phi) {
        super(3);
        double cos_phi = Math.cos(phi);
        double sin_phi = Math.sin(phi);
        double msin_phi = -1.0 * sin_phi;
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                this.A[i][j] = 0.0;
            }
            this.A[i][i] = cos_phi;
        }
        switch (n) {
            case 1:      // Rotate about 1st axis
                this.A[0][0] = 1.0;
                this.A[1][2] = sin_phi;
                this.A[2][1] = msin_phi;
                break;
            case 2:      // Rotate about 2nd axis
                this.A[1][1] = 1.0;
                this.A[0][2] = msin_phi;
                this.A[2][0] = sin_phi;
                break;
            case 3:      // Rotate about 3rd axis
                this.A[2][2] = 1.0;
                this.A[0][1] = sin_phi;
                this.A[1][0] = msin_phi;
                break;
            default:
                System.out.println("Only rotations about 1, 2 or 3 are allowed");
                break;
        }
    }
    
    /** Construct a two axis rotation matrix.
     * @param n      1st axis to rotate about.
     * @param theta  Angle to rotate through in radians.
     * @param m      2nd axis to rotate about.
     * @param phi    Angle to rotate through in radians.
     */
    public RotationMatrix(int n, double theta, int m, double phi) {
        super(3);
        RotationMatrix t1 = new RotationMatrix(n, theta);
        RotationMatrix t2 = new RotationMatrix(m, phi);
        Matrix temp = t2.times(t1);
        double [][] out = temp.getArrayCopy();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                this.A[i][j] = out[i][j];
            }
        }
    }
    
    /** Construct a three axis rotation matrix.
     * @param n      1st axis to rotate about.
     * @param psi    Angle to rotate through in radians.
     * @param m      2nd axis to rotate about.
     * @param theta  Angle to rotate through in radians.
     * @param l      3rd axis to rotate about.
     * @param phi    Angle to rotate through in radians.
     */
    
    public RotationMatrix(int n, double psi, int m, double theta, int l, double phi) {
        super(3);
        RotationMatrix t1 = new RotationMatrix(n, psi);
        RotationMatrix t2 = new RotationMatrix(m, theta);
        RotationMatrix t3 = new RotationMatrix(l, phi);
        Matrix temp1 = t2.times(t1);
        Matrix temp = t3.times(temp1);
        double [][] out = temp.getArrayCopy();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                this.A[i][j] = out[i][j];
            }
        }
    }
/* ------------------------
   public methods
 * ------------------------ */
    
    /** Perform the transformation on the vector.
     * @param in   Vector3 to rotate.
     * @return     Vector3 expressed in new frame.
     */
    
    public VectorN transform(VectorN in) {
        in.checkVectorDimensions(3);
        VectorN out = this.times(in);
        return out;
    }
    
    /** Find the inverse rotation matrix.
     * @return     Inverse of this rotation matrix.
     */
    
    public Matrix inverse() {
        Matrix out = this.transpose();
        return out;
    }
}
