/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2005 United States Government as represented by the
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
package jat.spacetime;

import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;

/**
 * Class used to perform translation or transformation between two
 * reference frames.  Whereas the ReferenceFrame class understands and 
 * simulates the complex dynamics of the reference frame, the translater
 * only works in a given time instant.
 * 
 * NOTE: This class is NOT units agnostic!  It assumes meters and meters/sec
 * in many of the translate*() methods.
 */
public class ReferenceFrameTranslater {
    
    /** The matrix to transform from the source reference frame
     * to the target reference frame. */
    private final Matrix xform;
    
    /** The origin of the target reference frame in terms of the
     * source reference frame (meters). */
    private final VectorN origin;
    
    /** The velocity of the target reference frame with respect to the
     * source reference frame in terms of the source reference frame
     * (meters/sec). */
    private final VectorN originVel;
    
    /** The rate of rotation of the target reference frame with respect
     * to the source reference frame in terms of the source reference frame. */
    private final VectorN omega;
    
    /** Whether this translater can be used to translate velocities. */
    private final boolean supportsVelocity;
    
    /**
     * Constructs a translater that translates from the source
     * reference frame to the target.
     * @param source the reference frame that you have coordinates in
     * @param target the reference frame you want to translate to
     * @param t the time at which translation will be done
     * @throws IllegalArgumentException if there is no translation
     * between the two reference frames
     */
    public ReferenceFrameTranslater(ReferenceFrame source, 
        ReferenceFrame target, Time t)
    {
      // We check with both the source and the target to see if
      // they know how to translate.
      ReferenceFrameTranslater xlater = source.getTranslater(target, t);
      if (xlater == null) {
        xlater = target.getTranslater(source, t);
        if (xlater == null) {
          throw new IllegalArgumentException("No known translation from " +
              source.getClass().getName() + " to " + 
              target.getClass().getName());
        }
        xlater = xlater.reverse();
      }
      
      // This is a constructor.  Copy the found values into this class.
      xform = xlater.xform;
      origin = xlater.origin;
      originVel = xlater.originVel;
      omega = xlater.omega;
      supportsVelocity = xlater.supportsVelocity;
    }
    
    /**
     * Construct a translater that does no translation.  Just spits
     * out what it is given.  This is useful when translating into
     * the same reference frame.
     */
    public ReferenceFrameTranslater()
    {
      xform = null;
      origin = null;
      originVel = null;
      omega = null;
      supportsVelocity = true;
    }

    /**
     * Construct a translater
     * @param transformationMatrix the transformation matrix for 
     * transforming directions from the source reference frame to
     * the target reference frame
     * @param originDifference the origin of the target reference frame
     * in terms of the source reference frame (meters)
     * @param originVelocity the velocity of the target reference frame with
     * respect to the source reference frame (meters/sec)
     */
    public ReferenceFrameTranslater(Matrix transformationMatrix,
        VectorN originDifference, VectorN originVelocity, VectorN rotation)
    {
      xform = transformationMatrix;
      origin = originDifference;
      originVel = originVelocity;
      omega = rotation;
      supportsVelocity = true;
    }
    
    /**
     * Construct a translater that does not support velocity transformations.
     * Useful when you know you will only do position transformations and
     * don't want to spend computation time figuring out
     * rotation and velocity dynamics.
     * @param transformationMatrix the transformation matrix for 
     * transforming directions from the source reference frame to
     * the target reference frame
     * @param originDifference the origin of the target reference frame
     * in terms of the source reference frame (meters)
     */
    public ReferenceFrameTranslater(Matrix transformationMatrix,
        VectorN originDifference)
    {
      xform = transformationMatrix;
      origin = originDifference;
      originVel = null;
      omega = null;
      supportsVelocity = false;
    }

    /**
     * Translates a given position's coordinates from the source 
     * reference frame to the coordinates in the target reference frame
     * @param coords x,y, and z coordinates in the source reference frame
     * (in meters)
     * @return x,y, and z coordinates in the target reference frame
     * (in meters)
     */
    public VectorN translatePoint(VectorN coords)
    {
      if ((origin == null) && (xform == null)) {
        coords = coords.copy();
      }
      else {
        coords = (origin == null ? coords : coords.minus(origin));
        coords = (xform == null ? coords : xform.times(coords));
      }
      return coords;
    }
    
    /**
     * Translates a given position's coordinates from the TARGET 
     * reference frame to the coordinates in the SOURCE reference frame
     * @param coords x,y, and z coordinates in the target reference frame
     * (in meters)
     * @return x,y, and z coordinates in the source reference frame
     * (in meters)
     */
    public VectorN translatePointBack(VectorN coords)
    {
      if ((origin == null) && (xform == null)) {
        coords = coords.copy();
      }
      else {
        coords = (xform == null ? coords : xform.transpose().times(coords));
        coords = (origin == null ? coords : coords.plus(origin));
      }
      return coords;
    }
    
    /**
     * Translates a given velocity's coordinates from the source 
     * reference frame to the coordinates in the target reference frame
     * @param velocity x,y, and z coordinates in the source reference frame
     * (in meters/sec)
     * @param position x,y, and z coordinates in the source reference frame
     * (in meters/sec)
     * @return x,y, and z coordinates in the target reference frame
     * (in meters/sec)
     */
    public VectorN translateVelocity(VectorN velocity, VectorN position)
    {
      if (!supportsVelocity) {
        throw new RuntimeException("Reference frames do not have a translation " +
                "that supports translating velocities.");
      }
      if ((originVel == null) && (xform == null) && (omega == null)) {
        velocity = velocity.copy();
      }
      else {
        velocity = (omega == null ? 
            velocity : velocity.minus(omega.crossProduct(position)));
        velocity = (xform == null ? velocity : xform.times(velocity));
        velocity = (originVel == null ? velocity : velocity.minus(originVel));
      }
      return velocity;
    }
    
    /**
     * Translates a given velocity's coordinates from the TARGET 
     * reference frame to the coordinates in the SOURCE reference frame
     * @param coords x,y, and z coordinates in the target reference frame
     * (in meters/sec)
     * @param position x,y, and z coordinates in the source reference frame
     * (in meters/sec)
     * @return x,y, and z coordinates in the source reference frame
     * (in meters/sec)
     */
    public VectorN translateVelocityBack(VectorN velocity, VectorN position)
    {
      if (!supportsVelocity) {
        throw new RuntimeException("Reference frames do not have a translation " +
                "that supports translating velocities.");
      }
      if ((originVel == null) && (xform == null) && (omega == null)) {
        velocity = velocity.copy();
      }
      else {
    	if (originVel != null) {
    		velocity = velocity.plus(originVel);
    	}
        if (xform != null) {
        	Matrix xformBack = xform.transpose();
            velocity = xformBack.times(velocity);
        	position = xformBack.times(position);       	
        }
        if (omega != null) {
        	velocity = velocity.plus(omega.crossProduct(position));
        }
      }
      return velocity;
    }
    
    /**
     * Transforms a given direction's coordinates in the source
     * reference frame to the coordinates in the target reference frame.
     * There is no translation to account for difference in
     * reference frames' origins  
     * @param coords x,y, and z coordinates in the source reference frame (meters)
     * @return x,y, and z coordinates in the target reference frame (meters)
     */
    public VectorN transformDirection(VectorN coords) 
    {
      return (xform == null ? coords.copy() : xform.times(coords));
    }
    
    /**
     * Transforms a given direction's coordinates from the TARGET
     * reference frame to the coordinates in the SOURCE reference frame.
     * There is no translation to account for difference in
     * reference frames' origins  
     * @param coords x,y, and z coordinates in the target reference frame (meters)
     * @return x,y, and z coordinates in the source reference frame (meters)
     */
    public VectorN transformDirectionBack(VectorN coords) 
    {
      return (xform == null ? coords.copy() : xform.transpose().times(coords));
    }
    
    /**
     * Transforms a given quaternian's angles in the source
     * reference frame to the angles in the target reference frame.
     * @param quat quaternian in the source reference frame
     * @return quaternian in the target reference frame
     */
    public double[] transformAngles(double[] quat)
    {
      // TODO: How do I transform quaternians?
      return null;
    }
    
    /**
     * Transforms a given quaternian's angles in the TARGET
     * reference frame to the angles in the SOURCE reference frame.
     * @param quat quaternian in the target reference frame
     * @return quaternian in the source reference frame
     */
    public double[] transformAnglesBack(double[] quat)
    {
      // TODO: How do I transform quaternians?
      return null;
    }   
    
    /**
     * Return a translater that translates from TARGET reference frame
     * to SOURCE reference frame, instread of source to target.
     * translater.translatePoint() == translater.reverse().translatePointBack()
     * @return a reverse translater
     */
    public ReferenceFrameTranslater reverse() {
      VectorN newOrigin = null;
      VectorN newVelocity = null;
      VectorN newRotation = null;
      Matrix newXform = (xform == null ? null : xform.transpose());
      if (origin != null) {
        // We have to reverse and transform the origin
        newOrigin = origin.times(-1);
        if (xform != null) {
          newOrigin = xform.times(newOrigin);
        }
      }
      if (originVel != null) {
        // We have to reverse and transform the origin velocity
        newVelocity = originVel.times(-1);
        if (xform != null) {
          newVelocity = xform.times(newVelocity);
        }
      }
      if (omega != null) {
        // We have to reverse and transform the rotation of the frame
        newRotation = omega.times(-1);
        if (xform != null) {
          newRotation = xform.times(newRotation);
        }
      }
      return new ReferenceFrameTranslater(newXform, newOrigin, newVelocity, newRotation);
    }
    
    public static void main(String[] args){
    	VectorN r = new VectorN(-7.0653447440413493e+002,  9.2387238910403778e+002, -1.4232051310790273e+003);
    	VectorN v = new VectorN(-7.6578459477044203e-001,  1.0064446005845396e+000,  1.0334979559135251e+000);
    	LunaFixedRef lfr = new LunaFixedRef();
		LunaRef lref = new LunaRef();
		Time t = new Time(58232);
		ReferenceFrameTranslater trans = new ReferenceFrameTranslater(lfr,lref,t);
		VectorN rLCI = trans.transformDirection(r);
		VectorN vLCI = trans.translateVelocity(v,r);
		System.out.println("rLCI: "+rLCI);
		System.out.println("vLCI: "+vLCI);
    }
}
