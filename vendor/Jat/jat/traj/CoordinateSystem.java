package jat.traj;

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
 import java.io.*;
 
/**
* <P>
* The CoordinateSystem.java Class provides the means for specifying the 
* coordinate system used in creating a trajectory.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/ 

public final class CoordinateSystem implements Serializable {

  private String name;

  private CoordinateSystem(String nm) { name = nm; }

  public String toString() { return name; }

  public final static CoordinateSystem
    INERTIAL = new CoordinateSystem("Inertial"),
    PLANETFIXED = new CoordinateSystem("PlanetFixed"),
    LLH = new CoordinateSystem("LLH"),
    OTHER = new CoordinateSystem("Other");


  public final static CoordinateSystem[] index =  {
    INERTIAL, PLANETFIXED, LLH, OTHER
  };


  public static void main(String[] args) {
    CoordinateSystem m = CoordinateSystem.INERTIAL;
    System.out.println(m);
    m = CoordinateSystem.index[1];
    System.out.println(m);
    System.out.println(m == CoordinateSystem.PLANETFIXED);
    System.out.println(m.equals(CoordinateSystem.INERTIAL));
  }
}