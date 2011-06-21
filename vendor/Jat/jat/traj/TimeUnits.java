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
* The TimeUnits.java Class provides the means for specifying the 
* time units used in creating a trajectory.
*
* @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
* @version 1.0
*/

public final class TimeUnits implements Serializable {

  private String name;

  private TimeUnits(String nm) { name = nm; }

  public String toString() { return name; }

  public final static TimeUnits
    SECONDS = new TimeUnits("s"),
    DAYS = new TimeUnits("days"),
    OTHER = new TimeUnits("Other");


  public final static TimeUnits[] index =  {
    SECONDS, DAYS, OTHER
  };


  public static void main(String[] args) {
    TimeUnits m = TimeUnits.SECONDS;
    System.out.println(m);
    m = TimeUnits.index[1];
    System.out.println(m);
    System.out.println(m == TimeUnits.SECONDS);
    System.out.println(m.equals(TimeUnits.DAYS));
  }
}