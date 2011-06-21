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
 /*
 *      EricTrajectory.java
 *
 *
 */
package jat.vr;

import jat.util.*;
import javax.media.j3d.*;
import javax.vecmath.*;


public class EricTrajectory extends Shape3D
{
    Color3f gray=new Color3f(0.3f,0.3f,0.3f);
    Color3f pink=new Color3f(0.9f,0.1f,0.1f);
    public double[] coords;
    public double[] time;

    public EricTrajectory(double[] coords)
    {
        this.coords=coords;
    }

    public EricTrajectory(String filename)
    {
        EasyReader inFile = new EasyReader(filename);
        if (inFile.bad())
        {
            System.err.println("Can't open " + filename);
            System.exit(1);
        }

        // find out how many numbers in file
        int num=0;
        while (!inFile.eof())
        {
            inFile.readWord();
            num++;
        }
        inFile.close();
        num--;

        // read the data from file
        inFile = new EasyReader(filename);
        time=new double[num/4];
        coords=new double[(num)*3/4];

        int i,line=0;
        for(i=0;i<num/4;i++)
        {
            time[i]        = inFile.readDouble();
            coords[line+0] = inFile.readDouble();
            coords[line+1] = inFile.readDouble();
            coords[line+2] = inFile.readDouble();
            //System.out.println(""+coords[line+2]);
            inFile.readLine();
            line+=3;
            //System.out.println(""+x);
        }

        //double x = inFile.readDouble();  // reads a double
        //System.out.println(""+x);
        inFile.close();

        draw_trajectory();

        //draw_lines_simple();

    }

    public EricTrajectory()
    {
        draw_trajectory();
    }

    private void draw_trajectory()
    {

        int num_vert=coords.length/3;
        //int[] stripLengths = { 200};
        int[] stripLengths = { num_vert};

        LineStripArray myLines = new LineStripArray
        (
            num_vert,
            GeometryArray.COORDINATES,
            stripLengths
        );
        myLines.setCoordinates( 0, coords );

        this.setGeometry(myLines);
    }

    private void draw_lines_simple()
    {
        double[] coords=
        {
            0.,0.,0.,
            100000.,0.,0.,
            1000.,1000.,0.,
            1000.,1000.,1000.
        };

        int[] stripLengths = { 4};

        LineStripArray myLines = new LineStripArray
        (
            4,
            GeometryArray.COORDINATES,
            stripLengths
        );
        myLines.setCoordinates( 0, coords );

        this.setGeometry(myLines);
    }

}
