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

// A simple plot application with two plots 



package jat.demo.PlotExample;

// This class is not in the ptolemy.plot package so that it is a
// more realistic example.
import ptolemy.plot.*;

import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import javax.swing.JFrame;

// The java.io imports are only necessary for the right hand plot.
//import java.io.File;
//import java.io.FileInputStream;
//import java.io.FileNotFoundException;
//import java.io.IOException;

//////////////////////////////////////////////////////////////////////////
//// PlotExample
/**
PlotExample is a simple example that uses displays two plots side by side
To compile and run this application, do the following:
<pre>
javac -classpath ../../.. PlotExample.java
java -classpath ../../.. ptolemy.plot.demo.PlotExample
</pre>

@author Christopher Hylands
@version $Id: PlotExample.java,v 1.3 2004/03/03 14:55:16 dgaylor Exp $
*/
public class PlotExample extends JFrame {

    /** We use a constructor here so that we can call methods
     *  directly on the Frame.  The main method is static
     *  so getting at the Frame is a little trickier.
     */
    PlotExample() {
        // Instantiate the two plots.
        Plot topPlot = new Plot();
        Plot bottomPlot = new Plot();
        Plot middlePlot = new Plot();

        // Set the size of the toplevel window.
        setSize(560, 700);            // width, height in pixels

        // Create the top plot by calling methods.
	topPlot.setSize(350,300);
        topPlot.setButtons(true);
        topPlot.setTitle("Top Plot");
//        topPlot.setYRange(-4, 4);
//        topPlot.setXRange(0, 100);
        topPlot.setXLabel("time");
        topPlot.setYLabel("value");
        topPlot.addYTick("-PI", -Math.PI);
        topPlot.addYTick("-PI/2", -Math.PI/2);
        topPlot.addYTick("0",0);
        topPlot.addYTick("PI/2", Math.PI/2);
        topPlot.addYTick("PI", Math.PI);
        topPlot.setMarksStyle("none");
        topPlot.setImpulses(false);
        topPlot.addLegend(0,"dataset 1");
        topPlot.addLegend(1,"dataset 2");

        boolean first = true;
        for (int i = 0; i <= 100; i++) {
            topPlot.addPoint(0, (double)i,
                    5 * Math.cos(Math.PI * i/20), !first);
            topPlot.addPoint(1, (double)i,
                    4.5 * Math.cos(Math.PI * i/25), !first);
            topPlot.addPoint(2, (double)i,
                    4 * Math.cos(Math.PI * i/30), !first);
            topPlot.addPoint(3, (double)i,
                    3.5* Math.cos(Math.PI * i/35), !first);
            topPlot.addPoint(4, (double)i,
                    3 * Math.cos(Math.PI * i/40), !first);
            topPlot.addPoint(5, (double)i,
                    2.5 * Math.cos(Math.PI * i/45), !first);
            topPlot.addPoint(6, (double)i,
                    2 * Math.cos(Math.PI * i/50), !first);
            topPlot.addPoint(7, (double)i,
                    1.5 * Math.cos(Math.PI * i/55), !first);
            topPlot.addPoint(8, (double)i,
                    1 * Math.cos(Math.PI * i/60), !first);
            topPlot.addPoint(9, (double)i,
                    0.5 * Math.cos(Math.PI * i/65), !first);
            first = false;

        }

        // Create the bottom plot by calling methods.
	middlePlot.setSize(350,300);
        middlePlot.setButtons(true);
        middlePlot.setTitle("Middle Plot");
//        middlePlot.setYRange(-4, 4);
//        middlePlot.setXRange(0, 100);
        middlePlot.setXLabel("time");
        middlePlot.setYLabel("value");
        middlePlot.addYTick("-PI", -Math.PI);
        middlePlot.addYTick("-PI/2", -Math.PI/2);
        middlePlot.addYTick("0",0);
        middlePlot.addYTick("PI/2", Math.PI/2);
        middlePlot.addYTick("PI", Math.PI);
        middlePlot.setMarksStyle("none");
        middlePlot.setImpulses(false);

        for (int i = 0; i <= 100; i++) {
            middlePlot.addPoint(0, (double)i,
                    5 * Math.sin(Math.PI * i/20), !first);
            middlePlot.addPoint(1, (double)i,
                    4.5 * Math.sin(Math.PI * i/25), !first);
            middlePlot.addPoint(2, (double)i,
                    4 * Math.sin(Math.PI * i/30), !first);
            middlePlot.addPoint(3, (double)i,
                    3.5* Math.sin(Math.PI * i/35), !first);
            middlePlot.addPoint(4, (double)i,
                    3 * Math.sin(Math.PI * i/40), !first);
            middlePlot.addPoint(5, (double)i,
                    2.5 * Math.sin(Math.PI * i/45), !first);
            middlePlot.addPoint(6, (double)i,
                    2 * Math.sin(Math.PI * i/50), !first);
            middlePlot.addPoint(7, (double)i,
                    1.5 * Math.sin(Math.PI * i/55), !first);
            middlePlot.addPoint(8, (double)i,
                    1 * Math.sin(Math.PI * i/60), !first);
            middlePlot.addPoint(9, (double)i,
                    0.5 * Math.sin(Math.PI * i/65), !first);
            first = false;

        }

        // Create the bottom plot by calling methods.
	bottomPlot.setSize(350,300);
        bottomPlot.setButtons(true);
        bottomPlot.setTitle("Bottom Plot");
//        bottomPlot.setYRange(-4, 4);
//        bottomPlot.setXRange(0, 100);
        bottomPlot.setXLabel("time");
        bottomPlot.setYLabel("value");
        bottomPlot.addYTick("-PI", -Math.PI);
        bottomPlot.addYTick("-PI/2", -Math.PI/2);
        bottomPlot.addYTick("0",0);
        bottomPlot.addYTick("PI/2", Math.PI/2);
        bottomPlot.addYTick("PI", Math.PI);
        bottomPlot.setMarksStyle("none");
        bottomPlot.setImpulses(false);

        for (int i = 0; i <= 100; i++) {
            bottomPlot.addPoint(0, (double)i,
                    5 * Math.sin(Math.PI * i/20), !first);
            bottomPlot.addPoint(1, (double)i,
                    4.5 * Math.sin(Math.PI * i/25), !first);
            bottomPlot.addPoint(2, (double)i,
                    4 * Math.sin(Math.PI * i/30), !first);
            bottomPlot.addPoint(3, (double)i,
                    3.5* Math.sin(Math.PI * i/35), !first);
            bottomPlot.addPoint(4, (double)i,
                    3 * Math.sin(Math.PI * i/40), !first);
            bottomPlot.addPoint(5, (double)i,
                    2.5 * Math.sin(Math.PI * i/45), !first);
            bottomPlot.addPoint(6, (double)i,
                    2 * Math.sin(Math.PI * i/50), !first);
            bottomPlot.addPoint(7, (double)i,
                    1.5 * Math.sin(Math.PI * i/55), !first);
            bottomPlot.addPoint(8, (double)i,
                    1 * Math.sin(Math.PI * i/60), !first);
            bottomPlot.addPoint(9, (double)i,
                    0.5 * Math.sin(Math.PI * i/65), !first);
            first = false;

        }


        // Layout the two plots
        GridBagLayout gridbag = new GridBagLayout();
        GridBagConstraints c = new GridBagConstraints();
        getContentPane().setLayout(gridbag);

        // Handle the topPlot
        c.gridx = 0;
        c.gridy = 0;
        c.gridwidth = 1;
        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1.0;
        c.weighty = 1.0;
        gridbag.setConstraints(topPlot, c);
        getContentPane().add(topPlot);

        // Handle the middlePlot
        c.gridx = 0;
        c.gridy = 1;
        c.gridwidth = 1;
        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1.0;
        c.weighty = 1.0;
        gridbag.setConstraints(middlePlot, c);
        getContentPane().add(middlePlot);

        // Handle the bottomPlot
        c.gridx = 0;
        c.gridy = 2;
        c.gridwidth = 1;
        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1.0;
        c.weighty = 1.0;
        gridbag.setConstraints(bottomPlot, c);
        getContentPane().add(bottomPlot);

        setVisible(true);
    }


    /** main method called in a standalone java application.
     *  We simple instantiate this class, most of the work
     *  happens in the constructor.
     */
    public static void main(String args[]) {
        PlotExample PlotExample = new PlotExample();
    }
}
