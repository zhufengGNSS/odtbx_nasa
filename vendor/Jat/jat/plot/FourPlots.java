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

package jat.plot;

import ptolemy.plot.*;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import javax.swing.JFrame;
/**
 * <P>
 * The ThreePlots Class provides a way to create a page with three plots.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class FourPlots extends JFrame {
	
	/** The uppermost plot */
    public Plot firstPlot = new Plot();
        
    /** The plot that is second from the top */
    public Plot secondPlot = new Plot();
    
    /** The plot that is third from the top */
    public Plot thirdPlot = new Plot();
    
    /** The bottom plot */
    public Plot fourthPlot = new Plot();


    /** Default constructor.
     */
    public FourPlots() {

        // Set the size of the toplevel window.
        setSize(560, 700);            // width, height in pixels

        // Create the top plot by calling methods.
        firstPlot.setSize(350,300);
        firstPlot.setButtons(true);

        // Create the second plot by calling methods.
        secondPlot.setSize(350,300);
        secondPlot.setButtons(true);

        // Create the third plot by calling methods.
        thirdPlot.setSize(350,300);
        thirdPlot.setButtons(true);

        // Create the bottom plot by calling methods.
        fourthPlot.setSize(350,300);
        fourthPlot.setButtons(true);


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
        gridbag.setConstraints(firstPlot, c);
        getContentPane().add(firstPlot);

        // Handle the secondPlot
        c.gridx = 0;
        c.gridy = 1;
        c.gridwidth = 1;
        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1.0;
        c.weighty = 1.0;
        gridbag.setConstraints(secondPlot, c);
        getContentPane().add(secondPlot);

        // Handle the thirdPlot
        c.gridx = 0;
        c.gridy = 2;
        c.gridwidth = 1;
        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1.0;
        c.weighty = 1.0;
        gridbag.setConstraints(thirdPlot, c);
        getContentPane().add(thirdPlot);

        // Handle the bottomPlot
        c.gridx = 0;
        c.gridy = 3;
        c.gridwidth = 1;
        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1.0;
        c.weighty = 1.0;
        gridbag.setConstraints(fourthPlot, c);
        getContentPane().add(fourthPlot);

    }

}
