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

package jat.cm;
import jat.plot.*;
import jat.util.*;
import jat.traj.*;
import jat.matvec.data.*;
import java.io.*;
import java.util.*;

/**
 * <P>
 * The GroundTrack Class provides the ability to plot ground tracks.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */

public class GroundTrack {

	SinglePlot gt = new SinglePlot();
	private int gt_index = 0;

	/** Creates a new instance of GroundTrack */
	public GroundTrack() {
		// set up the trajectory plot
		gt.setTitle("Ground Track");
		gt.plot.setXLabel("Long (deg)");
		gt.plot.setYLabel("Lat (deg)");

	}

	/** Read the trajectory data from a tab-delimited ASCII text file.
	 * The first line read in sets the number of columns.
	 * @param file filename and directory
	 */
	public void addCoasts() {
		String path=FileUtil.getClassFilePath("jat.cm","GroundTrack");
		String file = path + "new_coast_50.dat";

		try {
			FileReader fr = new FileReader(file);
			BufferedReader in = new BufferedReader(fr);
			String line;
			boolean first = false;

			int coast = 0;

			// loop through the file, one line at a time
			while ((line = in.readLine()) != null) {
				StringTokenizer tok = new StringTokenizer(line, "\t");
				String token = tok.nextToken();
				if (token.equals("#")) {
					coast = coast + 1;
					first = false;
				} else {
					double lon = Double.parseDouble(token);
					token = tok.nextToken();
					double lat = Double.parseDouble(token);
					this.gt.plot.addPoint(0, lon, lat, first);
					first = true;
				}
			}
			in.close();
			fr.close();
			this.gt_index++;
		} catch (IOException e) {
			System.err.println("Error opening:" + file);
			return;
		}
	}
	
	/**
	 * Add a trajectory to the ground track
	 * @param traj trajectory to be added
	 */
	public void addTrajectory(Trajectory traj){
		boolean first = false;
		double oldlon = 0.0;
		
		// loop through the trajectory
		while (traj.hasNext()) {

			// grab the next row and strip out t and inertial r vector
			double [] row = traj.next();
			double t = row[0];
			VectorN r_ijk = new VectorN(row[1], row[2], row[3]);
			
			// convert to r to ECEF
			double alpha = Constants.omega_e * t;
			RotationMatrix R = new RotationMatrix(3, alpha);
			VectorN r_ecef = R.times(r_ijk);
			
			// get lat, lon, alt
//			Geodetic geo = new Geodetic(r_ecef);
//			double lat = geo.getLatitude() * c.rad2deg;
//			double lon = geo.getLongitude() * c.rad2deg;

			double lat = Math.asin(r_ecef.x[2]/r_ecef.mag())* Constants.rad2deg;
			double lon = Math.atan2(r_ecef.x[1], r_ecef.x[0])* Constants.rad2deg;
			
			// ensure -180 < lon < 180			
			if (lon > 180.0) lon = lon - 360.0;
			if (lon < -180.0) lon = lon + 360.0;
			
			// check for boundary crossing and disconnect the curve
			if (Math.abs(lon - oldlon) > 180.0) {
				first = false;
			}
			
//			System.out.println(t+" "+lon+" "+lat+" "+first);
			
			// add to plot		
			this.gt.plot.addPoint(this.gt_index, lon, lat, first);
			
			oldlon = lon;
			first = true;
		}
		this.gt_index++;					
	}

	/** Runs the example.
	 * @param args Arguments.
	 */
	public static void main(String[] args) {
		GroundTrack gtrack = new GroundTrack();
		gtrack.addCoasts();
		gtrack.gt.plot.setVisible(true);
		
		// create a Trajectory
		Trajectory traj1 = new Trajectory();
		
        // create a TwoBody orbit using orbit elements
        TwoBody sat = new TwoBody(20000.0, 0.3, 45.0, 0.0, 0.0, 0.0);

        // find out the period of the orbit
        double period = sat.period();

        // set the final time = one orbit period
        double tf = 2.0 * period;

        // set the initial time to zero
        double t0 = 0.0;

        // propagate the orbit
        sat.propagate(t0, tf, traj1, true);
        
        // set the attributes
        traj1.setTitle("Orbit 1");
        traj1.setCentralBody(CentralBody.EARTH);
        traj1.setCoordinateSystem(CoordinateSystem.INERTIAL);
        traj1.setDistanceUnits(DistanceUnits.KILOMETERS);
        traj1.setEpoch(2003, 3, 27, 12, 16, 0.0);
        traj1.setTimeUnits(TimeUnits.SECONDS);
        String[] labels = {"t","x","y","z","xdot","ydot","zdot"};
        traj1.setLabels(labels);

		// create a Trajectory
		Trajectory traj2 = new Trajectory();
		
        // create a TwoBody orbit using orbit elements
        TwoBody sat2 = new TwoBody(26559.0, 0.75, 63.4, 30.0, 270.0, 0.0);

        // find out the period of the orbit
        period = sat2.period();

        // set the final time = one orbit period
        tf = 2.0 * period;

        // set the initial time to zero
        t0 = 0.0;

        // propagate the orbit
        sat2.propagate(t0, tf, traj2, true);
        
        // set the attributes
        traj2.setTitle("Molniya");
        traj2.setCentralBody(CentralBody.EARTH);
        traj2.setCoordinateSystem(CoordinateSystem.INERTIAL);
        traj2.setDistanceUnits(DistanceUnits.KILOMETERS);
        traj2.setEpoch(2003, 3, 27, 12, 16, 0.0);
        traj2.setTimeUnits(TimeUnits.SECONDS);
        traj2.setLabels(labels);

		// create a Trajectory
		Trajectory traj3 = new Trajectory();
		
        // create a TwoBody orbit using orbit elements
        TwoBody sat3 = new TwoBody(42166.26, 0.2, 45.0, 0.0, 0.0, 0.0);

        // find out the period of the orbit
        period = sat3.period();

        // set the final time = one orbit period
        tf = 2.0 * period;

        // set the initial time to zero
        t0 = 0.0;

        // propagate the orbit
        sat3.propagate(t0, tf, traj3, true);
        
        // set the attributes
        traj3.setTitle("Inclined Geosynch");
        traj3.setCentralBody(CentralBody.EARTH);
        traj3.setCoordinateSystem(CoordinateSystem.INERTIAL);
        traj3.setDistanceUnits(DistanceUnits.KILOMETERS);
        traj3.setEpoch(2003, 3, 27, 12, 16, 0.0);
        traj3.setTimeUnits(TimeUnits.SECONDS);
        traj3.setLabels(labels);

        
        // add trajectory to ground track
        gtrack.addTrajectory(traj1);
        gtrack.addTrajectory(traj2);
        gtrack.addTrajectory(traj3);
        
        gtrack.gt.plot.addLegend(1, traj1.getTitle());
        gtrack.gt.plot.addLegend(2, traj2.getTitle());
        gtrack.gt.plot.addLegend(3, traj3.getTitle());
        
        gtrack.gt.setVisible(true);	

	}
}
