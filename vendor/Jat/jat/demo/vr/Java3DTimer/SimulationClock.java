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

package jat.demo.vr.Java3DTimer;

import java.util.Enumeration;
import javax.media.j3d.*;
import javax.vecmath.*;

/**
 * @author Tobias Berthold
 *
 */
public class SimulationClock extends Behavior
{
	WakeupCriterion yawn;
	long SimTime;
	int DeltaT;
	Vector3d V = new Vector3d(0.f, 0.f, 0.0f);
	Java3DTimer JT;
	ControlPanel panel;
	int i;

	public SimulationClock(int ts, Java3DTimer JT, ControlPanel panel)
	{
		DeltaT = ts;
		this.JT = JT;
		this.panel = panel;
		yawn = new WakeupOnElapsedTime(DeltaT);
	}

	public void initialize()
	{
		wakeupOn(yawn);
		SimTime = 0;
	}

	public void processStimulus(Enumeration e)
	{
		// general animation
		SimTime += DeltaT;
		wakeupOn(yawn);

		// Update text in panel
		panel.label.setText("Time " + SimTime );

		// Rotate the cube
		i++;
		JT.cube.set_attitude(Math.PI / 2., i * 0.001f, 0);

	} //  end method

} //  end class
