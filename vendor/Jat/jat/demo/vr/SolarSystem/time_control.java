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
*  File        :   time_control.java
*  Author      :   Tobias Berthold
*  Date        :   11-13-2002
*  Change      :
*  Description :
*/

package jat.demo.vr.SolarSystem;

import java.awt.event.*;
import java.awt.AWTEvent;
import java.util.Enumeration;
import javax.media.j3d.*;

/**
* This class is a simple behavior that invokes the KeyNavigator
* to modify the time.
*/
public class time_control extends Behavior
{
	private WakeupOnAWTEvent wakeupOne = null;
	private WakeupCriterion[] wakeupArray = new WakeupCriterion[1];
	private WakeupCondition wakeupCondition = null;
	private final float TRANSLATE = 1.e6f;
	Constellation sr;

	public time_control(Constellation sr)
	{
		this.sr = sr;
		wakeupOne = new WakeupOnAWTEvent(KeyEvent.KEY_PRESSED);
		wakeupArray[0] = wakeupOne;
		wakeupCondition = new WakeupOr(wakeupArray);
	}

	/**
	*  Override Behavior's initialize method to setup wakeup criteria.
	*/
	public void initialize()
	{
		// Establish initial wakeup criteria
		wakeupOn(wakeupCondition);
	}

	/**
	*  Override Behavior's stimulus method to handle the event.
	*/
	public void processStimulus(Enumeration criteria)
	{
		WakeupOnAWTEvent ev;
		WakeupCriterion genericEvt;
		AWTEvent[] events;

		while (criteria.hasMoreElements())
		{
			genericEvt = (WakeupCriterion) criteria.nextElement();

			if (genericEvt instanceof WakeupOnAWTEvent)
			{
				ev = (WakeupOnAWTEvent) genericEvt;
				events = ev.getAWTEvent();
				processAWTEvent(events);
			}
		}

		// Set wakeup criteria for next time
		wakeupOn(wakeupCondition);
	}

	/**
	*  Process a keyboard event
	*/
	private void processAWTEvent(AWTEvent[] events)
	{
		for (int n = 0; n < events.length; n++)
		{
			if (events[n] instanceof KeyEvent)
			{
				KeyEvent eventKey = (KeyEvent) events[n];

				if (eventKey.getID() == KeyEvent.KEY_PRESSED)
				{

					int keyCode = eventKey.getKeyCode();
					int keyChar = eventKey.getKeyChar();

					switch (keyCode)
					{
						case KeyEvent.VK_ADD :
							//sr.SimClock.jd+=1.;
							//sr.SimClock.dd.setTime(sr.SimClock.dd.getTime() + 86400000);
							//sr.dd.setTime(sr.dd.getTime()+86400000);
							//sr.SimClock.dd+=86400.;
							break;
						case KeyEvent.VK_SUBTRACT :
							//sr.SimClock.dd.setTime(sr.SimClock.dd.getTime() - 86400000);
							//sr.dd.setTime(sr.dd.getTime()-86400000);
							break;
					}

					//System.out.println( "Steering: "  );
				}
			}
		}
	}
}
