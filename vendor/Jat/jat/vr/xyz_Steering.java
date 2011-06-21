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
*  File        :   xyz_Steering.java
*  Author      :   Tobias Berthold
*  Date        :   9-27-2002
*  Change      :
*  Description :   
*/

package jat.vr;

import java.awt.event.*;
import java.awt.AWTEvent;
import java.util.Enumeration;
import javax.vecmath.*;
import javax.media.j3d.*;

/**
 * This class is a simple behavior that invokes the KeyNavigator
 * to modify the view platform transform.
 *
 * @author Tobias Berthold
 * @version 1.0
 */
public class xyz_Steering extends Behavior
{
	private WakeupOnAWTEvent wakeupOne = null;
	private WakeupCriterion[] wakeupArray = new WakeupCriterion[1];
	private WakeupCondition wakeupCondition = null;
	private float TRANSLATE = 1.e7f;

	TransformGroup m_TransformGroup = null;

	public xyz_Steering(TransformGroup tg)
	{
		m_TransformGroup = tg;

		try
		{
			m_TransformGroup.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
			m_TransformGroup.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
		} catch (Exception e)
		{
		}

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

	public void set_translate(float translate)
	{
		this.TRANSLATE = translate;
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
			genericEvt = (WakeupCriterion)criteria.nextElement();

			if (genericEvt instanceof WakeupOnAWTEvent)
			{
				ev = (WakeupOnAWTEvent)genericEvt;
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
				KeyEvent eventKey = (KeyEvent)events[n];

				if (eventKey.getID() == KeyEvent.KEY_PRESSED)
				{
					int keyCode = eventKey.getKeyCode();
					int keyChar = eventKey.getKeyChar();

					Vector3f translate = new Vector3f();

					Transform3D t3d = new Transform3D();
					m_TransformGroup.getTransform(t3d);
					t3d.get(translate);

					switch (keyCode)
					{
						case KeyEvent.VK_PAGE_UP :
							translate.z += TRANSLATE;
							break;
						case KeyEvent.VK_PAGE_DOWN :
							translate.z -= TRANSLATE;
							break;
						case KeyEvent.VK_UP :
							translate.y += TRANSLATE;
							break;
						case KeyEvent.VK_DOWN :
							translate.y -= TRANSLATE;
							break;
						case KeyEvent.VK_LEFT :
							translate.x += TRANSLATE;
							break;
						case KeyEvent.VK_RIGHT :
							translate.x -= TRANSLATE;
							break;
					}

					// System.out.println( "Steering: " + translate.x );
					//translate.y = 0.5f;

					t3d.setTranslation(translate);
					m_TransformGroup.setTransform(t3d);
				}
			}
		}
	}
}
