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

package jat.vr;

import java.awt.event.*;
import java.awt.AWTEvent;
import java.util.Enumeration;
import javax.vecmath.*;
import javax.media.j3d.*;

/**
* This class is a simple behavior that invokes the KeyNavigator
* to modify the view platform transform.
*/
public class spherical_Steering extends Behavior
{
	private WakeupOnAWTEvent wakeupOne = null;
	private WakeupCriterion[] wakeupArray = new WakeupCriterion[1];
	private WakeupCondition wakeupCondition = null;
	private final float TRANSLATE_FROM_TO = 500.f;
	private final float	TRANSLATE_UP_DOWN = 0.01f;
	private final float	TRANSLATE_LEFT_RIGHT = 0.01f;
	static double phi,theta,r;
	TransformGroup TG_vp = null;
    Point3d origin=new Point3d(0.0f,0.0f,0.0f);

	public spherical_Steering( TransformGroup tg )
	{
        TG_vp = tg;
        theta=0.;
        phi=0.;

        /*
        try
        {
            m_TransformGroup.setCapability( TransformGroup.ALLOW_TRANSFORM_WRITE );
            m_TransformGroup.setCapability( TransformGroup.ALLOW_TRANSFORM_READ );
        }
        catch( Exception e )
        {
        }
        */
        wakeupOne = new WakeupOnAWTEvent(KeyEvent.KEY_PRESSED);
        wakeupArray[0] = wakeupOne;
        wakeupCondition = new WakeupOr( wakeupArray );
	}


	/**
	*  Override Behavior's initialize method to setup wakeup criteria.
	*/
	public void initialize( )
	{
		// Establish initial wakeup criteria
		wakeupOn( wakeupCondition );
	}

	/**
	*  Override Behavior's stimulus method to handle the event.
	*/
	public void processStimulus( Enumeration criteria )
	{
		WakeupOnAWTEvent ev;
		WakeupCriterion genericEvt;
		AWTEvent[] events;

		while (criteria.hasMoreElements( ))
		{
			genericEvt = (WakeupCriterion) criteria.nextElement( );

			if (genericEvt instanceof WakeupOnAWTEvent)
			{
				ev = (WakeupOnAWTEvent) genericEvt;
				events = ev.getAWTEvent( );
				processAWTEvent( events );
			}
		}

		// Set wakeup criteria for next time
		wakeupOn( wakeupCondition );
	}

	/**
	*  Process a keyboard event
	*/
	private void processAWTEvent( AWTEvent[] events )
	{
		for( int n = 0; n < events.length; n++)
		{
			if( events[n] instanceof KeyEvent)
			{
				KeyEvent eventKey = (KeyEvent) events[n];

				if( eventKey.getID( ) == KeyEvent.KEY_PRESSED )
				{
					int keyCode = eventKey.getKeyCode( );
					int keyChar = eventKey.getKeyChar( );

					//Vector3f translate = new Vector3f( );
					Vector3d translate = new Vector3d( );

					Transform3D t3d = new Transform3D( );
					TG_vp.getTransform( t3d );
					t3d.get( translate );
    			    r=Math.sqrt(translate.x*translate.x+translate.y*translate.y+translate.z*translate.z);


					switch (keyCode)
					{
					case KeyEvent.VK_ADD:
					    r+=TRANSLATE_FROM_TO;
						break;
					case KeyEvent.VK_SUBTRACT:
					    r-=TRANSLATE_FROM_TO;
						break;
					case KeyEvent.VK_UP:
					    phi+=TRANSLATE_UP_DOWN;
						break;
					case KeyEvent.VK_DOWN:
					    phi-=TRANSLATE_UP_DOWN;
						break;
					case KeyEvent.VK_LEFT:
					    theta+=TRANSLATE_LEFT_RIGHT;
					    //System.out.println(""+theta);
						break;
					case KeyEvent.VK_RIGHT:
					    theta-=TRANSLATE_LEFT_RIGHT;
						break;
					}

					// System.out.println( "Steering: " + translate.x );
					//translate.y = 0.5f;

					translate.x=(float)(r*Math.cos(theta)*Math.sin(phi));
					translate.y=(float)(r*Math.sin(theta)*Math.sin(phi));
					translate.z=(float)(r*Math.cos(phi));
                    //t3d=jat_view.lookat_T(translate, origin);

					t3d.setTranslation( translate );
					TG_vp.setTransform( t3d );
				}
			}
		}
	}
}
