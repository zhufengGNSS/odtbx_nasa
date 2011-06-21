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
 *  File        :   jat_behavior.java
 *  Author      :   Tobias Berthold
 *  Date        :   10-7-2002
 *  Change      :
 *  Description :   behavior related functions
 */
package jat.vr;

import com.sun.j3d.utils.behaviors.mouse.*;
import javax.media.j3d.*;

// TODO: Currently have to click into window first, before keyboard inputs are recognized 
public class jat_behavior
{
	public static xyz_Steering xyz_Behavior;
	private static float TRANSLATE = 1.e2f;

//	public  xyz_Steering xyz_Behavior;
//	private float TRANSLATE = 1.e2f;
   
	public static void behavior(BranchGroup BG_root, BranchGroup BG_vp, BoundingSphere bounds)
//	public void behavior(BranchGroup BG_root, BranchGroup BG_vp, BoundingSphere bounds)
    {

        xyz_Behavior = new xyz_Steering(jat_view.TG_vp);
        xyz_Behavior.setSchedulingBounds(bounds);
        BG_vp.addChild(xyz_Behavior);

        /*
        spherical_Steering sph_Behavior = new spherical_Steering(jat_view.TG_vp);
        sph_Behavior.setSchedulingBounds(bounds);
        BG_vp.addChild(sph_Behavior);
        */

        MouseRotate RotationBehavior = new MouseRotate(MouseBehavior.INVERT_INPUT);
        RotationBehavior.setSchedulingBounds(bounds);
        RotationBehavior.setFactor(0.01);
        RotationBehavior.setTransformGroup(jat_view.TG_vp);
        BG_vp.addChild(RotationBehavior);

        MouseZoom ZoomBehavior = new MouseZoom(MouseBehavior.INVERT_INPUT);
        ZoomBehavior.setSchedulingBounds(bounds);
        ZoomBehavior.setFactor(TRANSLATE);
        ZoomBehavior.setTransformGroup(jat_view.TG_vp);
        BG_vp.addChild(ZoomBehavior);

        /*
        SimulationClock SimClock = new SimulationClock(40, em.bod, em.panel, em.tra1);
        SimClock.setSchedulingBounds(bounds);
        BG_vp.addChild(SimClock);
        */
    }
    
	public void set_translate(float translate)
	{
//		this.TRANSLATE = translate;
	}

}
