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
 */

package jat.application.EarthSatellite;

import org.eclipse.swt.SWT;
import org.eclipse.swt.opengl.GLCanvas;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.Event;
import org.eclipse.swt.widgets.Listener;
import org.lwjgl.LWJGLException;
import org.lwjgl.opengl.GL11;
import org.lwjgl.opengl.GLContext;
import org.eclipse.swt.graphics.*;
import org.eclipse.swt.opengl.GLData;
import org.lwjgl.opengl.glu.GLU;

public class AnimationPane extends Composite implements Runnable
{
	Display display;
	final GLCanvas canvas;
	int rot;
	float frot = 1.0f;
	float faWhite[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	float faLightBlue[] = { 0.8f, 0.8f, .9f, 1f };
	float lightPosition[] = { -4f, 3f, 3f, 1.0f };
	int sphereTextureHandle = 0;

	public AnimationPane(Display display, Composite parent)
	{
		super(parent, SWT.NONE);
		this.display = display;
		GLData data = new GLData();
		data.doubleBuffer = true;
		canvas = new GLCanvas(this, SWT.NONE, data);
		canvas.setCurrent();
		try
		{
			GLContext.useContext(canvas);
		} catch (LWJGLException e)
		{
			e.printStackTrace();
		}
		canvas.addListener(SWT.Resize, new Listener()
		{
			public void handleEvent(Event event)
			{
				Rectangle bounds = canvas.getBounds();
				float fAspect = (float) bounds.width / (float) bounds.height;
				canvas.setCurrent();
				try
				{
					GLContext.useContext(canvas);
				} catch (LWJGLException e)
				{
					e.printStackTrace();
				}
				GL11.glViewport(0, 0, bounds.width, bounds.height);
				GL11.glMatrixMode(GL11.GL_PROJECTION);
				GL11.glLoadIdentity();
				GLU.gluPerspective(45.0f, fAspect, 0.5f, 400.0f);
				GL11.glMatrixMode(GL11.GL_MODELVIEW);
				GL11.glLoadIdentity();
			}
		});
		
		//GL11.glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
		//GL11.glColor3f(0.1f, 0.0f, 0.4f);
		GL11.glHint(GL11.GL_PERSPECTIVE_CORRECTION_HINT, GL11.GL_NICEST);
		GL11.glClearDepth(1.0);
		GL11.glLineWidth(1);
		GL11.glEnable(GL11.GL_DEPTH_TEST);
		// enable lighting and texture rendering
		GL11.glEnable(GL11.GL_LIGHTING);
		GL11.glEnable(GL11.GL_TEXTURE_2D);
		// Create sphere texture
		GLImage textureImg = GLlib.loadImage("images/earth.jpg");
		sphereTextureHandle = GLlib.makeTexture(textureImg);
		display.timerExec(50, this);
	}

	public void run()
	{
		{
			if (!canvas.isDisposed())
			{
				canvas.setCurrent();
				try
				{
					GLContext.useContext(canvas);
				} catch (LWJGLException e)
				{
					e.printStackTrace();
				}
				GL11.glClear(GL11.GL_COLOR_BUFFER_BIT | GL11.GL_DEPTH_BUFFER_BIT);
				//GL11.glClearColor(.1f, .1f, .1f, 1.0f);
				GL11.glLoadIdentity();
				GL11.glTranslatef(0.0f, 0.0f, -5.0f);
				//GL11.glRotatef(0.15f * rot, 2.0f * frot, 10.0f * frot, 1.0f);
				GL11.glRotatef(frot, 0.0f , 1.0f * frot, 1.0f);
				frot=frot+0.1f;
				GL11.glPolygonMode(GL11.GL_FRONT_AND_BACK, GL11.GL_SMOOTH);
				//GL11.glColor3f(0.9f, 0.5f, 0.3f);
				// Create a light (diffuse light, ambient light, position)
				GLlib.setLight(GL11.GL_LIGHT1, faWhite, faLightBlue, lightPosition);
				// drawTorus(1, 1.9f + ((float) Math.sin((0.004f * frot))), 15, 15);
				GLlib.renderSphere();
				
				GL11.glColor3f(1.0f, 0.0f, 0.0f);

				GLlib.drawLine(0.0f,0.0f,0.0f,0.0f,0.0f,1.5f);
				GLlib.drawLine(0.0f,0.0f,0.0f,0.0f,1.5f,0.0f);
				//GLlib.drawTriangle();
				canvas.swapBuffers();
				display.timerExec(5, this);
				// try {
				// //wait for Producer to put value
				// wait(50);
				// } catch (InterruptedException e) { }
				// display.asyncExec(this);
			}
		}
	}
}
