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

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import org.lwjgl.opengl.GL11;
import org.lwjgl.opengl.glu.Sphere;

public class GLlib
{
    // Used when allocating native buffers for data types
    public static final int SIZE_DOUBLE = 8;
    public static final int SIZE_FLOAT = 4;
    public static final int SIZE_INT = 4;
    public static final int SIZE_BYTE = 1;

    
    
    
	// Display List ID for a sphere (see makeSphere() and renderSphere())
	static int sphereID = -1;

	/**
	 * Render a unit sphere at the origin, using current color, texture and material. This function calls makeSphere()
	 * to create the sphere in a display list.
	 */
	public static void renderSphere()
	{
		if (sphereID < 0)
		{
			makeSphere();
		}
		GL11.glCallList(sphereID);
	}

	public static void makeSphere()
	{
		Sphere s = new Sphere(); // an LWJGL class for drawing sphere
		sphereID = GL11.glGenLists(1); // Create a Display List
		s.setTextureFlag(true); // generate texture coords
		GL11.glNewList(sphereID, GL11.GL_COMPILE); // Start building a list
		{
			GL11.glRotatef(-90, 1, 0, 0); // rotate around X axis
			s.draw(1, 24, 24); // run GL commands to draw sphere
		}
		GL11.glEndList(); // Done building the display list
	}
    
    
    
    
	public static void drawTorus(float r, float R, int nsides, int rings)
	{
		float ringDelta = 2.0f * (float) Math.PI / rings;
		float sideDelta = 2.0f * (float) Math.PI / nsides;
		float theta = 0.0f, cosTheta = 1.0f, sinTheta = 0.0f;
		for (int i = rings - 1; i >= 0; i--)
		{
			float theta1 = theta + ringDelta;
			float cosTheta1 = (float) Math.cos(theta1);
			float sinTheta1 = (float) Math.sin(theta1);
			GL11.glBegin(GL11.GL_QUAD_STRIP);
			float phi = 0.0f;
			for (int j = nsides; j >= 0; j--)
			{
				phi += sideDelta;
				float cosPhi = (float) Math.cos(phi);
				float sinPhi = (float) Math.sin(phi);
				float dist = R + r * cosPhi;
				GL11.glNormal3f(cosTheta1 * cosPhi, -sinTheta1 * cosPhi, sinPhi);
				GL11.glVertex3f(cosTheta1 * dist, -sinTheta1 * dist, r * sinPhi);
				GL11.glNormal3f(cosTheta * cosPhi, -sinTheta * cosPhi, sinPhi);
				GL11.glVertex3f(cosTheta * dist, -sinTheta * dist, r * sinPhi);
			}
			GL11.glEnd();
			theta = theta1;
			cosTheta = cosTheta1;
			sinTheta = sinTheta1;
		}
	}

	
	public static void drawLine(float x0, float y0,float z0,float x1,float y1,float z1)
	{
		GL11.glBegin(GL11.GL_LINES);
		GL11.glVertex3f(x0, y0, z0);
		GL11.glVertex3f(x1, y1, z1);
		GL11.glEnd();
	}

	public static void drawTriangle()
	{
		// GL11.glBegin(GL11.GL_TRIANGLES);
		GL11.glBegin(GL11.GL_LINES);
		GL11.glVertex3f(-1.0f, -1.0f, 0.0f);
		GL11.glVertex3f(1.0f, -1.0f, 0.0f);
		GL11.glVertex3f(0.0f, 1.0f, 0.0f);
		GL11.glEnd();
		}

	
	
	

	public static void setSpotLight(int GLLightHandle, float[] diffuseLightColor, float[] ambientLightColor,
			float[] position, float[] direction, float cutoffAngle)
	{
		FloatBuffer ltDirection = allocFloats(direction);
		setLight(GLLightHandle, diffuseLightColor, ambientLightColor, position);
		GL11.glLightf(GLLightHandle, GL11.GL_SPOT_CUTOFF, cutoffAngle); // width of the beam
		GL11.glLight(GLLightHandle, GL11.GL_SPOT_DIRECTION, ltDirection); // which way it points
		GL11.glLightf(GLLightHandle, GL11.GL_CONSTANT_ATTENUATION, 2F); // how light beam drops off
		// GL11.glLightf(GLLightHandle, GL11.GL_LINEAR_ATTENUATION, .5F); // how light beam drops off
		// GL11.glLightf(GLLightHandle, GL11.GL_QUADRATIC_ATTENUATION, .5F); // how light beam drops off
	}

	
	/**
	 * Set the color of a 'positional' light (a light that has a specific position within the scene). <BR>
	 * 
	 * Pass in an OpenGL light number (GL11.GL_LIGHT1), the 'Diffuse' and 'Ambient' colors (direct light and reflected
	 * light), and the position.<BR>
	 * 
	 * @param GLLightHandle
	 * @param diffuseLightColor
	 * @param ambientLightColor
	 * @param position
	 */
	public static void setLight(int GLLightHandle, float[] diffuseLightColor, float[] ambientLightColor,
			float[] position)
	{
		FloatBuffer ltDiffuse = allocFloats(diffuseLightColor);
		FloatBuffer ltAmbient = allocFloats(ambientLightColor);
		FloatBuffer ltPosition = allocFloats(position);
		GL11.glLight(GLLightHandle, GL11.GL_DIFFUSE, ltDiffuse); // color of the direct illumination
		GL11.glLight(GLLightHandle, GL11.GL_SPECULAR, ltDiffuse); // color of the highlight
		GL11.glLight(GLLightHandle, GL11.GL_AMBIENT, ltAmbient); // color of the reflected light
		GL11.glLight(GLLightHandle, GL11.GL_POSITION, ltPosition);
		GL11.glEnable(GLLightHandle); // Enable the light (GL_LIGHT1 - 7)
		// GL11.glLightf(GLLightHandle, GL11.GL_QUADRATIC_ATTENUATION, .005F); // how light beam drops off
	}
	
		
	

    /**
	 * Load an image from the given file and return a GLImage object.
	 * 
	 * @param imgFilename
	 * @return
	 */
    public static GLImage loadImage(String imgFilename) {
        GLImage img = new GLImage(imgFilename);
        if (img.isLoaded()) {
            return img;
        }
        return null;
    }

    /**
	 * Create a texture from the given image.
	 */
    public static int makeTexture(GLImage textureImg)
    {
        if ( textureImg == null ) {
            return 0;
        }
        else {
            return makeTexture(textureImg.pixelBuffer, textureImg.w, textureImg.h);
        }
    }

    /**
	 * Create a texture from the given pixels in RGBA format. Set the texture to repeat in both directions and use
	 * LINEAR for magnification.
	 * 
	 * @return the texture handle
	 */
    public static int makeTexture(ByteBuffer pixels, int w, int h)
    {
        // get a new empty texture
        int textureHandle = allocateTexture();
        // 'select' the new texture by it's handle
        GL11.glBindTexture(GL11.GL_TEXTURE_2D,textureHandle);
        // set texture parameters
        GL11.glTexParameteri(GL11.GL_TEXTURE_2D, GL11.GL_TEXTURE_WRAP_S, GL11.GL_REPEAT);
        GL11.glTexParameteri(GL11.GL_TEXTURE_2D, GL11.GL_TEXTURE_WRAP_T, GL11.GL_REPEAT);
        GL11.glTexParameteri(GL11.GL_TEXTURE_2D, GL11.GL_TEXTURE_MAG_FILTER, GL11.GL_LINEAR); // GL11.GL_NEAREST);
        GL11.glTexParameteri(GL11.GL_TEXTURE_2D, GL11.GL_TEXTURE_MIN_FILTER, GL11.GL_LINEAR); // GL11.GL_NEAREST);
        // Create the texture from pixels
        GL11.glTexImage2D(GL11.GL_TEXTURE_2D, 0, GL11.GL_RGBA, w, h, 0, GL11.GL_RGBA, GL11.GL_UNSIGNED_BYTE, pixels);
        return textureHandle;
    }

    /**
	 * Allocate a texture (glGenTextures) and return the handle to it.
	 */
    public static int allocateTexture()
    {
        IntBuffer textureHandle = allocInts(1);
        GL11.glGenTextures(textureHandle);
        return textureHandle.get(0);
    }
	

	
	
    // ========================================================================
    // Buffer allocation functions
    //
    // These functions create and populate the native buffers used by LWJGL.
    // ========================================================================

    public static ByteBuffer allocBytes(int howmany) {
        return ByteBuffer.allocateDirect(howmany * SIZE_BYTE).order(ByteOrder.nativeOrder());
    }

    public static IntBuffer allocInts(int howmany) {
        return ByteBuffer.allocateDirect(howmany * SIZE_INT).order(ByteOrder.nativeOrder()).asIntBuffer();
    }

    public static FloatBuffer allocFloats(int howmany) {
        return ByteBuffer.allocateDirect(howmany * SIZE_FLOAT).order(ByteOrder.nativeOrder()).asFloatBuffer();
    }

    public static DoubleBuffer allocDoubles(int howmany) {
        return ByteBuffer.allocateDirect(howmany * SIZE_DOUBLE).order(ByteOrder.nativeOrder()).asDoubleBuffer();
    }

    public static ByteBuffer allocBytes(byte[] bytearray) {
        ByteBuffer bb = ByteBuffer.allocateDirect(bytearray.length * SIZE_BYTE).order(ByteOrder.nativeOrder());
        bb.put(bytearray).flip();
        return bb;
    }

    public static IntBuffer allocInts(int[] intarray) {
        IntBuffer ib = ByteBuffer.allocateDirect(intarray.length * SIZE_FLOAT).order(ByteOrder.nativeOrder()).asIntBuffer();
        ib.put(intarray).flip();
        return ib;
    }

    public static FloatBuffer allocFloats(float[] floatarray) {
        FloatBuffer fb = ByteBuffer.allocateDirect(floatarray.length * SIZE_FLOAT).order(ByteOrder.nativeOrder()).asFloatBuffer();
        fb.put(floatarray).flip();
        return fb;
    }
		
	
	
	
}
