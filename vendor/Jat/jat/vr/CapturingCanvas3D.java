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
 File: CapturingCanvas3D.java

 University of Applied Science Berne,HTA-Biel/Bienne,
 Computer Science Department.

 Diploma thesis J3D Solar System Simulator
 Originally written by Marcel Portner & Bernhard Hari (c) 2000

 CVS - Information :

 $Header: /cvsroot/jat/jat/jat/vr/CapturingCanvas3D.java,v 1.8 2002/10/31 17:43:05 tberthold Exp $
 $Author: tberthold $
 $Date: 2005-01-11 13:33:39 -0500 (Tue, 11 Jan 2005) $
 $State: Exp $

*/

/**
 * Class CapturingCanvas3D, using the instructions from the Java3D
 * FAQ pages on how to capture a still image in jpeg format.
 *
 * A capture button would call a method that looks like
 *
 * <pre>
 *  public static void captureImage(CapturingCanvas3D MyCanvas3D) {
 *    MyCanvas3D.writeJPEG_ = true;
 *    MyCanvas3D.repaint();
 *  }
 * </pre>
 *
 * Peter Z. Kunszt
 * Johns Hopkins University
 * Dept of Physics and Astronomy
 * Baltimore MD
 *
 * @author Marcel Portner & Bernhard Hari
 * @version $Revision: 1 $
 */

package jat.vr;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
import javax.media.j3d.*;
import javax.vecmath.Point3f;
import com.sun.image.codec.jpeg.*;

public class CapturingCanvas3D extends Canvas3D
{
	//static String frames_path="frames/";
	public boolean writeJPEG_;
	private int postSwapCount_;
	String frames_path;

	/**
	* Constructor that generate a Canvas3D.
	*
	* @param gc  the GraphicsConfiguration
	*/
	public CapturingCanvas3D(GraphicsConfiguration gc, String frames_path)
	{
		super(gc);
		this.frames_path = frames_path;
		postSwapCount_ = 0;
		writeJPEG_ = false;
	}

	/**
	* Override Canvas3D's postSwap method to save a JPEG of the canvas.
	*/
	public void postSwap()
	{
		if (writeJPEG_)
		{
			System.out.println("Writing JPEG ");
			int dimX = this.getScreen3D().getSize().width;
			int dimY = this.getScreen3D().getSize().height;

			// The raster components need all be set!
			Raster ras =
				new Raster(
					new Point3f(-1.0f, -1.0f, -1.0f),
					Raster.RASTER_COLOR,
					0,
					0,
					dimX,
					dimY,
					new ImageComponent2D(
						ImageComponent.FORMAT_RGB,
						new BufferedImage(dimX, dimY, BufferedImage.TYPE_INT_RGB)),
					null);
			GraphicsContext3D ctx = getGraphicsContext3D();
			ctx.readRaster(ras);

			// Now strip out the image info
			BufferedImage img = ras.getImage().getImage();

			// write that to disk....
			try
			{
				FileOutputStream out = new FileOutputStream(frames_path + "Capture00" + postSwapCount_ + ".jpg");
				JPEGImageEncoder encoder = JPEGCodec.createJPEGEncoder(out);
				JPEGEncodeParam param = encoder.getDefaultJPEGEncodeParam(img);
				param.setQuality(1.0f, false); // x% quality JPEG
				encoder.setJPEGEncodeParam(param);
				encoder.encode(img);
				writeJPEG_ = false;
				out.close();
			} catch (IOException e)
			{
				System.out.println("I/O exception!Maybe path" + frames_path + "does not exist?");
			}
			postSwapCount_++;
		}
	}

	public void takePicture()
	{
		System.out.println("Writing JPEG "+postSwapCount_);
		int dimX = this.getScreen3D().getSize().width;
		int dimY = this.getScreen3D().getSize().height;

		
		// The raster components need all be set!
		Raster ras =
			new Raster(	new Point3f(-1.0f, -1.0f, -1.0f),
				Raster.RASTER_COLOR,
				0,
				0,
				dimX,
				dimY,
				new ImageComponent2D(
					ImageComponent.FORMAT_RGB,
					new BufferedImage(dimX, dimY, BufferedImage.TYPE_INT_RGB)),
				null);
		GraphicsContext3D ctx = getGraphicsContext3D();
		ctx.readRaster(ras);

		// Now strip out the image info
		BufferedImage img = ras.getImage().getImage();
//		BufferedImage img = new BufferedImage(800,600,1);

		
		// write that to disk....
		try
		{
			String filename="Capture";
			if(postSwapCount_<10) filename="Capture00";
			if(postSwapCount_>=10&&postSwapCount_<100) filename="Capture0";
			FileOutputStream out = new FileOutputStream(frames_path + filename + postSwapCount_ + ".jpg");
			JPEGImageEncoder encoder = JPEGCodec.createJPEGEncoder(out);
			JPEGEncodeParam param = encoder.getDefaultJPEGEncodeParam(img);
			param.setQuality(1.0f, false); // x% quality JPEG
			encoder.setJPEGEncodeParam(param);
			encoder.encode(img);
			out.close();
		} catch (IOException e)
		{
			System.out.println("I/O exception!Maybe path" + frames_path + "does not exist?");
		}
		
		postSwapCount_++;
	}

}
