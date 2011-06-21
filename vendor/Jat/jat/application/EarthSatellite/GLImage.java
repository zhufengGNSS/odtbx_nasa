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
import java.nio.*;
import java.io.*;
import java.awt.*;
import java.awt.image.*;
import java.net.URL;

/**
 * Loads an Image from file, stores pixels as ARGB int array, and RGBA ByteBuffer.
 * An alternate constructor creates a GLImage from a ByteBuffer containing
 * pixel data.
 * <P>
 * Static functions are included to load, flip and convert pixel arrays.
 * <P>
 * napier at potatoland dot org
 */

public class GLImage {
    public static final int SIZE_BYTE = 1;
    int h = 0;
    int w = 0;
    ByteBuffer pixelBuffer = null;   // to store bytes in GL_RGBA format
    int[] pixels = null;
    Image image = null;


    public GLImage() {
    }

    /**
     * Load pixels from an image file.
     * @param imgName
     */
    public GLImage(String imgName)
    {
        loadImage(imgName);
    }

    /**
     * Store pixels passed in a ByteBuffer.
     * @param pixels
     * @param w
     * @param h
     */
    public GLImage(ByteBuffer pixels, int w, int h) {
        this.pixelBuffer = pixels;
        this.pixels = null;
        this.image = null;  // image is not loaded from file
        this.h = h;
        this.w = w;
    }

    /**
     * return true if image has been loaded successfully
     * @return
     */
    public boolean isLoaded()
    {
        return (image != null);
    }

    /**
     * Flip the image pixels vertically
     */
    public void flipPixels()
    {
        pixels = GLImage.flipPixels(pixels, w, h);
    }

    /**
     * Scale the image and regenerate pixel array
     */
    /*
    public void resize(float scaleW, float scaleH)
    {
        if (image != null) {
            image = image.getScaledInstance(
                //(int)((float)w*scaleW),
                //(int)((float)h*scaleH),
                600, -1,
                Image.SCALE_SMOOTH);
            if (image != null) {
                w = image.getWidth(null);
                h = image.getHeight(null);
                pixels = getImagePixels(); // pixels in default Java ARGB format
                pixelBuffer = convertImagePixels(pixels, w, h, true); // convert to RGBA bytes
            }
        }
    }
    */

    /**
     * Load an image file and hold its width/height.
     * @param imgName
     */
    public void loadImage(String imgName) {
        Image tmpi = loadImageFromFile(imgName);
        if (tmpi != null) {
            w = tmpi.getWidth(null);
            h = tmpi.getHeight(null);
            image = tmpi;
            pixels = getImagePixels();  // pixels in default Java ARGB format
            pixelBuffer = convertImagePixels(pixels,w,h,true);  // convert to RGBA bytes
            System.out.println("GLImage: loaded " + imgName + ", width=" + w + " height=" + h);
        }
        else {
            System.out.println("GLImage: FAILED TO LOAD IMAGE " + imgName);
            image = null;
            pixels = null;
            pixelBuffer = null;
            h = w = 0;
        }
    }

    /**
     * Return the image pixels in default Java int ARGB format.
     * @return
     */
    public int[] getImagePixels()
    {
        if (pixels == null && image != null) {
            pixels = new int[w * h];
            PixelGrabber pg = new PixelGrabber(image, 0, 0, w, h, pixels, 0, w);
            try {
                pg.grabPixels();
            }
            catch (Exception e) {
                System.out.println("Pixel Grabbing interrupted!");
                return null;
            }
        }
        return pixels;
    }

    /**
     * return int array containing pixels in ARGB format (default Java byte order).
     */
    public int[] getPixelsARGB()
    {
        return pixels;
    }

    /**
     * return ByteBuffer containing pixels in RGBA format (commmonly used in OpenGL).
     */
    public ByteBuffer getPixelsRGBA()
    {
        return pixelBuffer;
    }

    //========================================================================
    //
    // Utility functions to prepare pixels for use in OpenGL
    //
    //========================================================================

    /**
     * Flip an array of pixels vertically
     * @param imgPixels
     * @param imgw
     * @param imgh
     * @return int[]
     */
    public static int[] flipPixels(int[] imgPixels, int imgw, int imgh)
    {
        int[] flippedPixels = null;
        if (imgPixels != null) {
            flippedPixels = new int[imgw * imgh];
            for (int y = 0; y < imgh; y++) {
                for (int x = 0; x < imgw; x++) {
                    flippedPixels[ ( (imgh - y - 1) * imgw) + x] = imgPixels[ (y * imgw) + x];
                }
            }
        }
        return flippedPixels;
    }

    /**
     * Convert ARGB pixels to a ByteBuffer containing RGBA pixels.<BR>
     * Can be drawn in ORTHO mode using:<BR>
     *         GL.glDrawPixels(imgW, imgH, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, byteBuffer); <BR>
     * If flipVertically is true, pixels will be flipped vertically (for OpenGL coord system).
     * @param imgFilename
     * @return ByteBuffer
     */
    public static ByteBuffer convertImagePixels(int[] jpixels, int imgw, int imgh, boolean flipVertically) {
        byte[] bytes;     // will hold pixels as RGBA bytes
        if (flipVertically) {
            jpixels = flipPixels(jpixels, imgw, imgh); // flip Y axis
        }
        bytes = convertARGBtoRGBA(jpixels);
        return allocBytes(bytes);  // convert to ByteBuffer and return
    }

    /**
     * Convert pixels from java default ARGB int format to byte array in RGBA format.
     * @param jpixels
     * @return
     */
    public static byte[] convertARGBtoRGBA(int[] jpixels)
    {
        byte[] bytes = new byte[jpixels.length*4];  // will hold pixels as RGBA bytes
        int p, r, g, b, a;
        int j=0;
        for (int i = 0; i < jpixels.length; i++) {
            int outPixel = 0x00000000; // AARRGGBB
            p = jpixels[i];
            a = (p >> 24) & 0xFF;  // get pixel bytes in ARGB order
            r = (p >> 16) & 0xFF;
            g = (p >> 8) & 0xFF;
            b = (p >> 0) & 0xFF;
            bytes[j+0] = (byte)r;  // fill in bytes in RGBA order
            bytes[j+1] = (byte)g;
            bytes[j+2] = (byte)b;
            bytes[j+3] = (byte)a;
            j += 4;
        }
        return bytes;
    }


    //========================================================================
    // Utility functions to load file into byte array
    // and create Image from bytes.
    //========================================================================

    /**
     * Same function as in GLApp.java.  Allocates a ByteBuffer to hold the given
     * array of bytes.
     *
     * @param bytearray
     * @return  ByteBuffer containing the contents of the byte array
     */
    public static ByteBuffer allocBytes(byte[] bytearray) {
        ByteBuffer bb = ByteBuffer.allocateDirect(bytearray.length * SIZE_BYTE).order(ByteOrder.nativeOrder());
        bb.put(bytearray).flip();
        return bb;
    }

    /**
     * Load an image from file.  Avoids the flaky MediaTracker/ImageObserver headache.
     * Assumes that the file can be loaded quickly from the local filesystem, so
     * does not need to wait in a thread.  If it can't find the file in the
     * filesystem, will try loading from jar file.  If not found will return
     * null.
     * <P>
     * @param imgName
     */
    public static Image loadImageFromFile(String imgName) {
        byte[] imageBytes = getBytesFromFile(imgName);
        Image tmpi = null;
        int numTries = 20;
        if (imageBytes == null) {
            // couldn't read image from file: try JAR file
            //URL url = getClass().getResource(filename);   // for non-static class
            URL url = GLImage.class.getResource(imgName);
            if (url != null) {
                // found image in jar: load it
                tmpi = Toolkit.getDefaultToolkit().createImage(url);
                // Wait up to two seconds to load image
                int wait = 200;
                while (tmpi.getHeight(null) < 0 && wait > 0) {
                    try {
                        Thread.sleep(10);
                    }
                    catch (Exception e) {}
                }
            }
        }
        else {
            tmpi = Toolkit.getDefaultToolkit().createImage(imageBytes, 0, imageBytes.length);
            while (tmpi.getWidth(null) < 0 && numTries-- > 0) {
                try { Thread.sleep(100); }
                catch( InterruptedException e ) {System.out.println(e);}
            }
            while (tmpi.getHeight(null) < 0 && numTries-- > 0) {
                try { Thread.sleep(100); }
                catch( InterruptedException e ) {System.out.println(e);}
            }
        }
        return tmpi;
    }

    public static Image loadImageFromFile_ORIG(String imgName) {
        byte[] imageBytes = getBytesFromFile(imgName);
        Image tmpi = null;
        int numTries = 20;
        if (imageBytes != null) {
            tmpi = Toolkit.getDefaultToolkit().createImage(imageBytes, 0, imageBytes.length);
            while (tmpi.getWidth(null) < 0 && numTries-- > 0) {
                try { Thread.sleep(100); }
                catch( InterruptedException e ) {System.out.println(e);}
            }
            while (tmpi.getHeight(null) < 0 && numTries-- > 0) {
                try { Thread.sleep(100); }
                catch( InterruptedException e ) {System.out.println(e);}
            }
        }
        return tmpi;
    }

    /**
     * Given name of file, return entire file as a byte array.
     * @param filename
     * @return
     */
    public static byte[] getBytesFromFile(String filename)
    {
        File f = new File(filename);
        byte[] bytes = null;
        if (f.exists()) {
            try {
                bytes = getBytesFromFile(f);
            }
            catch (Exception e) {
                System.out.println("getBytesFromFile() exception: " + e);
            }
        }
        return bytes;
    }

    /**
     * Given File object, returns the contents of the file as a byte array.
     */
    public static byte[] getBytesFromFile(File file) throws IOException {
        byte[] bytes = null;
        if (file != null) {
            InputStream is = new FileInputStream(file);
            long length = file.length();
            // Can't create an array using a long type.
            // Before converting to an int type, check
            // to ensure that file is not larger than Integer.MAX_VALUE.
            if (length > Integer.MAX_VALUE) {
                System.out.println("getBytesFromFile() error: File " + file.getName()+ " is too large");
            }
            else {
                // Create the byte array to hold the data
                bytes = new byte[ (int) length];
                int offset = 0;
                int numRead = 0;
                // Read in the bytes
                while (offset < bytes.length
                       && (numRead = is.read(bytes, offset, bytes.length - offset)) >= 0) {
                    offset += numRead;
                }
                // Ensure all the bytes have been read in
                if (offset < bytes.length) {
                    throw new IOException("getBytesFromFile() error: Could not completely read file " + file.getName());
                }
            }
            // Close the input stream and return bytes
            is.close();
        }
        return bytes;
    }

}