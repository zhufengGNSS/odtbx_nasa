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

package jat.util;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.net.*;

public class FileUtil {
  
    /** Directory where output files should go by default. */
    private static File outputDir = deduceOutputDir();
	
	public static String file_separator() {
		return File.separator;
	}
    
    public static String getClassFilePath(String package_name,
        String class_name) {
        
        String path = null;
        try {
          Class c = Class.forName(package_name + "." + class_name);
          path = getClassFilePath(c);
        }
        catch (ClassNotFoundException e) {
          System.err.println("Could not locate class " + package_name + "." + 
              class_name);
          path = null;
        }
        return path;
    }
    
    
    public static String getClassFilePath(Class c) {
		String out = null;

		// check for Windoze
		boolean windoze = (File.separatorChar == '\\');

		URL url = null;
        String fullyQualifiedName = c.getCanonicalName();
        String shortName = c.getSimpleName();
		String classRes = "/" + fullyQualifiedName.replace('.', '/') + ".class";
		url = c.getResource(classRes);

		// Make sure that the class was loaded from a class file
		// in a directory and not from a jar or from an applet
		if (url.getProtocol().equals("file")) {
		  try {
			// Get the path
            // Sometimes JDKs produce null for path and you have to
            // get it from the scheme specific part (which in a 
            // file URI should just be the path)
            out = url.toURI().getPath();
            if (out == null) {
              out = url.toURI().getSchemeSpecificPart();
            }

            // chop off the unneeded stuff
            out = out.substring(0, out.length() - (shortName + ".class").length());             
            if (windoze && out.startsWith("/")) {
                out = out.substring(1);
            } 
		  } catch (URISyntaxException e) {
            // Something has gone awry.  All valid URLs are URIs.
            throw new RuntimeException("Could not find " + fullyQualifiedName + 
                " in classpath.");
		  }

          out = out.replace('/', File.separatorChar);
        }
        else {
          // Not from a directory.  Return null.
          out = null;
        }
		return out;
	}
    
    /**
     * Returns the default directory to put output files.  Usually is jat/sim/output,
     * but when JAT is deployed as a jar, it is $TEMP/jat.  Not sure what it does as
     * an applet.
     * @return default directory to put output files
     */
    public static File getOutputDir() {
      return outputDir;
    }
    
    /**
     * Sets the default directory where output files should go.  This does not
     * need to be called, but can be called to overide the default location
     * of output files.
     * @param path a path to the directory.  If path is not absolute, will
     * use from the current working directory.
     */
    public static void setOutputDir(String path) {
      outputDir = new File(path);
    }
    
    /**
     * Determines what the default directory should be.  If JAT is setup
     * as a source/class tree, will just use a standard directory within
     * the tree.  If not, will try to use the temp directory.
     * @return the directory to store output files
     */
    private static File deduceOutputDir() {
      // First we need to see where JAT came from, so we find out where this
      // class came from and back up the tree to the root.
      File outputDir = null;
      String fileUtilDir = getClassFilePath(FileUtil.class);
      if (fileUtilDir != null) {
        // We back up one directory for every "." in the class name.
        File rootDir = new File(fileUtilDir);
        String name = FileUtil.class.getCanonicalName();
        for(int ct = name.indexOf("."); ct >= 0; ct = name.indexOf(".", ct+1)) {
          rootDir = rootDir.getParentFile();
        }
        outputDir = rootDir;
      }
      else {
        // JAT is not deployed in a directory tree.  Use the working directory.
        outputDir = new File(".");
      }
      return outputDir;
    }
    
    /**
     * Determines where the file should reside and opens an output stream
     * onto it.  If file path is relative, will treat the path as relative to
     * JAT root directory (or working directory if JAT is not a directory) 
     * @param path absolute or relative.  "-" is a special case indicating
     * standard out (its a Unix thing)
     * @return stream that can write to the output file
     */
    public static OutputStream openOutputFile(String path)
      throws IOException {
      // We open a writer to the output media.
      // If fname is null, we open a writer to System.out
      // If fname is a filename with no directory names at all
      //   (or if it specifies the default directory)
      //   we put it in the default directory.
      // Otherwise, we open a writer to the file and if not
      //   an absolute file name will create it relative to
      //   the working directory.
      
      OutputStream strm = null;
      if (path.equals("-")) {
        strm = System.out;
      }
      else if (new File(path).isAbsolute()) {
        strm = new FileOutputStream(path);
      }
      else {
        path = outputDir.getAbsolutePath() + File.separator + path;
        strm = new FileOutputStream(path);
      }
      return strm;
    }
    
    /**
	 * This method searches a given directory path for a case-insensitive match
	 * to the given filename.
	 * 
	 * @param path
	 *            The path to search (no case in-sensitive search is performed
	 *            on the path itself)
	 * @param filename
	 *            The filename to perform a case-insensitive search on
	 * @return The String that matches the found filename on the filesystem.
	 *         Note, this does not have to be either the given filename or the
	 *         string used during the search, only what is reported by the
	 *         filesystem. A null String indicates a failure to find a match at
	 *         all (a bad path or no matching filename).
	 */
	public static String findCaseInsensitiveFilename(final String path,
			final String filename) {

		String retval = null; // default to failure

		File p = new File(path);
		if (!p.exists()) {
			return retval; // no such path
		}

		File files[] = p.listFiles();
		if (files == null) {
			return retval; // no files in path
		}

		for (int i = 0; i < files.length; i++) {
			if (files[i] != null) {
				String f = files[i].getName();
				if (filename.compareToIgnoreCase(f) == 0) {
					retval = f;
					break;
				}
			}
		}
		return retval;
	}
}