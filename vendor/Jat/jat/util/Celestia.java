/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2005 United States Government as represented by the
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

import jat.traj.*;
//import jat.matvec.data.Matrix;
//import jat.matvec.data.VectorN;

//import java.lang.Math;
import java.io.*;
import java.text.*;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 * This class provides a way to create the spacecraft files and the position
 * information required to display a time-dependent spacecraft trajectory
 * in Celestia.
 * 
 * Celestia is an open source astronomical visualization application written
 * using C++ and is available at the following location.
 * 
 * The suggested use is as follows:
 *    Celestia celestia = new Celestia("C:/celestia/");
 *    // propagate orbit sim.traj
 *    celestia.set_trajectory(sim.get_traj());
 *    String name1 = "my_test_run_filename";
 *    String name2 = "my_test_run"
 *    celestia.write_trajectory(name,name,sim.mjd_utc_start+2400000.5);
 * 
 * Open Celestia, goto Earth, press enter and look for 'name2'
 * 
 * @link http://www.shatters.net/celestia
 * @link http://sourceforge.net/projects/celestia/
 * 
 * @author Richard C. Page III
 *
 */
public class Celestia {

    /**
     * Celestia's install directory
     */
    public static String celestia_home;
    /**
     * Double array holding the position vectors [km] as expressed in the 
     * inertial frame about the central body (presumably J2000 ECI).
     */
    private double[][] r;
    /**
     * Parallel array holding the time index for the corresponding
     * position vectors in Julian Date of UTC.
     */
    private double[] time;
    
    /**
     * Default Constructor.  Automatically searches for the Celestia home directory.
     * Currently, only windows is supported (sincere apologies).
     */
    public Celestia(){
    	this(null);
    }
    
    /**
     * Constructor
     * @param dir The install directory of Celestia (e.g. "C:/Celestia/")
     */
    public Celestia(String dir) {
        // We assume that every installation of Celestia has a
        // celestia.cfg in the root directory.
        final String targetName = "celestia.cfg";
        File targetFile = (dir == null ? null : new File(dir, targetName));
        if ((targetFile == null) || !targetFile.exists()) {
          
          // Try looking for it in the path
          targetFile = searchPath(targetName);
          if ((targetFile == null) || !targetFile.exists()) {
            
            // Other possibilities.
            String[] dirs = {System.getenv("CELESTIA_HOME"),
                "C:/Program Files/Celestia",
                "/opt/kde3/share/apps/celestia"};
            boolean found = false;
            for (int ctr=0; (ctr < dirs.length) && !found; ++ctr) {
              targetFile = (dirs[ctr] == null ?
                  null : new File(dirs[ctr], targetName));
              found = ((targetFile != null) && (targetFile.exists()));
            }
          }

          if ((targetFile == null) || (!targetFile.exists())) {
            String msg = "Cannot determine location of Celestia.";
            if (dir != null) {
              msg += "  Not found in " + dir;
            }
            throw new RuntimeException(msg);
          }
        }
    
        celestia_home = targetFile.getParentFile().getAbsolutePath() + 
          File.separator;
    }

    /**
     * Create an XYZ position file in Celestia's "data" directory using the
     * current values stored to the class members 'r' and 'time'.
     * @param filename The filename for the XYZ file (e.g. "STS_115")
     */
    public void write_xyz(String filename){
        String file = celestia_home+"data/"+filename+".xyz";
        File arg = new File(file);
        write_xyz(arg,this.time,this.r,this.time.length);
    }
    
    /**
     * Helper method for write_xyz(filename)
     * 
     * @param filename The filename for the XYZ file (e.g. "STS_115")
     * @param time The parallel array of time values
     * @param r The parallel array of position vectors
     * @param size The number of ephemeris points
     */
    public void write_xyz(File filename,double[] time, double[][] r, int size) {
		//fprintf(fid,'%12.8f %12.8f %12.8f %12.8f\n',t(k),a_xyz(1,k),a_xyz(2,k),a_xyz(3,k));
		try{
			NumberFormat formatter = new DecimalFormat();
			//* This format is to ensure that the correct precision is printed to Celestia
			//* as well as ensuring that there is no exponential notation
			formatter = new DecimalFormat("######.###############");
			
			BufferedWriter fout = new BufferedWriter(new FileWriter(filename));
			NumberFormat nf = NumberFormat.getInstance();
			nf.setMaximumIntegerDigits(12);
			
			int i=0;
			while(size>i){
				fout.write(""+formatter.format(time[i])+" "+formatter.format(r[i][0])+
							" "+formatter.format(r[i][1])+" "+formatter.format(r[i][2])+"\n");
				i++;
			}
			fout.close();
		}catch(java.io.IOException ioe){System.out.println("EX file io");}
	}

    /**
     * Write the Spacecraft information file in Celestia's "extras" directory.
     * @param name The name of the spacecraft to appear in Celestia (e.g. "STS-115")
     * @param filename The name of the file (e.g. "STS-115")
     * @param jd_begin The start date (julian date) of the spacecraft ephemerides
     * @throws IOException
     */
    public void write_ssc(String name, String filename, double jd_begin) throws IOException {
		File file = new File(celestia_home+"extras/"+filename+".ssc");
		BufferedWriter fout = new BufferedWriter(new FileWriter(file));
			fout.write("# "+name);
			fout.newLine();
			fout.write("#");fout.newLine();
			fout.write("\""+name+"\" \"Sol/Earth\"");
			fout.newLine();
			fout.write("{");fout.newLine();
			fout.write("Class \"spacecraft\"");
			fout.newLine();
			fout.write("Mesh \"cassini.3ds\"");
			fout.newLine();
			fout.write("Radius 0.01");
			fout.newLine();
			fout.newLine();
			fout.write("     Beginning "+(time[0]));//jd_begin);
			fout.newLine();
			fout.write("     Ending "+(time[time.length-1]));//(jd_begin+1000));
			fout.newLine();
			fout.newLine();
			fout.write("     SampledOrbit   \""+filename+".xyz\"");
			fout.newLine();
			fout.write("}");
			fout.newLine();
		fout.close();
	}
    
    /**
     * Write the Spacecraft information file in Celestia's "extras" directory.
     * Formats the file for a heliocentric trajectory.
     * 
     * @param name The name of the spacecraft to appear in Celestia (e.g. "STS-115")
     * @param filename The name of the file (e.g. "STS-115")
     * @param jd_begin The start date (julian date) of the spacecraft ephemerides
     * @throws IOException
     */
    public void write_ssc_heliocentric(String name, String filename, double jd_begin) throws IOException {
		File file = new File(celestia_home+"extras/"+filename+".ssc");
		BufferedWriter fout = new BufferedWriter(new FileWriter(file));
			fout.write("# "+name);
			fout.newLine();
			fout.write("#");fout.newLine();
			fout.write("\""+name+"\" \"Sol\"");
			fout.newLine();
			fout.write("{");fout.newLine();
			fout.write("Class \"spacecraft\"");
			fout.newLine();
			fout.write("Mesh \"cassini.3ds\"");
			fout.newLine();
			fout.write("Radius 0.01");
			fout.newLine();
			fout.newLine();
			fout.write("     Beginning "+(time[0]));//jd_begin);
			fout.newLine();
			fout.write("     Ending "+(time[time.length-1]));//(jd_begin+1000));
			fout.newLine();
			fout.newLine();
			fout.write("     SampledOrbit   \""+filename+".xyz\"");
			fout.newLine();
			fout.write("}");
			fout.newLine();
		fout.close();
	}
    
    /**
     * Reads a tab delimited state file into the class to then write to Celestia.
     * The file should contain the following in each row:
     * time(MJD_UTC)	position_x(km)	position_y(km)	position_z(km)		...other_stuff
     * 
     * @param filename Name of the file
     * @throws IOException
     */
    public void readStateFile(String filename) throws IOException {
        FileReader fr;
		BufferedReader in = null;
		int size = 0;
		try {
			fr = new FileReader(filename);
			in = new BufferedReader(fr);			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		ArrayList lineCollection = new ArrayList();

		// read in the file
		boolean eof = false;
		String str;
		StringTokenizer tok;

		while (!eof) {
			String line;
			try {
				if ((line = in.readLine()) == null) {
					eof = true;
				} else {
					// add to the collection
					lineCollection.add(line);
					size++;
				}
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		}
			
		r = new double[size][3];
		time = new double[size];
		for(int i=0; i<size; i++){
		    str = (String) lineCollection.get(i);
		    tok = new StringTokenizer(str, "\t");
		    time[i] = Double.parseDouble(tok.nextToken());
		    time[i] = time[i] + 2400000.5;
		    r[i][0] = Double.parseDouble(tok.nextToken());
		    r[i][1] = Double.parseDouble(tok.nextToken());
		    r[i][2] = Double.parseDouble(tok.nextToken());
		}
		
		
    }
    /**
     * Read a file into a trajectory and then store to the class for writing.
     * @see jat.traj.Trajectory
     * @param file
     */
    public void get_trajectory(String file){
        Trajectory traj = new Trajectory();
        traj.readFromFile(file);
        time = traj.timeArray();
        for(int i=0; i<traj.npts(); i++){
            time[i] = time[i] + 2400000.5;
        }
        r = traj.positionArray();
        
    }
    /**
     * Obtain time and position data from a Trajectory object and store into the
     * class to then format for Celestia.
     * 
     * Trajectory format should be as follows:
     * time(MJD_UTC)	position_x(km)	position_y(km)	position_z(km)		...other_stuff
     * 
     * @param traj Trajectory
     */
    public void set_trajectory(Trajectory traj_km){
        time = traj_km.timeArray();
        for(int i=0; i<time.length; i++){
            time[i] = time[i] + 2400000.5;
        }
        r = traj_km.positionArray();
    }
    /**
     * Obtain time and position data from a Trajectory object and store into the
     * class to then format for Celestia.  (Converts from meters to kilometers)
     * 
     * Trajectory format should be as follows:
     * time(MJD_UTC)	position_x(km)	position_y(km)	position_z(km)		...other_stuff
     * 
     * @param traj Trajectory
     */    
    public void set_trajectory_meters(Trajectory traj_m, double mjd0){
    	time = traj_m.timeArray();
    	r = traj_m.positionArray();
    	double units = 1;
    	double mjd = 0;
    	if(time[0]==0){
    		units = 86400;
    		mjd = mjd0;
    	}
    	for(int i=0; i<time.length; i++){
    		time[i] = time[i]/units + 2400000.5 + mjd;
    		for(int j=0; j<r[0].length; j++)
    			r[i][j] = r[i][j]/1000.0;
    	}
    }
    /**
     * Writes the appropriate files to Celestia.  The member data for 'r' and 'time' must be set.
     * 
     * @param filename Filename for the two files without extension (filename.xyz and filename.ssc)
     * @param pretty_name This string will appear in Celestia identifying the spacecraft
     * @param jd_begin The start date for the spacecraft.
     * @throws IOException
     */
    public void write_trajectory(String filename, String pretty_name, double jd_begin)
    	throws IOException
    {
        this.write_ssc(pretty_name,filename,jd_begin);
        this.write_xyz(filename);
    }
    
    /**
     * Writes the appropriate files to Celestia for a heliocentric trajectory.  
     * The member data for 'r' and 'time' must be set.
     * 
     * @param filename Filename for the two files without extension (filename.xyz and filename.ssc)
     * @param pretty_name This string will appear in Celestia identifying the spacecraft
     * @param jd_begin The start date for the spacecraft.
     * @throws IOException
     */
    public void write_heliocentric(String filename, String pretty_name, double jd_begin)
    	throws IOException
    {
        this.write_ssc_heliocentric(pretty_name,filename,jd_begin);
        this.write_xyz(filename);
    }

    /**
     * Search through the system PATH for Celestia
     * @param targetName the name of a file that will be found
     * in Celestia (and hopefully not elsewhere)
     * @return the target file.  If Celestia is not in the path
     * this may point to a non-existent file or null.
     */
    private File searchPath(String targetName){
      File targetFile = null;
      String pathStr = System.getenv("PATH");
      if (pathStr != null) {
        StringTokenizer tokenizer = new StringTokenizer(pathStr, File.pathSeparator);
        boolean found = false;
        while (!found && tokenizer.hasMoreTokens()) {
          targetFile = new File(tokenizer.nextToken(), targetName);
          found = targetFile.exists();
        }
      }
      return targetFile;
    }
    
	public static void main(String[] args) throws java.io.IOException {
	    //Celestia.run();
	    Celestia celestia = new Celestia();
	    //celestia.get_trajectory("C:/Code/Jat/jat/test/propagator/output/ISS_variable.txt");
	}
	
	public static void run() throws java.io.IOException {
	    Celestia celestia = new Celestia("C:/Celestia_dev/");
	    //String file[] = {"ISS6"};
	    String file[] = {"Sun-Sync12"};
	    //String file[] = {"GEO30","Molniya23"};
	    //String file[] = {"STK_moon_delim2"};
	    for(int i=0; i<file.length; i++){
		    celestia.readStateFile("C:/Code/Jat/jat/test/propagator/output/"+file[i]+".txt");
		    celestia.write_ssc(file[i]+"_jat",file[i]+"_jat",2453158);
		    celestia.write_xyz(file[i]+"_jat");
	        celestia.readStateFile("C:/STK_Test_Files/delim/"+file[i]+".txt");
		    celestia.write_ssc(file[i]+"_stk",file[i]+"_stk",2453158);
		    celestia.write_xyz(file[i]+"_stk");
//		    celestia.readStateFile("C:/ACME/"+file[i]+".txt");
//		    celestia.write_ssc(file[i],file[i],2453158);
//		    celestia.write_xyz(file[i]);
		    
	    }
	    
	    
	    
	    System.out.println("Finished");
	}


}
