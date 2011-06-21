/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2006 United States Government as represented by the
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
package jat.sim;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Properties;
import jat.util.FileUtil;

/**
 * TODO Javadoc
 * @author Richard C. Page III
 */
public class SimProperties {

	private Properties properties;
	private String propertiesFilename;
	private String propertiesDescription;
	private static String dir = "input";
	
//	  private long delay = 20;
//    private long countDown = 5*60;
//    private int maxPlayers = 100;
//    private int maxCards = 3;
//
//    private String delayName = "ball.delay";
//    private String countDownName = "count.down";
//    private String maxPlayersName = "max.players";
//    private String maxCardsName = "max.cards";
	
	//* Spacecraft
	public double[] r0 = new double[3];
	private String[] sr0 = {"x","y","z"};
	public double[] v0 = new double[3];
	private String[] sv0 = {"xdot","ydot","zdot"};
	public String sc_name;
	private String ssc_name = "spacecraft_name";
	public double mass;
	private String smass = "spacecraft_mass";
	public double coef_drag;
	private String scoef_drag = "spacecraft_drag_coefficient";
	public double coef_reflect;
	private String scoef_reflect = "spacecraft_reflectivity_coefficient";
	public double area;
	private String sarea = "spacecraft_cross-section";
	
	public double[] target_rf = new double[3];
	private String[] starget_rf = {"target_x","target_y","target_z"};
	public double[] target_vf = new double[3];
	private String[] starget_vf = {"target_xdot","target_ydot","target_zdot"};
	
	//* Simulation
	public double mjd_utc;
	private String smjd_utc = "mjd_utc";
	public double t0;
	private String st0 = "t0";
	public double tf;
	private String stf = "tf";
	public double stepsize;
	private String sstepsize = "stepsize";
	public double thinning;
	private String sthinning = "thinning";
	//* Force Models
	private boolean use_moon;
	private String suse_moon = "moon";
	private boolean use_sun;
	private String suse_sun = "sun";
	private boolean use_JGM3;
	private String suse_JGM3 = "JGM3";
	private boolean[] use_planet = new boolean[8];
	private String[] suse_planet = 
		{"mercury","venus","earth","mars","jupiter","saturn","uranus","neptune"};
		
	
	
	/**
	 * Default
	 */
	public SimProperties(){
		properties = new Properties();
		propertiesFilename = "default.dat";
		propertiesDescription = "default simulation file";
		setDefaults(properties);
		updateSettingsFromProperties();
	}
	/**
	 * Get Default properties and store to the given filename rather than default file.
	 * @param propertiesFilename filename to store data
	 * @param propertiesDescription description of the simulation data
	 */
	public SimProperties(String propertiesFilename, String propertiesDescription) {
		this.propertiesFilename = propertiesFilename;
		this.propertiesDescription = propertiesDescription;
		properties = new Properties();
		setDefaults(properties);
	}
	/**
	 * Get properties from a file.  Loads the properties from the given file.  The file
	 * must be in the jat.sim.input directory.
	 * @param propFile the filename of the data file
	 */
	public SimProperties(String propFile){
		propertiesFilename = propFile;
		loadProperties();
		
	}
	
	protected void setDefaults(Properties defaults) {
		defaults.put(sr0[0], "-4453783.586");
		defaults.put(sr0[1],"-5038203.756");
		defaults.put(sr0[2],"-426384.456");
		defaults.put(sv0[0],"3831.888");
		defaults.put(sv0[1],"-2887.221");
		defaults.put(sv0[2],"-6018.232");
		defaults.put(ssc_name,"spacecraft");
		defaults.put(smass,"1000");
		defaults.put(scoef_drag,"2.2");
		defaults.put(scoef_reflect,"1.2");
		defaults.put(sarea,"20");
		defaults.put(starget_rf[0],"0");
		defaults.put(starget_rf[1],"0");
		defaults.put(starget_rf[2],"0");
		defaults.put(starget_vf[0],"0");
		defaults.put(starget_vf[1],"0");
		defaults.put(starget_vf[2],"0");
		defaults.put(smjd_utc,"53157.5");  //* June 1, 2004 12:00 UTC
		defaults.put(st0,"0.0");
		defaults.put(stf,"86400.0");
		defaults.put(sstepsize,"60");
		defaults.put(sthinning,"12");
		defaults.put(suse_moon,"false");
		defaults.put(suse_sun,"false");
		defaults.put(suse_JGM3,"false");
		for(int i=0; i<8; i++){
			if(i==2)
				defaults.put(suse_planet[i],"true");
			else
				defaults.put(suse_planet[i],"false");
		}		
	}
	
	public void loadProperties(){
    	properties = new Properties();
    	FileInputStream in = null;
        try {
            String filesep = System.getProperty("file.separator");
            String basedir = FileUtil.getClassFilePath("jat.sim","SimProperties");
            in = new FileInputStream(basedir + filesep + dir + filesep + propertiesFilename);
            properties.load(in);

        } catch (java.io.FileNotFoundException e) {
            in = null;
            System.err.println("Can't find properties file. " +
                                "Using defaults.");
        } catch (java.io.IOException e) {
            System.err.println("Can't read properties file. " +
                                "Using defaults.");
        } finally {
            if (in != null) {
                try { in.close(); } catch (java.io.IOException e) { }
                in = null;
            }
        }

        updateSettingsFromProperties();
    }
    
    public void saveProperties() {
    	saveProperties(this.propertiesFilename);
    }
    
    public void saveProperties(String filename){
    	updatePropertiesFromSettings();

//      if (DEBUG) {
//          System.out.println("Just set properties: "
//                             + propertiesDescription);
//          System.out.println(toString());
//      }

      FileOutputStream out = null;

      try {
          String basedir = FileUtil.getClassFilePath("jat.sim","SimProperties");
          String filesep = System.getProperty("file.separator");
          out = new FileOutputStream(basedir +filesep+dir+ filesep + filename);
          properties.save(out, propertiesDescription);
      } catch (java.io.IOException e) {
          System.err.println("Can't save properties. " +
  			    "Oh well, it's not a big deal.");
      } finally {
          if (out != null) {
  	    try { out.close(); } catch (java.io.IOException e) { }
  	    out = null;
          }
      }
    	
    }

//    protected void getParameters() {
//    	Properties defaults = new Properties();
//    	FileInputStream in = null;
//    	
//    	setDefaults(defaults);
//    	
//    	properties = new Properties(defaults);
//    	
//    	try {
//    		String folder = System.getProperty("user.home");
//    		String filesep = System.getProperty("file.separator");
//    		in = new FileInputStream(folder
//    				+ filesep
//					+ propertiesFilename);
//    		properties.load(in);
//    		
//    	} catch (java.io.FileNotFoundException e) {
//    		in = null;
//    		System.err.println("Can't find properties file. " +
//    		"Using defaults.");
//    	} catch (java.io.IOException e) {
//    		System.err.println("Can't read properties file. " +
//    		"Using defaults.");
//    	} finally {
//    		if (in != null) {
//    			try { in.close(); } catch (java.io.IOException e) { }
//    			in = null;
//    		}
//    	}
//
//        updateSettingsFromProperties();
//
//    }
    
//    protected void saveParameters() {
//    	
//    	updatePropertiesFromSettings();
//    	
////    	if (DEBUG) {
////    		System.out.println("Just set properties: " + propertiesDescription);
////    		System.out.println(toString());
////    	}
//    	
//    	FileOutputStream out = null;
//    	
//    	try {
//    		String folder = System.getProperty("user.home");
//    		String filesep = System.getProperty("file.separator");
//    		out = new FileOutputStream(folder 
//    				+ filesep
//					+ propertiesFilename);
//    		properties.save(out, propertiesDescription);
//    	} catch (java.io.IOException e) {
//    		System.err.println("Can't save properties. " +
//    		"Oh well, it's not a big deal.");
//    	} finally {
//    		if (out != null) {
//    			try { out.close(); } catch (java.io.IOException e) { }
//    			out = null;
//    		}
//    	}
//    }
    

    private void updatePropertiesFromSettings(){
    	for(int i=0; i<3; i++){
    		properties.put(sr0[i], new Double(r0[i]).toString());
    		properties.put(sv0[i], new Double(v0[i]).toString());
    		properties.put(starget_rf[i], new Double(target_rf[i]).toString());
    		properties.put(starget_vf[i], new Double(target_vf[i]).toString());
    	}
    	properties.put(ssc_name,sc_name);
    	properties.put(smass, new Double(mass).toString());
    	properties.put(scoef_drag, new Double(coef_drag).toString());
    	properties.put(scoef_reflect, new Double(coef_reflect).toString());
    	properties.put(sarea, new Double(area).toString());
    	properties.put(smjd_utc, new Double(mjd_utc).toString());  
		properties.put(st0, new Double(t0).toString());
		properties.put(stf, new Double(tf).toString());
		properties.put(sstepsize, new Double(stepsize).toString());
		properties.put(sthinning, new Double(thinning).toString());
		properties.put(suse_moon, new Boolean(use_moon).toString());
		properties.put(suse_sun, new Boolean(use_sun).toString());
		properties.put(suse_JGM3, new Boolean(use_JGM3).toString());
		for(int i=0; i<8; i++){
			properties.put(suse_planet[i], new Boolean(use_planet[i]).toString());			
		}		
    }
    
    private void updateSettingsFromProperties(){
    	long ONE_SECOND = 1000;
    	try {
    		for(int i=0; i<3; i++){
    			r0[i] = Double.parseDouble(properties.getProperty(sr0[i]));
    			v0[i] = Double.parseDouble(properties.getProperty(sv0[i]));
    			target_rf[i] = Double.parseDouble(properties.getProperty(starget_rf[i]));
    			target_vf[i] = Double.parseDouble(properties.getProperty(starget_vf[i]));
    		}
    		sc_name = ssc_name;
    		mass = Double.parseDouble(properties.getProperty(smass));
    		coef_drag = Double.parseDouble(properties.getProperty(scoef_drag));
    		coef_reflect = Double.parseDouble(properties.getProperty(scoef_reflect));
    		area = Double.parseDouble(properties.getProperty(sarea));
    		mjd_utc = Double.parseDouble(properties.getProperty(smjd_utc));
    		t0 = Double.parseDouble(properties.getProperty(st0));
    		tf = Double.parseDouble(properties.getProperty(stf));
    		stepsize = Double.parseDouble(properties.getProperty(sstepsize));
    		thinning = Double.parseDouble(properties.getProperty(sthinning));
    		use_moon = Boolean.parseBoolean(properties.getProperty(suse_moon));
    		use_sun = Boolean.parseBoolean(properties.getProperty(suse_sun));
    		use_JGM3 = Boolean.parseBoolean(properties.getProperty(suse_JGM3));
    		for(int i=0; i<8; i++)
    			use_planet[i] = Boolean.parseBoolean(properties.getProperty(suse_planet[i]));    		
    	} catch (NumberFormatException e) {
    		System.err.println("NumberFormatException ");
    	    // we don't care if the property was of the wrong format,
    	    // they've all got default values. So catch the exception
    	    // and keep going.
    	}
    }
    
    public void toggle_use_moon(boolean b){
    	this.use_moon = b;
    }
    public void toggle_use_sun(boolean b){
    	this.use_sun = b;
    }
    public void toggle_use_JGM3(boolean b){
    	this.use_JGM3 = b;
    }
    public void toggle_use_planet(int i, boolean b){
    	this.use_planet[i] = b;
    }
    public void toggle_all_planet(boolean b){
    	for(int i=0; i<8; i++){
    		this.use_planet[i] = b;
    	}
    }
    
	public static void main(String[] args) {
		System.out.println("test");
		SimProperties test = new SimProperties();
		System.out.println("x:    "+test.r0[0]);
		System.out.println("z:    "+test.r0[2]);
		System.out.println("ydot: "+test.v0[1]);
		System.out.println("zdot: "+test.v0[2]);
		System.out.println("jptr: "+test.suse_planet[4]+" "+test.use_planet[4]);
		//test.r0[0] = 42000.2;
		//test.use_planet[2] = true;
		test.saveProperties("default.prop");
		System.out.println("finished");
	}
}
