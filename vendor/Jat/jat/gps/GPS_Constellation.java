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
 *
 */
 
package jat.gps;

//import java.io.*;
import java.util.*;
import jat.matvec.data.*;
import jat.spacetime.*;
import jat.util.*;
/**
 * <P>
 * The GPS_Constellation Class provides a model of the GPS constellation
 * derived from a RINEX nav file.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class GPS_Constellation {
    private ArrayList SVList = new ArrayList();
    private GPSTimeFormat earliestEpoch;    
    private GPSTimeFormat latestEpoch;
    private Map SVMap = new HashMap();
    private double [] alpha;
    private double [] beta;
    
    /**
     * Constructor
     * @param filename path and filename for the RINEX nav file
     */
    public GPS_Constellation(String filename) {
        RINEXnav rinex = new RINEXnav(filename);
        SVList = rinex.process();
        this.alpha = rinex.alpha();
        this.beta = rinex.beta();
        earliestEpoch = findEarliestEpoch();
        latestEpoch = findLatestEpoch();
        
        // create the SVMap
        for (int i = 0; i < SVList.size(); i++) {
        	GPS_SV sv = this.getSV(i);
        	Integer prn = new Integer(sv.prn());
        	Integer index = new Integer(i);
        	SVMap.put(prn, index);
        }
    }
    
    /**
     * Print out the GPS constellation characteristics
     */
    public void print(){
        System.out.println("GPS Constellation Report");
        System.out.println("Total Number of SVs: "+this.size());
        earliestEpoch.print("Earliest Epoch");        
        latestEpoch.print("Latest Epoch");       
        for(int i = 0; i < this.size(); i++){
            this.getSV(i).print("SV["+i+"]");
        }
        System.out.println("Iono correction terms");
        VectorN a = new VectorN(this.alpha);
        VectorN b = new VectorN(this.beta);
        System.out.println("alpha: "+a.toString());
        System.out.println("beta:  "+b.toString());
    }
    
    /**
     * Get a GPS SV from the constellation
     * @param index SV index
     * @return GPS_SV corresponding to the index
     */
    public GPS_SV getSV(int index){
        return (GPS_SV)SVList.get(index);
    }
    
    /**
     * Get the SV index for a particular PRN
     * @param prn int containing the PRN
     * @return int containing the SV index
     */
    public int getIndex(int prn){
    	Integer PRN = new Integer(prn);
    	Integer ind = (Integer)SVMap.get(PRN);
    	int index = ind.intValue();
    	return index;
    }   

    /**
     * Get a GPS SV by PRN
     * @param prn int containing the PRN
     * @return GPS_SV corresponding to the prn
     */
    public GPS_SV getPRN(int prn){
    	int index = this.getIndex(prn);
        return (GPS_SV)SVList.get(index);
    }
    
	/**
	 * Add a GPS_SV to the Constellation
	 * @param sv GPS_SV to be added      
     */
    public void addSV(GPS_SV sv){
        SVList.add(sv);
    }
    
    /**
     * Return the size of the constellation
     * @return number of SVs in the constellation
     */
    public int size(){
        return SVList.size();
    }
    
    public double[] alpha() {
    	return this.alpha;
    }
    
    public double[] beta() {
    	return this.beta;
    }
    
    private GPSTimeFormat findLatestEpoch(){
        double out = 0.0;
        for (int i = 0; i < SVList.size(); i++){
            double temp = getSV(i).getTOEmjd();
            if (temp > out) {
                out = temp;
            }
        }
        GPSTimeFormat output = new GPSTimeFormat(out);
        return new GPSTimeFormat(out);
    }
    
    private GPSTimeFormat findEarliestEpoch(){
        double out = 1.0E30;
        for (int i = 0; i < SVList.size(); i++){
            double temp = getSV(i).getTOEmjd();
            if (temp < out) {
                out = temp;
            }
        }
        GPSTimeFormat output = new GPSTimeFormat(out);
        return new GPSTimeFormat(out);
    }
    
    /**
     * Return the earliest epoch
     * @return GPSTimeFormat containing the earliest epoch
     */
    public GPSTimeFormat getEarliestEpoch(){
        return this.earliestEpoch;
    }
    
    /**
     * Return the latest epoch
     * @return GPSTimeFormat containing the latest epoch
     */
    public GPSTimeFormat getLatestEpoch(){
        return this.latestEpoch;
    }
    
    public static void main(java.lang.String args[]) {
    	String fs = FileUtil.file_separator();
    	String dir = FileUtil.getClassFilePath("jat.gps","GPS_Constellation");
        String directory = dir+fs+"navfiles"+fs;
        String filename = "rinex.n";
//        String filename = "predata_2.n";        
        String file = directory + filename;
        GPS_Constellation gps = new GPS_Constellation(file);
        gps.print();

//        double t_mjd = (352769.0)/86400.0;
//        for (int i = 0; i < 29; i++){
//            GPS_SV sv1 = gps.getSV(i);
//            sv1.print("sv1");
//            VectorN r = sv1.rECI(t_mjd);
//            int prn = sv1.getPRN();
//            r.print("GPS position vector for "+prn);
//            t_mjd = t_mjd + 1.0/86400.0;
//        }
    }
}



