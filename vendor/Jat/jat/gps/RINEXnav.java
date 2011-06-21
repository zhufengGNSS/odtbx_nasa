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
import java.io.*;
import java.util.*;
import jat.spacetime.*;

/**
 * <P>
 * The RINEXnav Class reads in a RINEX format GPS navigation file.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class RINEXnav {
    
    private String filename;
    private ArrayList PRNlist = new ArrayList();
    private double [] alpha = new double[4];
    private double [] beta = new double[4];
    
    /** Creates a new instance of RINEXparser */
    public RINEXnav(String fname) {
        this.filename = fname;
    }
    
    private double parse(String str){
//        System.out.println("Parsing: "+str);
    	try{
        StringTokenizer tok = new StringTokenizer(str, "D");
        double mantissa = Double.parseDouble(tok.nextToken());
        double exponent = Double.parseDouble(tok.nextToken());
        double multiplier = Math.pow(10.0, exponent);
        double out = mantissa * multiplier;
//        System.out.println("Returning: "+out);
        return out;
    	}catch(Exception e){
    		StringTokenizer tok = new StringTokenizer(str, "E");
            double mantissa = Double.parseDouble(tok.nextToken());
            double exponent = Double.parseDouble(tok.nextToken());
            double multiplier = Math.pow(10.0, exponent);
            double out = mantissa * multiplier;
//            System.out.println("Returning: "+out);
            return out;
    	}
    }
    
    private double[] parseLine(String str){
        double [] out = new double[4];
        for (int i = 0; i < 4; i++) {
            int start = i*19 + 3;
            int end = i*19 + 22;
            String string = str.substring(start, end);
            out[i] = parse(string);
        }
        return out;
    }
    
    private double[] parseFirstLine(String str){
        double [] out = new double[3];
        for (int i = 0; i < 3; i++) {
            int start = (i+1)*19 + 3;
            int end = (i+1)*19 + 22;
            String string = str.substring(start, end);
            out[i] = parse(string);
        }
        return out;
    }
    
    private int parsePRN(String str){
        StringTokenizer tok = new StringTokenizer(str, " ");
        String prns = tok.nextToken();
        int prn = Integer.parseInt(prns);  // get the PRN
        return prn;
    }
    
    
    private GPSTimeFormat parseTOC(String str){
        
        // grab just the first part of the 1st line before the clock bias
        String string = str.substring(0, 21);
        
        StringTokenizer tok = new StringTokenizer(string, " ");
        String prn = tok.nextToken();
        
        String yrs = tok.nextToken();
        String mons = tok.nextToken();
        String days = tok.nextToken();
        String hrs = tok.nextToken();
        String mins = tok.nextToken();
        String secs = tok.nextToken();
        int year = Integer.parseInt(yrs);
        //* The following line is a "fix" for the Y2K bug
        //* The code will break again after 2079
        if(year>80) year = year + 1900;
        else year = year + 2000;
        int month = Integer.parseInt(mons);
        int day = Integer.parseInt(days);
        int hour = Integer.parseInt(hrs);
        int min = Integer.parseInt(mins);
        double sec = Double.parseDouble(secs);
        CalDate toc_dt = new CalDate(year, month, day, hour, min, sec);
        GPSTimeFormat toc = new GPSTimeFormat(toc_dt);
        return toc;
    }
    
    public ArrayList process() {
        ArrayList out = new ArrayList();
        try {
            
            System.out.println("Processing File: "+filename);
            FileReader fr = new FileReader(filename);
            BufferedReader in = new BufferedReader(fr);
            String s;
            
            // read first line and get version
//            System.out.println("Reading first line");
            s = in.readLine();
            StringTokenizer t1 = new StringTokenizer(s, " ");
            double version = Double.parseDouble(t1.nextToken());
            if (version < 2.0){
                System.out.println("RINEX Version < 2.0. I can't handle this.");
                System.exit(99);
            }
            
//            System.out.println("RINEX Version "+version);
            
            // read second and third lines
            s = in.readLine();
            //s = in.readLine();
            
            // read fourth line
            try{
            //s = in.readLine();
            //t1 = new StringTokenizer(s, " ");
            //for (int i = 0; i < 4; i++) {
            	//this.alpha[i] = this.parse(t1.nextToken());
            //}
            
            // read fifth line
            //s = in.readLine();
            //t1 = new StringTokenizer(s, " ");
            //for (int i = 0; i < 4; i++) {
            	//this.beta[i] = this.parse(t1.nextToken());
            //}
            
            // find the end of header
            String end = "END";
            String header = "HEADER";
            int check = 0;
            
            boolean eoh = false;
            while (!eoh) {
                s = in.readLine();
//                System.out.println(s);
                StringTokenizer t2 = new StringTokenizer(s, " ");
                while (t2.hasMoreTokens()) {
                    String temp = t2.nextToken();
                    if ((temp.equals(end))||(temp.equals(header))) check++;
                    if (check > 1) eoh = true;
                }
            }
            }catch(Exception e){
            	e.printStackTrace();
            	 String end = "END";
                 String header = "HEADER";
                 int check = 0;
                 
                 boolean eoh = false;
                 while (!eoh) {                 
//                     System.out.println(s);
                     StringTokenizer t2 = new StringTokenizer(s, " ");
                     while (t2.hasMoreTokens()) {
                         String temp = t2.nextToken();
                         if ((temp.equals(end))||(temp.equals(header))) check++;
                         if (check > 1) eoh = true;
                     }
                     s = in.readLine();
                 }
            }
            
//            System.out.println("got to end of header");
            
            // read in all the data records
            
            boolean eof = false;
            int numLines = 8;
            
            String [] line = new String[numLines];
            
            while (!eof){
                
                // read in a group of 8 lines per SV
                
                boolean completeRecord = true;
                
                for (int i = 0; i < numLines; i++){
                    if ((line[i] = in.readLine()) == null){
                        eof = true;
                        completeRecord = false;
                        //                 System.out.println("Incomplete record");
                    }
                    //             else {
                    //                 System.out.println("line "+i+": "+line[i]);
                    //             }
                }
                
                //process the record
                if ((!eof)&&(completeRecord)){
                    int prn = parsePRN(line[0]);  // get the PRN
                    Integer PRN = new Integer(prn);
//                    System.out.println("Processing PRN "+prn);
                    boolean duplicate = PRNlist.contains(PRN);
//                    if (duplicate){
//                        System.out.println("Already processed PRN "+prn+". Skipping this record.");
//                    }
                    
                    // check for complete record and redundant PRN
                    if (!duplicate) {
                        
                        // grab data from first line
                        GPSTimeFormat toc = parseTOC(line[0]);
                        double [] temp1 = parseFirstLine(line[0]);
                        
                        double bias = temp1[0];
                        double drift = temp1[1];
                        double driftrate = temp1[2];
                        
                        // second line
                        double [] temp = parseLine(line[1]);
                        double iode = temp[0];
                        double crs = temp[1];
                        double deltan = temp[2];
                        double m0 = temp[3];
                        
                        //third line
                        temp = parseLine(line[2]);
                        double cuc = temp[0];
                        double ecc = temp[1];
                        double cus = temp[2];
                        double sqrta = temp[3];
                        
                        //fourth line
                        temp = parseLine(line[3]);
                        double toe_sow = temp[0];
                        double cic = temp[1];
                        double omega = temp[2];
                        double cis = temp[3];
                        
                        //fifth line
                        temp = parseLine(line[4]);
                        double inc = temp[0];
                        double crc = temp[1];
                        double w = temp[2];
                        double omegadot = temp[3];
                        
                        //sixth line
                        temp = parseLine(line[5]);
                        double idot = temp[0];
                        double codes = temp[1];
                        double gpsweek = temp[2];
                        long toe_gpsweek = (long)gpsweek;                       // truncate to get a long
                        
                        //seventh line
                        temp = parseLine(line[6]);
                        double accuracy = temp[0];
                        double health = temp[1];
                        double tgd = temp[2];
                        double iodc = temp[3];
                        
                        boolean goodhealth = true;
                        if (health != 0.0) {
//                        	System.out.println("Bad Health: PRN: "+prn+" health = "+health);
                        	goodhealth = false;
                        }
                                                
                        // create a new SV and initialize it
                        if (goodhealth) {
                        	
                        	PRNlist.add(PRN);                        // add this PRN to the list
                        	                        	
	                        GPS_SV sv = new GPS_SV(prn);            // create a new SV
	                        sv.setTOC(toc);
	                        sv.setClockParams(bias, drift, driftrate);
	                        sv.setHarmonicCorrections(crc, crs, cuc, cus, cic, cis);
	                        sv.setOrbit(sqrta, ecc, inc, omega, w, m0);
	                        sv.setOrbitCorrections( deltan, omegadot, idot);
	                        GPSTimeFormat toe = new GPSTimeFormat(toe_gpsweek, toe_sow);
	                        sv.setTOE(toe);
	                        
	                        // add the new SV to the SVList
	                        out.add(sv);
                        }
                    }
                }
            }
            in.close();
        }
        catch (IOException e)                 // what to do if an i/o error occurs
        {
            System.out.println("Error: "+e);
            System.exit(1);
        }
        System.out.println("RINEX file processing completed");
        System.out.println("Total number of SVs = "+out.size());
        System.out.println("==================================");
        return out;
    }
    
    public double[] alpha(){
    	return this.alpha;
    }
    
    public double[] beta(){
    	return this.beta;
    }
        
}
