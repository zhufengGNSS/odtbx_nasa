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
package jat.forces.harrispriester;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import jat.util.FileUtil;

/**
 * @author Richard C. Page III
 *
 */
public class HarrisPriesterData {

    private String filename;
    private int[] line_ref;
    private int[] data_n;
    private double[] f107;
    
    protected double[][] data1;
    protected double[][] data2;
    protected double[][] data;
    protected int ref = 4;
    protected double ref_f107=150;
    protected boolean interp = false;
   
    public HarrisPriesterData(){
        String fs = FileUtil.file_separator();
		String directory = 
		    FileUtil.getClassFilePath("jat.forces.harrispriester","HarrisPriesterData")+fs;
		filename = directory+"harrispriester.txt";
		f107 = new double[10];
		line_ref = new int[10];
		data_n = new int[10];
    }
    
    public void process(double arg_f107){
        this.ref_f107 = arg_f107;
        FileReader fr;
		BufferedReader in = null;
		try {
			fr = new FileReader(filename);
			in = new BufferedReader(fr);			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
//		ArrayList lineCollection = new ArrayList();

		// read in the file
		boolean eof = false;
		String line = "";
		try {   
		    while(line.indexOf("****")<0){
		        line = in.readLine();
		    }
		    line = in.readLine();
		    line = in.readLine();
		    for(int i=0; i<10; i++){
		        line = in.readLine();
		        line = in.readLine();
		        this.data_n[i] = Integer.parseInt(line.substring(36));
		        line = in.readLine();
		        this.f107[i] = Double.parseDouble(line.substring(67));
		        line = in.readLine();
		        this.line_ref[i] = Integer.parseInt(line.substring(7));
		        line = in.readLine();
		    }
		    for(int i=0; i<10; i++){
		        if(arg_f107 >= this.f107[i]) ref = i;
		        else break;
		    }
		    this.data = new double[this.data_n[ref]][3];
		    in.close();
			try {
				fr = new FileReader(filename);
				in = new BufferedReader(fr);			
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		    for(int i=0; i<this.line_ref[ref]; i++) line = in.readLine();
		    data1 = new double[data_n[ref]][3];
		    for(int i=0; i<this.data_n[ref]; i++){
		        line = in.readLine();
		        this.data1[i][0] = Double.parseDouble(line.substring(7,20));
		        this.data1[i][1] = Double.parseDouble(line.substring(29,42));
		        this.data1[i][2] = Double.parseDouble(line.substring(51,64));
		    }
		    if(arg_f107>this.f107[ref]){
		        interp = true;
		        data2 = new double[data_n[ref+1]][3];
		        int tmp = this.data_n[ref]+this.line_ref[ref];
		        for(int i=tmp; i<this.line_ref[ref+1]; i++) in.readLine();
			    for(int i=0; i<this.data_n[ref]; i++){
			        line = in.readLine();
			        this.data2[i][0] = Double.parseDouble(line.substring(7,20));
			        this.data2[i][1] = Double.parseDouble(line.substring(29,42));
			        this.data2[i][2] = Double.parseDouble(line.substring(51,64));
			    }
			    data = data2;
		    }else{
		        data = data1;
		    }
		    in.close();
		} catch (IOException e1) {
		    e1.printStackTrace();
		}
		//if(interp) interpolate(arg_f107);
    }
    
    protected void interpolate(double arg_f107){
        for(int i=0; i<data_n[ref]; i++){
            //*assume altitude is parallel
            data[i][0] = data1[i][0];
            data[i][1] = data1[i][1]+(data2[i][1]-data1[i][1])*(arg_f107-f107[ref])/
            						 (f107[ref+1]-f107[ref]);
            data[i][2] = data1[i][2]+(data2[i][2]-data1[i][2])*(arg_f107-f107[ref])/
			 						(f107[ref+1]-f107[ref]);
        }
    }

	public int get_data_n(){ return this.data_n[ref];}
	public double[][] get_data(){ return this.data;}
    
    public static void main(String[] args) {
        HarrisPriesterData d = new HarrisPriesterData();
        d.process(150);
        double[][] data = d.get_data();
        int n = d.get_data_n();
        for(int i=0; i<n; i++){
            System.out.println(""+data[i][0]+"   "+data[i][1]+"   "+data[i][2]);
        }
    }
}
