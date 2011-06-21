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
 * 
 * File Created on 3 October 2005
 * Created by Kathryn Bradley, Emergent Space Technologies
 * */
package jat.measurements;

import jat.gps.GPS_Measurement;
import jat.matvec.data.Matrix;
import jat.matvec.data.RotationMatrix;
import jat.matvec.data.VectorN;
import jat.sim.initializer;
import jat.spacetime.CalDate;
import jat.spacetime.EarthRef;
import jat.spacetime.FitIERS;
import jat.spacetime.GPSTimeFormat;
import jat.spacetime.Time;
import jat.traj.Trajectory;
import jat.util.FileUtil;

import java.io.BufferedReader;
import java.io.EOFException;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.Vector;

public class ObservationMeasurementList 
{
	int count=0;
	int numOfSatellites=0;
	int timeStep = 0;
	int currentIndex = 0;
	int observationNumber=0;
	int dataCount=0;
	
	double version;
	double timeMjd;
	double current_mjd;
	
	boolean notDone = true;
	boolean timeStamp = false;
	
	String line;
	String lineNew;
	
	Vector data = new Vector();
	Vector observationType= new Vector();
	Vector prns = new Vector();
	ArrayList list = new ArrayList();
	
	GPSmeasurementModel model_gps;
	GPSStateMeasurementModel model_gpsstate;
	rangeMeasurementModel model_range;
	stateMeasurementModel model_state;
	stateUpdateMeasurementModel model_stateupdate;
	
	private HashMap input;
	
	public ObservationMeasurementList(){}
	
	public ObservationMeasurementList(HashMap hm){
		this.input = hm;
		for(int i=0; i<initializer.parseInt(input, "MEAS.types"); i++){
			if(initializer.parseString(input, "MEAS."+i+".desc").equalsIgnoreCase("GPS")){
				model_gps = new GPSmeasurementModel(hm);
			}else if(initializer.parseString(input, "MEAS."+i+".desc").equalsIgnoreCase("pseudoGPS")){
				model_gpsstate = new GPSStateMeasurementModel(hm);
			}else if(initializer.parseString(input, "MEAS."+i+".desc").equalsIgnoreCase("range")){
				model_range = new rangeMeasurementModel(hm);
			}else if(initializer.parseString(input, "MEAS."+i+".desc").equalsIgnoreCase("stateUpdate")){
				model_state = new stateMeasurementModel(hm);
			}else{
				model_stateupdate = new stateUpdateMeasurementModel(hm);
			}
		}
	}
	
	public static void main(String[] args) // Main Method
	throws IOException{
		ObservationMeasurementList x = new ObservationMeasurementList();
		String path = FileUtil.getClassFilePath("jat.measurements","ObservationMeasurementList");
		String fs = FileUtil.file_separator();
		//x.processRINEX(path+fs+"ExampleRINEXGEONS.rnx");
		//x.processRINEX(path+fs+"Case-820.rnx");
		x.processStateUpdateFile(path+fs+"test1_8.rnx",50985,50999);
	}
	
	public void processStateUpdateFile(String fileName, double mjd0, double mjdf) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(new File(fileName)));
		boolean loop = true;
		Time time;
		EarthRef earth;
		RotationMatrix rot;
		VectorN r = new VectorN(3);
		int year, day,i;
		double sec, mjd;
		GPSTimeFormat date;
		double x,y,z,clock;
		String line = "start";
		StringTokenizer tok;
		boolean processedOther = false;
		ObservationMeasurement om;
		Vector v_type = new Vector();
		v_type.add("S3");
		Vector v_data;
		if(this.list.size()>0) processedOther = true;
		line = in.readLine();
		while(!line.equals("") && loop){
			try{
				tok = new StringTokenizer(line, " ");
				year = Integer.parseInt(tok.nextToken());
				day = Integer.parseInt(tok.nextToken());
				sec = Double.parseDouble(tok.nextToken());
				date = new GPSTimeFormat(year,day,0,0,sec);
				mjd = date.mjd();
				if(mjd>=mjd0 && mjd<= mjdf){
					x = Double.parseDouble(tok.nextToken());
					y = Double.parseDouble(tok.nextToken());
					z = Double.parseDouble(tok.nextToken());
					clock = Double.parseDouble(tok.nextToken())/jat.cm.Constants.c;
					mjd = mjd - (clock)/86400;
					r = new VectorN(x,y,z);
					time = new Time(mjd);
					//*TODO watch
					earth = new EarthRef(time);
					
					FitIERS iers = new FitIERS();
					//iers.process(); //* don't need this
					double[] param = iers.search(time.mjd_tt());
					earth.setIERS(param[0],param[1]);
					time.set_UT1_UTC(param[2]);
					//earth.computePole(time);
					
					Matrix tmp = earth.eci2ecef(time);
					rot = new RotationMatrix(tmp.transpose());
					r = rot.transform(r);
					v_data = new Vector();
					v_data.add(""+r.x[0]);
					v_data.add(""+r.x[1]);
					v_data.add(""+r.x[2]);
					if(processedOther){
						i = searchListMJD(mjd);
						for(int snum=0; snum<3; snum++){
							om = new ObservationMeasurement(mjd,v_type,v_data,"0",this.input,
									this.model_gpsstate,ObservationMeasurement.TYPE_GPSSTATE);
							om.set_whichState(snum);
							System.out.println(om.toString());
							list.add(i,om);
						}
					}else{
						for(int snum=0; snum<3; snum++){
							om = new ObservationMeasurement(mjd,v_type,v_data,"0",this.input,
									this.model_gpsstate,ObservationMeasurement.TYPE_GPSSTATE);
							om.set_whichState(snum);
							System.out.println(om.toString());
							list.add(om);
						}
					}
				}
				line = in.readLine();
			}catch(EOFException eof){
				loop = false;
			}
		}
	}
	/**		Method Class which accesses the other classes to read the RINEX observation files.*/
	public void processRINEX(String fileName)	//Needs the fileName 
	throws IOException{
		File inputFile = new File(fileName);
		FileReader fr = new FileReader(inputFile);
		BufferedReader in = new BufferedReader(fr);
		
		while (notDone) //Reading the Header information in this While loop
		{
			count++;
			line = in.readLine();
			StringTokenizer tok = new StringTokenizer(line, " ");
			if (count == 1)
			{
				version = Double.parseDouble(tok.nextToken());
			}
			else if (count == 2)
			{
				observationNumber = Integer.parseInt(tok.nextToken());
				for (int k =0; k < observationNumber; k++)
				{
					observationType.add(tok.nextToken());
				}
			}
			else 
			{
				ObservationMeasurementList headerFinished = new ObservationMeasurementList();
				notDone = headerFinished.setHeaderReader(tok, count); //Returns false if header has been read
			}
		}
		while (in.ready()) //Reading the file blocks with the measurement data
		{
			/**dataCount is invented for first line of data. The timeStamp 
			 is true if the data from the previous time has already been read.
			 Initially timeStamp is set to false and dataCount is set to zero*/
			
			dataCount++;
			if ((dataCount==1) || (timeStamp == true)) 
			{
				lineNew = in.readLine();
				//ObservationMeasurementList date = new ObservationMeasurementList();
				GPSTimeFormat newDate =  ObservationMeasurementList.setDateStuff(lineNew);
				timeMjd =newDate.mjd();
				//if(timeMjd > 52189.06403935186){
					//int stop_to_debug = 0;
				//}
				if(currentIndex <0){
					current_mjd = timeMjd;
					currentIndex = 0;
				}
				//ObservationMeasurementList num = new ObservationMeasurementList();
				numOfSatellites = ObservationMeasurementList.setNumberOfSatellites(lineNew);
				//data.add(new Integer(numOfSatellites));
				int newNumOfSatellites = numOfSatellites;
				int start = 32;
				if (numOfSatellites >12)
				{
					newNumOfSatellites =12;
				}
				//prns.clear();
				prns = new Vector();
				for (int t=0; t<newNumOfSatellites; t++)
				{
					String prnIds = new String(lineNew.substring(start, start+3));
					prns.add(prnIds);
					start = start + 3;
				}
				if (numOfSatellites > 12){
					String secondLine = in.readLine();
					
					int secondLineSatellites = numOfSatellites-12;
					int secondStart = 33;
					for (int n=0; n<secondLineSatellites; n++)
					{
						String prnIds = new String (secondLine.substring(secondStart, secondStart+2));
						prns.add(prnIds); 
						secondStart = secondStart +3;
					}
				}
				timeStamp = false; //Resets timeStamp to read the data for the time stamp
			}
			else 
			{
				for (int j=0; j<numOfSatellites; j++ )
				{
					String lineNew = in.readLine();
					StringTokenizer tokNew = new StringTokenizer(lineNew, " ");
					//data.clear();
					data = new Vector();
					for (int t=0; t<observationNumber; t++)
					{
						String dataNew = tokNew.nextToken();
						data.add(dataNew);					
					}
					ObservationMeasurement info;
					switch(ObservationMeasurement.chooseType((String)prns.get(j))){
					case ObservationMeasurement.TYPE_GPS:
						info = new ObservationMeasurement(timeMjd,observationType, data, prns.get(j),
								input,model_gps,ObservationMeasurement.TYPE_GPS);
						EarthRef earth = new EarthRef(new Time(timeMjd));
						RotationMatrix R = new RotationMatrix(earth.ECI2ECEF());
						info.set_ECF2ECI(R);
						break;
					case ObservationMeasurement.TYPE_STATEUPDATE:
						info = new ObservationMeasurement(timeMjd,observationType, data,prns.get(j),
								input,model_stateupdate,ObservationMeasurement.TYPE_STATEUPDATE);
						break;
					case ObservationMeasurement.TYPE_RANGE:
						info = new ObservationMeasurement(timeMjd,observationType, data, prns.get(j),
								input,model_range, ObservationMeasurement.TYPE_RANGE);
						break;
//* TODO Include the rest of the model types	
					default:
						info = new ObservationMeasurement(
								timeMjd,observationType, data, prns.get(j),input,model_gps,
								ObservationMeasurement.TYPE_GPS);
					}
					this.add(info);
					String meas = info.toString();
					System.out.println(meas);
					
				}
				timeStamp = true; //Resets timeStamp to read the next time stamp line
				timeStep++; //keeps track of the number of time steps, initially zero
			}
		}
	}
	public void processRINEX(String fileName,boolean useGPS, boolean useCross)	//Needs the fileName 
	throws IOException{
		File inputFile = new File(fileName);
		FileReader fr = new FileReader(inputFile);
		BufferedReader in = new BufferedReader(fr);
		double MJDF = initializer.parseDouble(input, "init.MJDF")+initializer.parseDouble(input, "init.TF")/86400.0;
		
		while (notDone) //Reading the Header information in this While loop
		{
			count++;
			line = in.readLine();
			StringTokenizer tok = new StringTokenizer(line, " ");
			try{
			if (count == 1)
			{
				version = Double.parseDouble(tok.nextToken());
			}
			else if (count == 2)
			{
				observationNumber = Integer.parseInt(tok.nextToken());
				for (int k =0; k < observationNumber; k++)
				{
					observationType.add(tok.nextToken());
				}
			}
			else 
			{
				ObservationMeasurementList headerFinished = new ObservationMeasurementList();
				notDone = headerFinished.setHeaderReader(tok, count); //Returns false if header has been read
			}
			}catch(Exception e){
				ObservationMeasurementList headerFinished = new ObservationMeasurementList();
				notDone = headerFinished.setHeaderReader(tok, count); //Returns false if header has been read
			}
		}
		while (in.ready()) //Reading the file blocks with the measurement data
		{
			/**dataCount is invented for first line of data. The timeStamp 
			 is true if the data from the previous time has already been read.
			 Initially timeStamp is set to false and dataCount is set to zero*/
			
			dataCount++;
			if ((dataCount==1) || (timeStamp == true)) 
			{
				lineNew = in.readLine();
				//ObservationMeasurementList date = new ObservationMeasurementList();
				GPSTimeFormat newDate =  ObservationMeasurementList.setDateStuff(lineNew);
				//timeMjd =newDate.mjd();
				//* TODO Watch this
				timeMjd =newDate.mjd_utc();//-13.0/86400.0;
				if(timeMjd > MJDF) 
					break;
				//if(timeMjd > 52189.06403935186){
					//int stop_to_debug = 0;
				//}
				if(currentIndex <0){
					current_mjd = timeMjd;
					currentIndex = 0;
				}
				//ObservationMeasurementList num = new ObservationMeasurementList();
				numOfSatellites = ObservationMeasurementList.setNumberOfSatellites(lineNew);
				//data.add(new Integer(numOfSatellites));
				int newNumOfSatellites = numOfSatellites;
				int start = 32;
				if (numOfSatellites >12)
				{
					newNumOfSatellites =12;
				}
				//prns.clear();
				prns = new Vector();
				for (int t=0; t<newNumOfSatellites; t++)
				{
					String prnIds = new String(lineNew.substring(start, start+3));
					prns.add(prnIds);
					start = start + 3;
				}
				if (numOfSatellites > 12){
					String secondLine = in.readLine();
					
					int secondLineSatellites = numOfSatellites-12;
					int secondStart = 33;
					for (int n=0; n<secondLineSatellites; n++)
					{
						String prnIds = new String (secondLine.substring(secondStart, secondStart+2));
						prns.add(prnIds); 
						secondStart = secondStart +3;
					}
				}
				timeStamp = false; //Resets timeStamp to read the data for the time stamp
			}
			else 
			{
				for (int j=0; j<numOfSatellites; j++ )
				{
					String lineNew = in.readLine();
					StringTokenizer tokNew = new StringTokenizer(lineNew, " ");
					//data.clear();
					data = new Vector();
					for (int t=0; t<observationNumber; t++)
					{
						String dataNew = tokNew.nextToken();
						data.add(dataNew);					
					}
					ObservationMeasurement info;
					String meas;
					switch(ObservationMeasurement.chooseType((String)prns.get(j))){
					case ObservationMeasurement.TYPE_GPS:
						if(useGPS){
							info = new ObservationMeasurement(timeMjd,observationType, data, prns.get(j),
								input,model_gps,ObservationMeasurement.TYPE_GPS);
							EarthRef earth = new EarthRef(new Time(timeMjd));
							RotationMatrix R = new RotationMatrix(earth.ECI2ECEF());
							info.set_ECF2ECI(R);
							if(!(info.prn.equals("G17") || info.prn.equals("G24") || info.prn.equals("G30"))){
								this.add(info);
							}							
							meas = info.toString();
							System.out.println(meas);
						}
						break;
					case ObservationMeasurement.TYPE_STATEUPDATE:
						info = new ObservationMeasurement(timeMjd,observationType, data,prns.get(j),
								input,model_stateupdate,ObservationMeasurement.TYPE_STATEUPDATE);
						this.add(info);
						meas = info.toString();
						System.out.println(meas);
						break;
					case ObservationMeasurement.TYPE_RANGE:
						if(useCross){
							info = new ObservationMeasurement(timeMjd,observationType, data, prns.get(j),
									input,model_range, ObservationMeasurement.TYPE_RANGE);
							this.add(info);
							meas = info.toString();
							System.out.println(meas);
						}
						break;
//* TODO Include the rest of the model types	
					default:
						info = new ObservationMeasurement(
								timeMjd,observationType, data, prns.get(j),input,model_gps,
								ObservationMeasurement.TYPE_GPS);
						this.add(info);
						meas = info.toString();
						System.out.println(meas);
					}
					
				}
				timeStamp = true; //Resets timeStamp to read the next time stamp line
				timeStep++; //keeps track of the number of time steps, initially zero
			}
		}
	}
	public void assignModels(){
		for(int i=0; i<this.list.size(); i++){
			ObservationMeasurement o = get(i);
			
			
		}
	}
	
	/** Add an Observation Measurement to the collection
	 * @param meas ObservationMeasurementList object
	 */
	public void add(ObservationMeasurement meas) {
		list.add(meas);
	}
	
	/** Get an ObservationMeasurementList out of the collection
	 * @param index Index of the measurement
	 * @return the ObservationMeasurementList
	 */
	public ObservationMeasurement get(int index) {
		return (ObservationMeasurement) list.get(index);
	}
	public ObservationMeasurement getCurrent(){
		ObservationMeasurement om = (ObservationMeasurement)list.get(currentIndex);
		return om;
	}
	public ObservationMeasurement getNext(){
		if(currentIndex>=(list.size()-1)){
			//currentIndex=(list.size()-1);
			return null;
		}else
			currentIndex++;
		ObservationMeasurement om = (ObservationMeasurement)list.get(currentIndex);
		return om;
	}
	public ObservationMeasurement getPrev(){
		if(currentIndex>0)
			currentIndex--;
		else
			currentIndex=0;
		ObservationMeasurement om = (ObservationMeasurement)list.get(currentIndex);
		return om;
	}
	public ObservationMeasurement getFirst(){
		currentIndex = 0;
		return (ObservationMeasurement)list.get(currentIndex);
	}
	public int getIndex(){
		return currentIndex;
	}
	public Trajectory get_meas_traj(){
		Trajectory traj = new Trajectory();
		double t;
		double[] data = new double[6];
		ObservationMeasurement tmp;
		for(int i=0; i<this.list.size(); i++){
			tmp = (ObservationMeasurement)list.get(i);
			if(tmp.get_type()==ObservationMeasurement.TYPE_GPSSTATE){
				t = tmp.get_mjd();
				for(int j=0; j<3; j++)
					data[j] = tmp.get_state(3).x[j];
				for(int j=3; j<6; j++)
					data[j] = 0.0;
				traj.add(t,data);
			}
		}
		return traj;
	}
	
	/**
	 * Search the list and return the first index of the element which is later than or equal to mjd.  
	 * @param mjd Modified Julian Date to search for
	 * @return Index
	 */
	public int searchListMJD(double mjd){
		int i;
		for(i=0; i<list.size(); i++){
			if(((ObservationMeasurement)list.get(i)).time_mjd()>=mjd)
				return i;
		}
		return i;
	}
	
	/** This method setNumberOfSatellites gets the number of satellites from the
	 * file and returns it to the process method. */
	private static int setNumberOfSatellites(String lineNew)
	{
		int num;
		StringTokenizer satelliteToken = new StringTokenizer(lineNew.substring(29,32), " ");
		String numDum = satelliteToken.nextToken();
		num = Integer.parseInt(numDum);
		return num;
	}
	
	/** This method setDateStuff reads the time stamp from the line being fed 
	 * in and then calls two methods GPSTimeFormat and CalDate to determine the 
	 * time in GPS time. This time is then sent to the Process method to be 
	 * stored. */
	private static GPSTimeFormat setDateStuff(String lineNew)
	{
		StringTokenizer newTok = new StringTokenizer(lineNew, " ");
		String yrsString = newTok.nextToken();
		String monString = newTok.nextToken();
		String dayString = newTok.nextToken();
		String hrString  = newTok.nextToken();
		String minString = newTok.nextToken();
		String secString = newTok.nextToken();
		int month = Integer.parseInt(monString);
		int day = Integer.parseInt(dayString);
		int hour  =Integer.parseInt(hrString);
		int min = Integer.parseInt(minString);
		double sec = Double.parseDouble(secString);
		int years = Integer.parseInt(yrsString);
		int year = 0;
		if (years < 50)
		{
			year = (years + 2000);
		}
		else 
		{
			year = (years + 1900);
		}
		CalDate toc_dt = new CalDate(year, month, day, hour, min, sec);
		GPSTimeFormat date = new GPSTimeFormat(toc_dt); 
		return date;
	}	
	/** This method setHeaderRead determines if the header has finished being 
	 * read. If the header is finished being read then it returns a false 
	 * to the Process method. */
	public boolean setHeaderReader(StringTokenizer tok, int count)
	{
		String end = "END";
		String header = "HEADER";
		String of = "OF";
		int check=0;
		boolean headerFinished = true;
		while(tok.hasMoreTokens())
		{
			String temp = tok.nextToken();
			if ((temp.equals(end))||(temp.equals(of))||(temp.equals(header))) check++;
			{
				if (check > 2) 
				{
					headerFinished = false; 
				}
			}	
		}
		return headerFinished;
	}
}

