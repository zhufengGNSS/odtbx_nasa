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
 */

package jat.sim;



import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;


public class initializer{

    String generic_prefix;
    static HashMap hm = new HashMap();
    
    public initializer(){
    	generic_prefix = new String("FFTB.");	
    }
    
    public HashMap configure_components(String hardfile, String scenfile){
    	HashMap temp;
    
    	temp = parse_file(hardfile);
    	hm.putAll(temp);
    	
    	temp = parse_file(scenfile);
    	hm.putAll(temp);
    	
    	return hm;
    }
    
    public HashMap get_config(){
    	return hm;
    }
    
    public String gen_truth_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.STARS.STATE"));
    }

    public String gen_control_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.CTRL.DELTA-V"));
    }

    public String gen_deltav_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.FE.DELTA-V"));
    }

    public String gen_spirent_control_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.SPIRENT.CTRL"));
    }

    public String gen_f14_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.GPSA.F14"));
    }

    public String gen_f40_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.GPSA.F40"));
    }

    public String gen_f42_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.GPSA.F42"));
    }

    public String gen_f44_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.GPSA.F44"));
    }

    public String gen_f45_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.GPSA.F45"));
    }

    public String gen_f46_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.GPSA.F46"));
    }

    public String gen_f62_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.GPSA.F62"));
    }

    public String gen_xlink_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.CROSSLINK.2WAY"));
    }

    public String gen_est_state_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.FE.STATE"));
    }

    public String gen_rec_cmd_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.DIR.ORION"));
    }

    public String gen_rtmon_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.FE.RTMON"));
    }

    public String gen_timer_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.TIMER.TICK"));
    }

    public String gen_count_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.TIMER.COUNT"));
    }

    public String gen_ephem_request_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.GPSA.EPHEMREQ"));
    }

    public String gen_error_control_grp(String mission, Integer sat){
    	return formulate_name(mission, sat, generic_prefix, new String("MSG.RTSD.ERROR.CONTROL"));
    }

    private static String formulate_name(String mission, Integer sat, String prefix, String postfix){
    	StringBuffer temp = new StringBuffer();
    	temp.append(prefix);
    	temp.append(mission);
    	temp.append(".");
    	temp.append(sat.intValue());
    	temp.append(".");
    	temp.append(postfix);
    	return temp.toString();
    }
    
    // nowhere near as good as Jason's, but should work. -DMZ
    public static HashMap parse_file(String filename){
    	BufferedReader br = null;
    	HashMap hm = new HashMap();
    	String line = null, key = new String(""), value = new String(""), temp=null;
    	StringTokenizer st = null;
    	boolean first=true;
    	
    	//open the file
    	try {
			br = new BufferedReader(new FileReader(filename));
		} catch (FileNotFoundException e) {
			System.out.println("Sorry - could not find the file "+filename+".");
			System.exit(1);
			e.printStackTrace();
		}

		try {
			while((line = br.readLine()) != null){

				value=new String("");
				first=true;
				
			    // ignore comments
			    if(line.startsWith("#")){
			    	continue;
			    }

				// tokenize the line (regex might be better?)
				st = new StringTokenizer(line);
				
				// blank line
				if(!st.hasMoreTokens()){
					continue;
				}
				
				key = st.nextToken();
				
				while(st.hasMoreTokens()){
					temp=st.nextToken();
					if(temp.substring(0,0).equals("#")){
						break;
					}
					else{
						if(first){
							value=value.concat(temp);
							first=false;
						}
						else{
							value=value.concat(" ");
							value=value.concat(temp);
						}

					}

				}
				
				// strip out leading and trailing spaces in key and value
				key=key.trim();
				value=value.trim();
				
				System.out.println("Key is: \""+key+"\" Value is: \"" + value +"\"");
				
				// stash the key and value into the HashMap
				hm.put(key, value);
			
			}
		} catch (IOException e) {
			System.out.println("Sorry - caught an IOException while reading file.");
			e.printStackTrace();
			System.exit(1);
		}
		
		return hm;
    }
    
    public static Double parseDouble (HashMap hm, String key)
    {
    	Double out = null;
    	
   	 	Object ob =  hm.get(key);
        if (ob != null) {
          String ss = ob.toString();
          out = new Double(Double.parseDouble(ss));
        }
        
    	return out;
    }
    
    public static Integer parseInt (HashMap hm, String key)
    {
    	Integer out = null;
    	
   	 	Object ob =  hm.get(key);
        if (ob != null) {
          String ss = ob.toString();
          out = new Integer(Integer.parseInt(ss));
        }
    	return out;
    }
    
     
    public static Boolean parseBool (HashMap hm, String key)
    {
    	Boolean out = null;
    	
   	 	Object ob =  hm.get(key);
        if (ob != null) {
          String ss = ob.toString();
          if(ss.equals("true")||ss.equals("True"))
            out = Boolean.TRUE;
          else
   	 		out = Boolean.FALSE;
        }
    	
    	return out;
    }
    
    public static String parseString(HashMap hm, String key)
    {
        String out = null;
        
   	 	Object ob =  hm.get(key);
        if (ob != null) {
          out = ob.toString();
        }
        
    	return out;
    }
}
