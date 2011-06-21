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
package jat.spacecraft;

import java.util.*;

import jat.matvec.data.Matrix;
import jat.matvec.data.VectorN;

/**
 * This class represents a spacecraft formation.  It contains any number of spacecraft
 * with one designated as the primary spacecraft.
 * 
 * @author Richard C. Page III
 *
 */
public class Formation {

    /**
     * String ID for the formation
     */
    public String id;
    /**
     * List of member spacecraft
     */
    protected ArrayList member_spacecraft;
    /**
     * Parallel list of member spacecraft IDs
     */
    protected ArrayList member_id;
    /**
     * Primary spacecraft of the formation.
     */
    protected SpacecraftModel primary;
    /**
     * Connectivity table.  (Graph) Spacecraft represent nodes and their 
     * relative distances represent edges retrieved by calls to the member
     * spacecraft.
     */
    protected Matrix connect;
    
    /**
     * Default Constructor.
     */
    public Formation(){
        connect = new Matrix(2);
        member_spacecraft = new ArrayList();
        member_id = new ArrayList();
    }
    /**
     * Constructor - initializes for a certain number of spacecraft
     * @param num Number of spacecraft.
     */
    public Formation(int num){
        connect = new Matrix(num);
        member_spacecraft = new ArrayList(num);
        member_id = new ArrayList(num);
    }
    /**
     * Constructor - initializes the primary spacecraft.
     * @param p Primary spacecraft.
     */
    public Formation(SpacecraftModel p){
        primary = p;
        connect = new Matrix(2);
        member_spacecraft = new ArrayList();
        member_id = new ArrayList();
    }
    /**
     * Add a spacecraft to the formation.  If no primary has been designated,
     * the first spacecraft added will become the primary.
     * @param s Spacecraft model.
     */
    public void add_spacecraft(SpacecraftModel s){
        if(primary == null){
            primary = s;
        } else{            
            member_spacecraft.add(s);
            member_id.add(s.get_id());
        }
    }
    /**
     * Add a spacecraft to the formation.  If no primary has been designated,
     * the first spacecraft added will become the primary.
     * @param s Spacecraft used to create a generic spacecraft model.
     */
    public void add_spacecraft(Spacecraft s){
        if(primary == null){
            primary = new SpacecraftModel(s);
        } else {
            SpacecraftModel sc = new SpacecraftModel(s);
            member_spacecraft.add(sc);
            member_id.add(s.id);
        }
    }
    /**
     * Return the spacecraft model for a given index (in order of addition).
     * @param i Index
     * @return Spacecraft model
     */
    public SpacecraftModel get_spacecraftmodel(int i){
        if(i>0)
            return (SpacecraftModel)member_spacecraft.get(i-1);
        else if(i==0)
            return primary;
        else
            return null;
    }
    /**
     * Return the spacecraft object for a given index (in order of addition).
     * @param i Index
     * @return Spacecraft
     */
    public Spacecraft get_spacecraft(int i){
        if(i>0){
            SpacecraftModel sc = (SpacecraftModel)member_spacecraft.get(i-1);
            return sc.get_spacecraft();
        } else if(i==0)
            return primary.get_spacecraft();
        else
            return null;
    }
    /**
     * Search the list for a spacecraft corresponding to the given ID. 
     * @param id
     * @return Spacecraft Model or null if not found
     */
    public SpacecraftModel get_spacecraftmodel(String id){
        int i = search(id);
        if(i>0)
            return (SpacecraftModel)member_spacecraft.get(i-1);
        else if(id.equalsIgnoreCase(primary.get_id()))
            return primary;
        else
            return null;
    }
    /**
     * Search the list for a spacecraft corresponding to the given ID.
     * @param id
     * @return Spacecraft or null if not found
     */
    public Spacecraft get_spacecraft(String id){
        int i = search(id);
        if(i>0){
            SpacecraftModel sc = (SpacecraftModel)member_spacecraft.get(search(id));
            return sc.get_spacecraft();
        }else if(id.equalsIgnoreCase(primary.get_id()))
            return primary.get_spacecraft();
        else
            return null;

    }
    /**
     * Return the primary spacecraft model.
     * @return Spacecraft model
     */
    public SpacecraftModel get_primarymodel(){
        return primary;
    }
    /**
     * Return the primary spacecraft.
     * @return Spacecraft.
     */
    public Spacecraft get_primary(){
        return primary.get_spacecraft();
    }
    /**
     * Get the total number of spacecraft (including the primary).
     * @return Number of spacecraft.
     */
    public int get_num_sc(){
        return member_spacecraft.size()+1;
    }
    
    /**
     * Connect two spacecraft in the formation.  
     * If 'directed' is true then "x tracks y", otherwise both are connected.
     * @param x String ID for the first spacecraft
     * @param y String ID for the second spacecraft
     * @param directed
     */
    //* if directed==true connect scx to scy - "x tracks y"
    //* else connect both
    public void connect(String x, String y, boolean directed){
        int xnum = search(x);
        int ynum = search(y);
        if(xnum > 0 && ynum > 0){
            if(directed)
                connect.set(xnum,ynum,1.0);
            else{
                connect.set(xnum,ynum,1.0);
                connect.set(ynum,xnum,1.0);
            }
                
        } else {
            if(xnum < 0)
                System.err.println("Error: '"+x+"' does not appear in formation '"+this.id+"'");
            if(ynum < 0)
                System.err.println("Error: '"+y+"' does not appear in formation '"+this.id+"'");
        }
    }
    
    /**
     * Connect two spacecraft in the formation.  
     * If 'directed' is true then "x tracks y", otherwise both are connected.
     * @param x Index for the first spacecraft
     * @param y Index for the second spacecraft
     * @param directed
     */
    //* if directed==true connect scx to scy - "x tracks y"
    //* else connect both
    public void connect(int xnum, int ynum, boolean directed){
        if(xnum > 0 && ynum > 0){
            if(directed)
                connect.set(xnum,ynum,1.0);
            else{
                connect.set(xnum,ynum,1.0);
                connect.set(ynum,xnum,1.0);
            }
                
        } else {
            if(xnum > get_num_sc())
                System.err.println("Error: '"+xnum+"' does not appear in formation '"+this.id+"'");
            if(ynum > get_num_sc())
                System.err.println("Error: '"+ynum+"' does not appear in formation '"+this.id+"'");
        }
    }

    //    private void connect(String x){
//        int xnum = search(x);
//    }

    /**
     * Searches the formation for a spacecraft with the given ID.
     * @param x The spacecraft ID.
     */
    public int search(String x){
        int out = -1;
        for(int i=0; i<member_spacecraft.size(); i++){
            SpacecraftModel tmp = (SpacecraftModel)member_spacecraft.get(i);
            if(tmp.get_id().equalsIgnoreCase(x)){
                if(out < 0)
                    out = i;
                else
                    System.err.println("Warning: multiple matching spacecraft id's");                   
            }
        }
        // else return -1
        return out;
    }
    
    public void compute_control(double t){
        primary.compute_control(t,primary.get_spacecraft().toStateVector());
        double[] X;
        double[] Xnew;
        for(int i=0; i<this.get_num_sc()-1; i++){
            SpacecraftModel sm = ((SpacecraftModel)member_spacecraft.get(i));
            X = sm.get_spacecraft().toStateVector(true);
            Xnew = sm.get_rel_state(primary);
            sm.compute_control(t,X, Xnew);	
        }
    }
}
