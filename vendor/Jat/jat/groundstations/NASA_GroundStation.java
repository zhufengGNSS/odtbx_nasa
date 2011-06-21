/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2007 United States Government as represented by the
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
 * File Created on May 20, 2007
 */
package jat.groundstations;

import jat.cm.Constants;
import jat.math.MathUtils;
import jat.matvec.data.Matrix;
import jat.matvec.data.RotationMatrix;
import jat.matvec.data.VectorN;

public class NASA_GroundStation extends GroundStation{
	
	private static final long serialVersionUID = 1L;

	private String _STDN = null;
	
	private String _location = null;
	
	private String _NASA = null;
	
    /** latitude in radians.
     */
	private double _lat = -9999.9;
	
    /** longitude in radians.
     */
	private double _lon = -9999.9;
	
    /** height above the ellipsoid in m.
     */	
	private double _alt = -9999.9;
	
    /** minimum elevation for viewing in radians.
     */		
	private double _minElevation = 0.0;
	
	
	// WGS-84
    private double R_equ = 6378.137e3;

    // WGS-84
    private double f = 1.0/298.257223563;	
	
    /** create a Ground Station from latitude, longitude, HAE using WGS-84 Geoid.
     * @param stdn String containing the STDN code.
     * @param loc String containing the location.
     * @param nasa String containing the NASA number.
     * @param lat  double containing latitude in radians.
     * @param lon  double containing longitude in radians.
     * @param alt  double containing height above ellipsoid in meters.
     */
	public NASA_GroundStation(String stdn, String loc, String nasa, double lat, double lon, double alt){
		super(loc, lat, lon, alt);
		_STDN = stdn;
		_location = loc;
		_NASA = nasa;
	}
	
    /** create a Ground Station from an ECEF Position Vector.
     * @param stdn String containing the STDN code.
     * @param loc String containing the location.
     * @param nasa String containing the NASA number.
     * @param r Geocentric position in m.
     */
    public NASA_GroundStation(String stdn, String loc, String nasa, VectorN r){
    	super(loc, r);
    	_STDN = stdn;
    	_location = loc;
		_NASA = nasa;
    	
    }
	   /** ODToolbox interface 
	    * create a Ground Station from ECEF Position Vector.
  * @param name String containing the name.
  * @param staPos VectorN containing ECEF position
  */
	public NASA_GroundStation(String stdn, String loc, String nasa, double[] staPos){
    	super(loc, staPos);
    	_STDN = stdn;
    	_location = loc;
		_NASA = nasa;
	}    

    /** create a Ground Station from latitude, longitude, HAE using WGS-84 Geoid.
     * @param stdn String containing the STDN code.
     * @param loc String containing the location.
     * @param nasa String containing the NASA number.
     * @param lat  double containing latitude in radians.
     * @param lon  double containing longitude in radians.
     * @param alt  double containing height above ellipsoid in meters.
     * @param minElevation double containing minimum elevation in radians.
     */
	public NASA_GroundStation(String stdn, String loc, String nasa, double lat, double lon, double alt, double minElevation){
		super(loc, lat, lon, alt, minElevation);
		_STDN = stdn;
		_location = loc;
		_NASA = nasa;
	}
	
    /** create a Ground Station from an ECEF Position Vector.
     * @param stdn String containing the STDN code.
     * @param loc String containing the location.
     * @param nasa String containing the NASA number.
     * @param r Geocentric position in m.
     * @param minElevation double containing minimum elevation in radians.
     */
    public NASA_GroundStation(String stdn, String loc, String nasa, VectorN r, double minElevation){
    	super(loc, r, minElevation);
    	_STDN = stdn;
    	_location = loc;
		_NASA = nasa;
    	
    }
	   /** ODToolbox interface 
	    * create a Ground Station from ECEF Position Vector.
  * @param name String containing the name.
  * @param staPos VectorN containing ECEF position
  * @param minElevation double containing minimum elevation in radians.
  */
	public NASA_GroundStation(String stdn, String loc, String nasa, double[] staPos, double minElevation){
    	super(loc, staPos, minElevation);
    	_STDN = stdn;
    	_location = loc;
		_NASA = nasa;
	}    
	
	
    /** Return the STDN code.
     * @return STDN code.
     */
    public String getSTDN(){
        return _STDN;
    }    

    /** Return the NASA number.
     * @return NASA number.
     */
    public String getNASA(){
        return _NASA;
    }        

    /** Return the location.
     * @return location.
     */
    public String getLocation(){
        return _location;
    }      





}
