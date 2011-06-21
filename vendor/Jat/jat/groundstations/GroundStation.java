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

import java.io.Serializable;
import jat.cm.Constants;
import jat.constants.IERS_1996;
import jat.math.MathUtils;
import jat.matvec.data.Matrix;
import jat.matvec.data.RotationMatrix;
import jat.matvec.data.VectorN;
import jat.spacetime.ITRF;
import jat.spacetime.Time;

/**
 * Defines an ECEF stationary Ground Station using the WGS-84 Geoid.
 */
public class GroundStation implements Serializable {

	private static final long serialVersionUID = 1L;

	/**
	 * The speed of light in a vacuum, m/s.
	 */
	private static final double c = IERS_1996.c;

	private String _name = null;

	/**
	 * Geodetic latitude in radians.
	 */
	private double _lat = -9999.9;

	/**
	 * Geodetic longitude in radians.
	 */
	private double _lon = -9999.9;

	/**
	 * Height above the ellipsoid in m.
	 */
	private double _alt = -9999.9;

	/**
	 * ECEF Position Vector in m.
	 */
	private VectorN _ecefPosition = null;

	/**
	 * ECEF Velocity Vector in m/s, defaults to zero for a GroundStation.
	 */
	private VectorN _ecefVelocity = new VectorN(3);

	/**
	 * minimum elevation for viewing in radians.
	 */
	private double _minElevation = 0.0;

	/** 
	 * WGS-84 radius of the Earth in meters.
	 */
	private static double R_static = 6378.137e3;

	/**
	 * WGS-84
	 */
	private static double f_static = 1.0 / 298.257223563;


	/**
	 * create a Ground Station from latitude, longitude, HAE using WGS-84 Geoid.
	 * 
	 * @param name
	 *            String containing the name.
	 * @param lat
	 *            double containing getodetic latitude in radians.
	 * @param lon
	 *            double containing getodetic longitude in radians.
	 * @param alt
	 *            double containing height above ellipsoid in meters.
	 */
	public GroundStation(String name, double lat, double lon, double alt) {
		super();
		
		_name = name;
		_lat = lat;
		_lon = lon;
		_alt = alt;
		_ecefPosition = computeECEFposition(_lat, _lon, _alt);
	}

	/**
	 * create a Ground Station from ECEF Position Vector.
	 * 
	 * @param name
	 *            String containing the name.
	 * @param staPos
	 *            VectorN containing ECEF position
	 */
	public GroundStation(String name, VectorN staPos) {
		super();
		
		_name = name;
		// set ECEF position and velocity
		_ecefPosition = new VectorN(staPos.getArray());

		// Set the Latitude, Longitude and Height
//		double[] llh = ECEF2LLH(staPos.getArray());
		double[] llh = getLLA(staPos.getArray());
		_lat = llh[0];
		_lon = llh[1];
		_alt = llh[2];
	}

	/**
	 * ODToolbox interface create a Ground Station from ECEF Position Vector.
	 * 
	 * @param name
	 *            String containing the name.
	 * @param staPos
	 *            VectorN containing ECEF position (m)
	 */
	public GroundStation(String name, double[] staPos) {
		super();
		
		_name = name;
		// set ECEF position and velocity
		_ecefPosition = new VectorN(staPos);
		// Set the Latitude, Longitude and Height
//		double[] llh = ECEF2LLH(staPos);
		double[] llh = getLLA(staPos);
		_lat = llh[0];
		_lon = llh[1];
		_alt = llh[2];
	}

	/**
	 * create a Ground Station from latitude, longitude, HAE using WGS-84 Geoid.
	 * 
	 * @param name
	 *            String containing the name.
	 * @param lat
	 *            double containing getodetic latitude in radians.
	 * @param lon
	 *            double containing getodetic longitude in radians.
	 * @param alt
	 *            double containing height above ellipsoid in meters.
	 * @param minElevation
	 *            double containing minimum elevation in radians.
	 */
	public GroundStation(String name, double lat, double lon, double alt,
			double minElevation) {
		super();
		
		_name = name;
		_lat = lat;
		_lon = lon;
		_alt = alt;
		_minElevation = minElevation;
		_ecefPosition = computeECEFposition(_lat, _lon, _alt);
	}

	/**
	 * create a Ground Station from ECEF Position Vector.
	 * 
	 * @param name
	 *            String containing the name.
	 * @param staPos
	 *            VectorN containing ECEF position
	 * @param minElevation
	 *            double containing minimum elevation in radians.
	 */
	public GroundStation(String name, VectorN staPos, double minElevation) {
		super();
		
		_name = name;
		// set ECEF position and velocity
		_ecefPosition = new VectorN(staPos);
		// Set the Latitude, Longitude and Height
//		double[] llh = ECEF2LLH(staPos.getArray());
		double[] llh = getLLA(staPos.getArray());
		_lat = llh[0];
		_lon = llh[1];
		_alt = llh[2];
		_minElevation = minElevation;
	}

	/**
	 * ODToolbox interface create a Ground Station from ECEF Position Vector.
	 * 
	 * @param name
	 *            String containing the name.
	 * @param staPos
	 *            VectorN containing ECEF position
	 * @param minElevation
	 *            double containing minimum elevation in radians.
	 * 
	 */
	public GroundStation(String name, double[] staPos, double minElevation) {
		super();
		
		_name = name;
		// set ECEF position and velocity
		_ecefPosition = new VectorN(staPos);
		// Set the Latitude, Longitude and Height
//		double[] llh = ECEF2LLH(staPos);
		double[] llh = getLLA(staPos);
		_lat = llh[0];
		_lon = llh[1];
		_alt = llh[2];
		_minElevation = minElevation;
	}
	
	/**
	 * Copy constructor.  This performs a "deep copy" and preserves the name.
	 * @param gs Source instance of the copy.
	 */
	public GroundStation(final GroundStation gs) {
		super();
		
		this._alt = gs._alt;
		this._lat = gs._lat;
		this._lon = gs._lon;
		this._minElevation = gs._minElevation;
		this._name = gs._name;
		this._ecefPosition = gs._ecefPosition.copy();
		this._ecefVelocity = gs._ecefVelocity.copy();
	}
	
	/**
	 * Copy constructor that also sets the name of the new instance.
	 * This performs a "deep copy".
	 * @param gs Source instance of the copy.
	 */
	public GroundStation(final GroundStation gs, final String newName) {
		super();
		
		this._alt = gs._alt;
		this._lat = gs._lat;
		this._lon = gs._lon;
		this._minElevation = gs._minElevation;
		this._name = newName;
		this._ecefPosition = gs._ecefPosition.copy();
		this._ecefVelocity = gs._ecefVelocity.copy();
	}

	/**
	 * Return the name.
	 * 
	 * @return name.
	 */
	public String getName() {
		return _name;
	}

	/**
	 * Return the latitude. Note: This latitude may contain some error if this
	 * GroundStation was constructed from ECEF position due to ECEF2LLH
	 * iteration.
	 * 
	 * @return Geodetic latitude in radians.
	 */
	public double getLatitude() {
		return _lat;
	}

	/**
	 * Return the longitude. Note: This longitude may contain some error if this
	 * GroundStation was constructed from ECEF position due to ECEF2LLH
	 * iteration.
	 * 
	 * @return Geodetic longitude in radians.
	 */
	public double getLongitude() {
		return _lon;
	}

	/**
	 * return the height above the ellipsoid. Note: This height may contain some
	 * error if this GroundStation was constructed from ECEF position due to
	 * ECEF2LLH iteration.
	 * 
	 * @return height above the ellipsoid in m.
	 */
	public double getHAE() {
		return _alt;
	}
	
	/**
	 * @return The ECEF position in meters.
	 */
	public VectorN getECEFPositionVector() {
		return _ecefPosition;
	}

	/**
	 * @return The ECEF velocity in meters/second.
	 */
	public VectorN getECEFVelocityVector() {
		return _ecefVelocity;
	}

	/**
	 * @return The ECEF position in meters.
	 */
	public double[] getECEFPosition() {
		return _ecefPosition.getArray();
	}
	
	/**
	 * @return The ECEF velocity in meters/second.
	 */
	public double[] getECEFVelocity() {
		return _ecefVelocity.getArray();
	}

	/**
	 * computes the ECEF position vector.
	 * 
	 * @return ECEF position in m.
	 */
	static private VectorN computeECEFposition(final double lat, final double lon, final double alt) {
		double e2 = f_static * (2.0 - f_static); // Square of eccentricity
		double CosLat = Math.cos(lat); // (Co)sine of geodetic latitude
		double SinLat = Math.sin(lat);
		double CosLon = Math.cos(lon); // (Co)sine of geodetic latitude
		double SinLon = Math.sin(lon);
		double N;
		VectorN r = new VectorN(3);

		// Position vector

		N = R_static / Math.sqrt(1.0 - e2 * SinLat * SinLat);
		double Nh = N + alt;
		r.x[0] = Nh * CosLat * CosLon;
		r.x[1] = Nh * CosLat * SinLon;
		r.x[2] = ((1.0 - e2) * Nh) * SinLat;
		return r;
	}

	/**
	 * Computes the ECI position vector.  No separate handling of polar 
	 * motion is performed unless included in the eci2ecef argument.
	 * 
	 * @param eci2ecef
	 *            transformation matrix from ECI to ECEF
	 * @return ECI position in m.
	 */
	public VectorN getECIPosition(Matrix eci2ecef) {
		VectorN ecef = this.getECEFPositionVector();
		VectorN out = eci2ecef.transpose().times(ecef);
		return out;
	}

	/**
	 * Computes the ECI position vector.  No separate handling of polar 
	 * motion is performed unless included in the eci2ecef argument.
	 * 
	 * @param eci2ecef
	 *            transformation matrix from ECI to ECEF
	 * @return ECI position in m.
	 */
	public VectorN getECIPosition(ITRF ref) {
		VectorN ecef = this.getECEFPositionVector();
		VectorN out = ref.ECI2ECEF().times(ecef);
		return out;
	}

	/**
	 * Computes the ECI velocity vector.  This function does not include the
	 * effects of polar motion or variation in the length of the day.
	 * 
	 * @param eci2ecef
	 *            transformation matrix from ECI to ECEF
	 * @return ECI position in m.
	 */
    public VectorN getECIVelocity(final Matrix eci2ecef) {
    	Matrix ecef2eci = eci2ecef.transpose();
    	
    	VectorN recef = this.getECEFPositionVector();
    	VectorN vecef = this.getECEFVelocityVector();

    	double w = Constants.omega_e;
    	
    	// matrix equivalent of cross product of [0 0 w]:
    	Matrix wCross = new Matrix(3,3);
    	wCross.set(0, 1,-w);
    	wCross.set(1, 0, w);
    	
    	VectorN x1 = wCross.times(recef);
    	VectorN x2 = vecef.plus(x1);
    	
    	VectorN out =  ecef2eci.times(x2);
    	
    	return out;
    }


	/**
	 * Compute the transformation from ECEF to SEZ (topocentric horizon)
	 * reference frame.
	 * 
	 * @return ECEF to SEZ transformation matrix
	 */
	public Matrix ECEF2SEZ() {
		double lambda = _lon;
		double phi = _lat;

		RotationMatrix M = new RotationMatrix(3, lambda, 2,
				(Constants.pi / 2.0 - phi));
		Matrix out = new Matrix(M);
		return out;
	}



	/**
	 * Returns the azimuth and elevation angles in radians of a spacecraft
	 * relative to a ground station.  If the spacecraft and ground station are
	 * somehow coincident then the azimuth and elevation are set to 0.0.
	 * 
	 * @param satPos
	 *            double array containing the spacecraft position vector in ECEF
	 *            frame (m)
	 * @return azimuth angle (from North, positive towards East, 0 to 2pi radians), 
	 *         and elevation angle (positive up, radians)
	 */
	public double[] azEl(double[] satPos) {
		// compute vector from station to satellite in ECEF
		VectorN rho = new VectorN(satPos).minus(this.getECEFPositionVector());
		// compute transformation from ECEF to NED
		Matrix N = new RotationMatrix(2, Math.PI);
		Matrix T = N.times(ECEF2SEZ());
		// transform rho vector to NED from ECEF
		VectorN Rho = T.times(rho);
		// get components
		double x = Rho.get(0);
		double y = Rho.get(1);
		double z = Rho.get(2);
		// compute azimuth
		double az = Math.atan2(y, x);
		if (az < 0.0)
			az = az + 2.0 * Math.PI;
		// compute elevation
		double el = 0.0;
		if (Rho.mag() > 0.0) {
			el = Math.asin(-1.0*z / Rho.mag());
		}
		// create output
		double[] out = new double[2];
		out[0] = az;
		out[1] = el;
		return out;
	}


	/**
	 * Determine if a spacecraft is visible (above a given min elevation)
	 * 
	 * @param satPos
	 *            VectorN containing the spacecraft position vector in ECEF
	 *            frame
	 * @param minElevation
	 *            minimum elevation for visibility in radians
	 * @return boolean, true if elevation >= minElevation, otherwise false
	 */
	public boolean isVisible(VectorN satPos) {
		boolean out = false;
		double [] azel = this.azEl(satPos.getArray());
		if (azel[1] >= _minElevation)
			out = true;
		return out;
	}

	/**
	 * Compute the Geodetic Latitude, Longitude and Height using the algorithm
	 * from Montebruck's book "Satellite Oribts".
	 * 
	 * @param rho
	 *            double containing the ground station position vector in ECEF
	 *            frame
	 * @return double - [lat ;lon;alt] in radians, radians, m or km (depending
	 *         on input)
	 */
	public static double[] ECEF2LLH(double[] rho) {
		double theta, newRange;
		double N = 0.0;

		double x = rho[0];
		double y = rho[1];
		double z = rho[2];

		double[] out = new double[3];

		// Set the initial conditions for the iteration
		double e = Math.sqrt(1 - (1 - f_static) * (1 - f_static));
		double deltaZ = 1000; // Set high to force into loop
		double deltaZNew = e * e * z;

		// Iterate until the error is small
		while (Math.abs(deltaZNew - deltaZ) > 1e-12) {
			deltaZ = deltaZNew;

			newRange = Math.sqrt(x * x + y * y + (z + deltaZ) * (z + deltaZ));
			theta = Math.asin((z + deltaZ) / newRange);

			N = R_static
					/ Math.sqrt(1 - e * e * Math.sin(theta) * Math.sin(theta));

			deltaZNew = N * e * e * Math.sin(theta);

		}

		// Compute the Latitude, Longitude and Heigth
		out[1] = Math.atan2(y, x);
		out[0] = Math.atan2((z + deltaZNew), Math.sqrt(x * x + y * y));
		out[2] = Math.sqrt(x * x + y * y + (z + deltaZ) * (z + deltaZ)) - N;

		return out;
	}

	/**
	 * Returns inertial satellite position at time of signal transmission,
	 * r(t-tau) computed via Taylor series expansion (to avoid interpolation)
	 * 
	 * @param tau
	 *            double containing signal travel time (s)
	 * @param r
	 *            VectorN containing satellite inertial position at signal
	 *            reception time (t)
	 * @param v
	 *            VectorN containing satellite inertial velocity at signal
	 *            reception time (t)
	 * @param a
	 *            VectorN containing satellite inertial acceleration at signal
	 *            reception time (t)
	 * @return VectorN containing satellite inertial position at time of signal
	 *         transmission, r(t-tau)
	 */
	public static VectorN r_trans_time(double tau, VectorN r, VectorN v,
			VectorN a) {
		VectorN term1 = v.times(-1.0 * tau);
		VectorN term2 = a.times(0.5 * tau * tau);
		VectorN dr = term1.plus(term2);
		return r.plus(dr);
	}

	public double two_way_range(Time t, VectorN r, VectorN v, VectorN a) {
		int maxit = 500;
		double eps = 1.0E-16;

		// downlink light time iteration
		int i = 0;
		double tau_d_new = 0.0;
		double tau_d_old = t.mjd_tt();
		double rho_down = 0.0;
		double diff = Math.abs(tau_d_new - tau_d_old);
		VectorN r_trans = new VectorN(3);
		while ((diff > eps) && (i < maxit)) {
			r_trans = r_trans_time(tau_d_old, r, v, a);
			ITRF itrf = new ITRF(t, 0.0, 0.0, 0.0);
			VectorN R = this.getECIPosition(itrf);
			VectorN los = r_trans.minus(R);
			double rho = los.mag();
			tau_d_new = rho / c;
			diff = Math.abs(tau_d_new - tau_d_old);
			tau_d_old = tau_d_new;
			rho_down = rho;
			i = i + 1;
			if (i >= maxit) {
				System.out.println("transmitTime too many iterations, diff = "
						+ diff);
				if (diff > 1.0E-10) {
					throw new RuntimeException(
							"GPS_Utils.transmitTime too many iterations at t_mjd = "
									+ t.mjd_tt() + ", diff = " + diff);
				}
			}
		}
		Time tt = new Time(t.mjd_utc() + tau_d_new / 86400.0);

		// uplink time iteration
		i = 0;
		double tau_u_new = 0.0;
		double tau_u_old = t.mjd_tt();
		double rho_up = 0.0;
		diff = Math.abs(tau_u_new - tau_u_old);
		while ((diff > eps) && (i < maxit)) {
			Time ttt = tt.plus(-1.0 * tau_u_old);
			ITRF itrf = new ITRF(ttt, 0.0, 0.0, 0.0);
			VectorN R = this.getECIPosition(itrf);
			VectorN los = r_trans.minus(R);
			double rho = los.mag();
			tau_u_new = rho / c;
			diff = Math.abs(tau_u_new - tau_u_old);
			tau_u_old = tau_u_new;
			rho_up = rho;
			i = i + 1;
			if (i >= maxit) {
				System.out.println("transmitTime too many iterations, diff = "
						+ diff);
				if (diff > 1.0E-10) {
					throw new RuntimeException(
							"GPS_Utils.transmitTime too many iterations at t_mjd = "
									+ t.mjd_tt() + ", diff = " + diff);
				}
			}
		}
		double range = 0.5 * (rho_up + rho_down);
		return range;
	}

	/**
	 * Calculates the Geodetic latitude, longitude, and altitude from the input
	 * ECEF Position (m). <br>
	 * This method starts with the output from the geod class and uses the
	 * Jacobian matrix (partial derivatives of the exact solution for the
	 * ECEF position vector with respect to latitude, longitude, and height)
	 * to find an accurate solution using Newton's iteration method.<br>
	 * Note: coded by Keith Speckman July 1, 2008
	 * @param rho ECEF Position in m
	 * @return Vector containing the latitude (rad), longitude (rad), and 
	 * altitude (m), followed by the iteration error of the x y z ECEF 
	 * position (m).  This error is the difference between the equivalent
	 * ECEF position from the iteration and rho.
	 */
	public static double[] getLLA(double [] rho){
		double rtd = Constants.rad2deg;    // radians to degrees conversion
		double dtr = Constants.deg2rad;    // degress to radians conversion
		double LatLonAltArray[] = geod(rho); // m in, rad & m out
		double[] out = new double[6];
		Matrix Jinv = new Matrix(3);
		VectorN L = new VectorN(3);
		VectorN L2 = new VectorN(3);

		L.x[0] = LatLonAltArray[0] * rtd; // in deg
		L.x[1] = LatLonAltArray[1] * rtd; // in deg
		L.x[2] = LatLonAltArray[2] / 1000.0; // in km

		//Matrix Jinv = new getJacobian(L.x);

		VectorN Xinit = new VectorN(rho); // in m
		Xinit = Xinit.divide(1000.0); // now in km
		VectorN X = Xinit; // in km
		GroundStation GS = new GroundStation("tmp",LatLonAltArray[0],LatLonAltArray[1],
											 LatLonAltArray[2]); 
		VectorN X2 = GS.getECEFPositionVector().divide(1000.0); // in km

		int i = 0;
		while (Math.abs(X2.x[0]-X.x[0]) > 1e-10 ||       
		       Math.abs(X2.x[1]-X.x[1]) > 1e-10 ||
		       Math.abs(X2.x[2]-X.x[2]) > 1e-10  ) {  // km tol check
			i = i + 1;
			if (i == 12) {
				L.x[0] = LatLonAltArray[0] * rtd;
				L.x[1] = LatLonAltArray[1] * rtd;
				L.x[2] = LatLonAltArray[2] / 1000.0;
			        GS = new GroundStation("tmp",LatLonAltArray[0],LatLonAltArray[1],
						             LatLonAltArray[2]);
			        X2 = GS.getECEFPositionVector().divide(1000.0); // in m
				System.out.println("\nWARNING - LLA calculation failed to converge.");
				break;
			}
			X = X2;
			Jinv = getJacobian(L); // i/o: deg, km
			L2 = Jinv.times(Xinit.minus(X));
			L = L2.plus(L);
			GS = new GroundStation("tmp",L.x[0]*dtr,L.x[1]*dtr,L.x[2]*1000.0);
			X2 = GS.getECEFPositionVector().divide(1000.0); // in km
		}
		out[0] = L.x[0]*dtr; // back to rad
		out[1] = L.x[1]*dtr; // back to rad
		out[2] = L.x[2]*1000.0; // back to m
		out[3] = X2.x[0]*1000.0 - rho[0];
		out[4] = X2.x[1]*1000.0 - rho[1];
		out[5] = X2.x[2]*1000.0 - rho[2];

		return out;
	}

	/**
	 * Compute the Geodetic Latitude, Longitude and Height using K.M. Borkowski's
	 * algorithm from " Accurate Algorithms to Transform Geocentric to Geodetic 
	 * Coordinates," Bulletin Geodesique, Vol 63, pp. 50-56, 1989. Equations
	 * can be found in Seidelmann, Explanatory Supplement to the Astronomical
	 * Almanac, 1992, pp. 206-207. 
	 * Note: this solution is exact (no iteration) but contains a small error.
	 * 
	 * @param rho
	 *            double containing the ground station position vector in ECEF
	 *            frame (m)
	 * @return double - [lat lon alt] in radians, radians, m.
	 */
	public static double[] geod(double [] rho){
		double a = R_static; // meters
		double X = rho[0];
		double Y = rho[1];
		double Z = rho[2];
		
		double lat = 0.0;
		double lon = Math.atan2(Y, X);
		double height = 0.0;
		double [] out = new double[3];
		
		double r = Math.sqrt(X*X + Y*Y);
		if (r == 0.0) {
			// no X or Y?, then latitude is +/- PI/2 by definition
			out[0] = (Z >= 0) ? Math.PI/2 : -Math.PI/2;
			out[1] = 0.0; // longitude is 0 by convention
			out[2] = Math.abs(Z)-a; 
			return out;
		}
		double a2 = a*a;
		double b = MathUtils.sign(a-a*f_static, Z);
		double b2 = b*b;
		double E = (b*Z - (a2 - b2))/(a*r);
		double F = (b*Z + (a2 - b2))/(a*r);
		double P = (E*F+1.0)*4.0/3.0;
		double Q = (E*E-F*F)*2.0;
		double D = P*P*P+Q*Q;
		double v = Math.pow(Math.sqrt(D) - Q, (1.0/3.0)) - Math.pow(Math.sqrt(D) + Q, (1.0/3.0));
		if (D < 0.0) {
			v = 2.0 * Math.sqrt(-P) * Math.cos(Math.acos(Q*Math.sqrt(-P)/P)/3.0);
		}
		double G = 0.5*(E+Math.sqrt(E*E+v));
		double t = Math.sqrt(G*G+(F-v*G)/(G+G-E))-G;
		lat = Math.atan(((1.-t*t)*a)/(2*b*t));
		height = (r-a*t)*Math.cos(lat)+(Z-b)*Math.sin(lat);
		out[0] = lat;
		out[1] = lon;
		out[2] = height;
		return out;	
	}

	/**
	 * Compute the Jacobian (partial derivatives matrix) inverse.  The Jabian (J)
	 * is the 3x3 matrix of partial derivatives of the ECEF state vector with respect
	 * to latitude, longitude, and height.  The Jacobian is defined as follows:
	 *        |  drho1/dlat     drho1/dlon    drho1/dheight |
	 *        |  drho2/dlat     drho2/dlon    drho2/dheight |
	 *        |  drho3/dlat     drho3/dlon    drho3/dheight |
	 * This algorithm outputs the invers of that matrix
	 * Note: this solution is exact.
	 * Note: coded by Keith Speckman - June 30, 2008
	 * 
	 * @param L
	 *            VectorN [3] containing the [lat, lon, alt] (deg,km)
	 * @return Jinv
	 *            Matrix [3][3] - inverse of Jacobian matrix
	 */
	private static Matrix getJacobian(VectorN L) {
		double[][] J = new double[3][3];
		double R_equ=6378.137;             // Mean Earth Radius (km)
		double f=1.0 / 298.257223563;
		double e2= f*(2.0-f);              // Square of eccentricity
		double dtr = Constants.deg2rad;    // degrees to radians conversion

		double lat = L.x[0];
		double lon = L.x[1];
		double height = L.x[2];

		double CosLat = Math.cos(lat*dtr); // (Co)sine of geodetic latitude
		double SinLat = Math.sin(lat*dtr);
		double CosLon = Math.cos(lon*dtr); // (Co)sine of geodetic latitude
		double SinLon = Math.sin(lon*dtr);

		double N = R_equ / Math.sqrt(1.0-e2*SinLat*SinLat);
		double Nh = N + height;

		// Calculate Jacobian (partial Xecef / partial Lat,Lon,Height)
		J[0][0] = (-SinLat*CosLon*Nh+R_equ*e2*SinLat*CosLat*CosLat*CosLon*
			   Math.pow(1-e2*SinLat*SinLat,-3/2))*dtr;
		J[0][1] = -SinLon*CosLat*Nh*dtr;
		J[0][2] = CosLat*CosLon;
		J[1][0] = (-SinLat*SinLon*Nh+R_equ*e2*SinLat*CosLat*CosLat*SinLon*
			   Math.pow(1-e2*SinLat*SinLat,-3/2))*dtr;
		J[1][1] = CosLon*CosLat*Nh*dtr;
		J[1][2] = CosLat*SinLon;
		J[2][0] = (1-e2)*(CosLat*Nh+R_equ*e2*SinLat*SinLat*CosLat*
			   Math.pow(1-e2*SinLat*SinLat,-3/2))*dtr;
		J[2][1] = 0;
		J[2][2] = (1-e2)*SinLat;

		Matrix JMat = new Matrix(J,3,3);
		Matrix Jinv = JMat.inverse();

		return Jinv;
	}
	
	public static void main(String[] args) {

		// Create a Ground station using latitude longitude and height
		double staLat = 55 * MathUtils.DEG2RAD;
		double staLon = -120 * MathUtils.DEG2RAD;
		double staHeight = 25;
		GroundStation sta1 = new GroundStation("tmp1", staLat, staLon,
				staHeight);

		// Get the ECEF coordinates
		VectorN sta1Pos = sta1.getECEFPositionVector();
		sta1Pos.print("sta1Pos");

		// Convert it back to make sure we get the same answer
		// VectorN sta1Pos = new VectorN(1917032.190, 6029782.349, -801376.113);
		GroundStation station = new GroundStation("tmp", sta1Pos);
		double lat = station.getLatitude();
		double lon = station.getLongitude();
		double height = station.getHAE();
		System.out.println("Latitude: " + lat * MathUtils.RAD2DEG);
		System.out.println("Longitude: " + lon * MathUtils.RAD2DEG);
		System.out.println("Height: " + height);
		
		double[] xtest = new double[3];
		xtest[0] = 2050728.95875797;
		xtest[1] = 5634331.50760922;
		xtest[2] = 2167730.76088181;
		double[] llh = getLLA(xtest);
		System.out.println("Latitude2: " + llh[0] * MathUtils.RAD2DEG);
		System.out.println("Longitude2: " + llh[1] * MathUtils.RAD2DEG);
		System.out.println("Height2: " + llh[2]);

	}

}
