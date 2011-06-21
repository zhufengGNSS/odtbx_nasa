package jat.ground_tracking;

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

import java.io.FileOutputStream;
import java.io.PrintStream;

import jat.math.MathUtils;
import jat.matvec.data.VectorN;
import jat.spacetime.Time;
import jat.groundstations.GroundStation;

/**
 * IonosphericDelayModel is a model of ranging signal delay due to ionospheric 
 * plasma effects.  An instance of this class represents a single GroundStation
 * site that all calculations are based upon.  The model's coefficients can be
 * tailored for each calculation call to computeDelay().<br><br>
 * 
 * This is a first-order approximation model from Burkhart:
 * <b><i>
 * Burkhart, P. D., "Adaptive Orbit Determination For Interplanetary Spacecraft",
 * PhD Dissertation, The University Of Texas At Austin, May 1995, section 2.2.3
 * "Ionosphere Model", pp. 22-27.
 * </i></b>
 * It only considers spacecraft elevation angle and diurnal ionospheric 
 * fluctuations.  
 * <br><br> 
 * Burkhart notes:
 * <br><i>
 * The model used is designed to best fit the ionospheric delay in the 
 * afternoon period, when the delay is largest, and to be able to more 
 * closely model the characteristics for specific latitudes of interest 
 * by varying the coefficients in the model.
 * </i><br><br>
 * This is not a worldwide time-delay model.
 */
public class IonosphericDelayModel {

	/**
	 * Mean height of the ionosphere in meters
	 */
	private double _h = 350e3;

	/**
	 * Nighttime ionospheric behavior coefficient (Hz^2 * m)
	 */
	private double _N = 3.1928e18;

	/**
	 * Daytime (diurnal) ionospheric behavior coefficient (Hz^2 * m)
	 */
	private double _D = 3.1928e19;

	/**
	 * Period of Diurnal Variation (hours)
	 */
	private double _Pbar = 32;

	/**
	 * Phase of Diurnal Variation (hours)
	 */
	private double _thetaD = 14;

	/**
	 * Frequency of signal (Hz)
	 */
	private double _freq = 1.6e9;

	/**
	 * The GroundStation that represents the basis for this instance's
	 * calculations.
	 */
	private final GroundStation _gs;

	/**
	 * Creates an ionospheric model for a GroundStation. A default frequency of
	 * 1.6Ghz is used.
	 * 
	 * @param gs
	 *            GroundStation object
	 */
	public IonosphericDelayModel(final GroundStation gs) {
		super();
		// Use a deep copy of the given GroundStation to preserve semantics
		this._gs = new GroundStation(gs);
	}

	/**
	 * Creates an ionospheric model for a GroundStation. A default frequency of
	 * 1.6Ghz is used.
	 * 
	 * @param staPos
	 *            double array containing the ECEF ground station position (m)
	 */
	public IonosphericDelayModel(double[] staPos) {
		super();
		this._gs = new GroundStation("tmp", staPos);
	}

	/**
	 * Creates an ionospheric model from the geodetic latitude, longitude and
	 * height. A default frequency of 1.6Ghz is used.
	 * 
	 * @param lat
	 *            double geodetic latitude (rad)
	 * @param lon
	 *            double geodetic longitude (rad)
	 * @param height
	 *            double geodetic height (m)
	 */
	public IonosphericDelayModel(double lat, double lon, double height) {
		super();
		this._gs = new GroundStation("tmp", lat, lon, height);
	}

	/**
	 * Creates an ionospheric model for a GroundStation with a specified 
	 * frequency.
	 * 
	 * @param gs
	 *            GroundStation object
	 * @param f
	 *            double containing frequency in Hz
	 */
	public IonosphericDelayModel(final GroundStation gs, double f) {
		super();
		// Use a deep copy of the given GroundStation to preserve semantics
		this._gs = new GroundStation(gs);
		this._freq = f;
	}

	/**
	 * Creates an ionospheric model for a GroundStation with a specified 
	 * frequency.
	 * 
	 * @param staPos
	 *            double array containing the ECEF ground station position (m)
	 * @param f
	 *            double containing frequency in Hz
	 */
	public IonosphericDelayModel(double[] staPos, double f) {
		super();
		this._gs = new GroundStation("tmp", staPos);
		this._freq = f;
	}

	/**
	 * Creates an ionospheric model from the geodetic latitude, longitude and
	 * height of a ground station with a specified frequency.
	 * 
	 * @param lat
	 *            double ground station geodetic latitude (radians)
	 * @param lon
	 *            double ground station geodetic longitude (radians)
	 * @param height
	 *            double ground station geodetic height (m)
	 * @param f
	 *            double containing frequency in Hz
	 */
	public IonosphericDelayModel(double lat, double lon, double height, double f) {
		super();
		this._gs = new GroundStation("tmp", lat, lon, height);
		this._freq = f;
	}

	/**
	 * Computes the slant angle and ionospheric delay between a satellite and
	 * the ground station.
	 * 
	 * @param T
	 *            double local time of the subionospheric point
	 * @param satPos
	 * 			  double array containing the spacecraft position vector in the
	 *            ECEF frame (m)
	 * @returns double array containing slant angle (radians) and ionospheric
	 *          delay (m)
	 */
	public double[] computeDelay(final Time T, final double[] satPos) {
		double[] azel = this._gs.azEl(satPos);

		// azimuth and elevation angles (rad) from the tracking station
		// to the satellite:
		double az = azel[0];
		double el = azel[1];

		// tracking station geodetic latitude and longitude (rad):
		double staLat = this._gs.getLatitude();
		double staLon = this._gs.getLongitude();

		// Earth radius (m)
		final double re = jat.constants.WGS84.R_Earth;

		// Determine Slant Height
		// Burkhart Eq 2.13, "zenith angle", p. 26 (radians):
		double slantAngle = Math.asin(re / (re + this._h) * Math.cos(el));

		// Determine the latitude of the sub-ionospheric point (rad)
		// Burkhart (unnumbered eqn), p. 25:
		double ionoLat = Math.asin(Math.sin(staLat) * Math.sin(el + slantAngle)
				+ Math.cos(staLat) * Math.cos(az) * Math.cos(el + slantAngle));

		// Determine the longitude of the sub-ionospheric point (rad)
		// Burkhart (unnumbered eqn), p. 25:
		double ionoLong = staLon
				+ Math.asin(Math.sin(az) * Math.cos(el + slantAngle)
						/ Math.cos(ionoLat));

		// Determine the local time, number of hours into current day
		double UTCHours = T.secOfDay() / 3600;

		// local time of the subionospheric point (hours),
		// Burkhart (unnumbered eqn), p. 25, constrained to be 0-24 hours:
		double tLocal = MathUtils.mod((MathUtils.RAD2DEG * ionoLong) / 15 + UTCHours, 24);

		// Burkhart Eq 2.12, p. 22:
		double chi = Math.abs((2 * Math.PI / this._Pbar) * (tLocal - this._thetaD));
		if (chi > (Math.PI) / 2) {
			chi = Math.PI / 2;
		}

		// Burkhart Eq 2.11, "1st order delay approximation", p. 22, except
		// the div by freq^2 is not explicitly in that section. Campbell
		// references Burkhart but explicitly includes the frequency and
		// N & D coefficients' units.
		double ionoDelay = (1 / Math.cos(slantAngle))
				* (this._N + this._D * Math.cos(chi)) / (this._freq * this._freq);

		double[] out = new double[2];
		out[0] = slantAngle;
		out[1] = ionoDelay;
		return out;
	}

	/**
	 * @return The mean height of the ionosphere in meters
	 */
	public double get_h() {
		return this._h;
	}

	/**
	 * @return The nighttime ionospheric behavior coefficient (Hz^2 * m)
	 */
	public double get_N() {
		return this._N;
	}

	/**
	 * @return The daytime ionospheric behavior coefficient (Hz^2 * m)
	 */
	public double get_D() {
		return this._D;
	}

	/**
	 * @return The Period of Diurnal Variation (hours)
	 */
	public double get_Pbar() {
		return this._Pbar;
	}

	/**
	 * @return The Phase of Diurnal Variation (hours)
	 */
	public double get_thetaD() {
		return this._thetaD;
	}

	/**
	 * @return The Frequency of signal (Hz)
	 */
	public double get_freq() {
		return this._freq;
	}

	/**
	 * @return A copy of the GroundStation instance that represents the basis
	 *         for calculations of this model instance.
	 */
	public final GroundStation get_gs() {
		return new GroundStation(this._gs);
	}

	/**
	 * Sets the mean height of the ionosphere in meters.
	 * 
	 * @param h
	 *            double The new mean ionospheric height (m)
	 */
	public void set_h(double h) {
		this._h = h;
	}

	/**
	 * Sets the Nighttime ionospheric behavior coefficient.
	 * 
	 * @param N
	 *            double The new nighttime ionospheric behavior coefficient
	 *            (Hz^2 * m)
	 */
	public void set_N(double N) {
		this._N = N;
	}

	/**
	 * Sets the Daytime ionospheric behavior coefficient.
	 * 
	 * @param D
	 *            double The new Daytime ionospheric behavior coefficient
	 *            (Hz^2 * m)
	 */
	public void set_D(double D) {
		this._D = D;
	}

	/**
	 * Sets the period of the dinural variation in hours.
	 * 
	 * @param Pbar
	 *            double The new period of the dinural variation (hours)
	 */
	public void set_Pbar(double Pbar) {
		this._Pbar = Pbar;
	}

	/**
	 * Sets the phase of the dinural variation in hours.
	 * 
	 * @param thetaD
	 *            double The new phase of the dinural variation (hours)
	 */
	public void set_thetaD(double thetaD) {
		this._thetaD = thetaD;
	}

	/**
	 * Sets the frequency of the signal in Hz.
	 * 
	 * @param freq
	 *            double The new frequency of the signal (Hz)
	 */
	public void set_freq(double freq) {
		this._freq = freq;
	}

	/**
	 * Example routine.
	 * 
	 * This routine sets a ground station at 0 lat, 0 lon, and 10m alt, then
	 * varies a satellite position in the ECEF X-Y plane.  The satellite
	 * x,y,z position and corresponding delay are output to the file 
	 * "iono.txt".  The defaulted class values are used.
	 * 
	 * @param args (not used)
	 */
	public static void main(String[] args) {

		double x, y, z;
		VectorN state;

		// Set the Ground Station Parameters
		double lat = 0; // rad
		double lon = 0; // rad
		double height = 10; // m

		// Set some variables for file output
		FileOutputStream out;
		PrintStream p;

		try {
			out = new FileOutputStream("iono.txt");
			p = new PrintStream(out);

			for (double deg = -90; deg < 90; deg = deg + .01) {

				// Set the spacecraft's state (m)
				x = jat.constants.WGS84.R_Earth * Math.cos(deg * Math.PI / 180);
				y = jat.constants.WGS84.R_Earth * Math.sin(deg * Math.PI / 180);
				z = 0;
				state = new VectorN(x, y, z);

				// Set the time, assuming a 90 minute orbit
				Time T = new Time(54103.0 + 1.7361e-006);

				// Create the Ionospheric Delay Model
				IonosphericDelayModel iono = new IonosphericDelayModel(lat,
						lon, height);
				double[] temp = iono.computeDelay(T, state.getArray());
				double delay = temp[1];

				p.println(x + " " + y + " " + z + " " + delay);

			}
			p.close();
		}

		catch (Exception e) {
			System.err.println("Error writing to file");
		}

	}

}
