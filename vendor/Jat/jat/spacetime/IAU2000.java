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
 */

package jat.spacetime;

import jat.math.MathUtils;
import jat.matvec.data.*;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 * A class for coordinate system conversions to accomplish the IAU 2000
 * resolutions. Adapted from C++ and MATLAB code from David Vallado.<br>
 * References include:
 * <ul>
 * <li><i>Fundamentals of Astrodynamics and Applicaitons, Third Ed.</i>,
 * Vallado, David A., Microcosm Press and Springer, 2007. (See sections 3.7,
 * 3.7.1 and 3.7.2)</li>
 * <li>The original code is at: http://celestrak.com/software/vallado-sw.asp.</li>
 * <li>http://www.centerforspace.com/downloads/files/pubs/AAS-06-134.pdf</li>
 * </ul>
 * 
 * @author abrown
 * 
 */
public class IAU2000 {

	/**
	 * IAU Theory Enumeration
	 * Enumeration options (from ASTREDUC.H, line 40)
	 */
	public enum eOpt {
		/**
		 * e80 is the IAU-1980 theory
		 */
		e80,

		/**
		 * e96 is the IAU-1996 Theory
		 */
		e96,

		/**
		 * e00a is the IAU-2000A Theory
		 */
		e00a,

		/**
		 * e00b is the IAU-2000B Theory
		 */
		e00b
	}

	/**
	 * IAU80 constants (rad) (from ASTREDUC.H, line 56)
	 */
	public static class iau80data {
		
		/**
		 * Integers for fk5 1980
		 */
		public int iar80[][] = new int[107][6];
		
		/**
		 * Reals for fk5 1980
		 */
		public double rar80[][] = new double[107][5];
	}

	/**
	 * A class to hold the outputs from precess() since Java is pass-by-value.
	 * This class is not immutable because its members are used in internal
	 * computation by precess().  Note that as designed precess() can create
	 * two different sets of outputs.
	 * 
	 * The convention is:
	 * 
	 * If the type is e80 then the Matrix prec will be set but the other
	 * members will be invalid.
	 * 
	 * If the type is not e80, then the Matrix prec will be invalid but the
	 * other members will be set.
	 */
	public static class precessOut {
		
		/**
		 * The algorithm type used to create this output.
		 */
		public eOpt type = null;

		/**
		 * cannonical precession angle rad (00 only)
		 */
		public double psia;

		/**
		 * cannonical precession angle rad (00 only)
		 */
		public double wa;

		/**
		 * cannonical precession angle rad (00 only)
		 */
		public double epsa;

		/**
		 * cannonical precession angle rad (00 only)
		 */
		public double chia;

		/**
		 * transformation matrix for mod - j2000 (80 only) (matrix converting
		 * from "mod" to gcrf)
		 */
		public Matrix prec = new Matrix(3, 3);
	}

	/**
	 * A class to hold the outputs from nutation() since Java is pass-by-value.
	 * This class is not immutable because its members are used in internal
	 * computation by nutation().
	 */
	public static class nutationOut {

		/**
		 * nutation in longiotude angle rad
		 */
		public double deltapsi;

		/**
		 * (undocumented in C++ code)
		 */
		public double deltaeps;

		/**
		 * true obliquity of the ecliptic rad
		 */
		double trueeps;

		/**
		 * mean obliquity of the ecliptic rad
		 */
		double meaneps;

		/**
		 * (rad)
		 */
		double omega;

		/**
		 * transform matrix for tod
		 */
		Matrix nut = new Matrix(3, 3);
	}

	/**
	 * A class to hold the outputs from fundarg() since Java is pass-by-value.
	 * This class is not immutable because its members are used in internal
	 * computation by fundarg().
	 */
	public static class fundargOut {

		/**
		 * Maximum length of planetlon.
		 */
		public static final int PLANETLONLENGTH = 8;

		/**
		 * delaunay element (rad)
		 */
		public double l;

		/**
		 * delaunay element (rad)
		 */
		public double l1;

		/**
		 * delaunay element (rad)
		 */
		public double f;

		/**
		 * delaunay element (rad)
		 */
		public double d;

		/**
		 * delaunay element (rad)
		 */
		public double omega;

		/**
		 * Calculated planetary longitude (rad). The index order is planetary
		 * order from the sun, e.g. [0]=Mercury, [1]=Venus, ... [7]=Neptune.
		 * Length of PLANETLONLENGTH.
		 */
		public double[] planetlon = new double[PLANETLONLENGTH];

		/**
		 * (undocumented)
		 */
		public double precrate;
	}

/**
 *  This function transforms a vector between the earth fixed (ITRF) frame, and
 *  the GCRF mean equator mean equinox.  This is the preferrred method to
 *  accomplish the new IAU 2000 resolutions and uses the eop corrections.
 *  This implementation can transform from ITRF to GCRF or from GCRF to ITRF.
 *  The direction is determined by the itrf2gcrf boolean parameter.  The
 *  VectorN parameters are either inputs or outputs accordingly while the
 *  trans Matrix parameter is always output. <br>
 *  
 *  This behavior is defined as:<br>
 *  <ol>
 *  <li>When <b>itrf2gcrf is true</b> the ritrf, virtf, and aitrf VectorN parameters
 *  are inputs and not altered.  The rgcrf, vgcrf, and agcrf VectorN parameters
 *  are altered and the trans Matrix is set from ITRF to GCRF.</li>
 *  <li>When <b>itrf2gcrf is false</b> the rgcrf, vgcrf, and agcrf VectorN parameters
 *  are inputs and not altered.  The ritrf, virtf, and aitrf VectorN parameters
 *  are altered and the trans Matrix is set from GCRF to ITRF.</li>
 *  </ol> <br>
 *
 * @param iau80rec iau80data containing the IAU80 constants (rad)
 * @param ttt Julian centuries of TT (centuries)
 * @param jdut1 Julian Date of UT1 (days from 4713 BC)
 * @param lod excess length of day (sec)
 * @param xp polar motion coefficient (rad)
 * @param yp polar motion coefficient (rad)
 * @param eqeterms  terms for ast calculation (0 or 2)
 * @param ddpsi delta psi correction to gcrf (rad)
 * @param ddeps delta eps correction to gcrf (rad)
 * @param itrf2gcrf If true then the transformation is from the ITRF to GCRF
 * systems, if false, then the transformation is from the GCRF to ITRF systems;
 * Note the VectorN parameters are adjusted accordingly.
 * @param ritrf (input or output) 3x1 position vector earth fixed (km)
 * @param vitrf (input or output) 3x1 velocity vector earth fixed (km/s)
 * @param aitrf (input or output) 3x1 acceleration vector earth fixed (km/s^2)
 * @param rgcrf (input or output) 3x1 position vector GCRF (km)
 * @param vgcrf (input or output) 3x1 velocity vector GCRF (km/s)
 * @param agcrf (input or output) 3x1 acceleration vector GCRF (km/s^2)
 * @param trans (output) 3x3 coordinate transformation matrix from one system to
 * the other, this transformation corresponds to the itrf2gcrf boolean setting.
 * @throws IllegalArgumentException If eqetemrms is out of range of if a VectorN
 * or Matrix parameter is not 3x1 or 3x3.
 */
	public static void iau76fk5_itrf_gcrf(final iau80data iau80rec,
			final double ttt, final double jdut1, final double lod,
			final double xp, final double yp, final int eqeterms,
			final double ddpsi, final double ddeps, final boolean itrf2gcrf,
			VectorN ritrf, VectorN vitrf, VectorN aitrf, VectorN rgcrf,
			VectorN vgcrf, VectorN agcrf, Matrix trans)
			throws IllegalArgumentException {


		/* Note, portions of this code come from ASTREDUC.CPP line 58 and 
		 * MATLAB's ecef2eci.m and eci2efec.m.
		 * ----------------------------------------------------------------------------
		 *
		 *                           function iau76fk5_itrf_gcrf
		 *
		 *  this function transforms a vector between the earth fixed (itrf) frame, and
		 *    the gcrf mean equator mean equinox. this is the preferrred method to
		 *    accomplish the new iau 2000 resolutions and uses the eop corrections
		 *
		 *  author        : david vallado                  719-573-2600   23 nov 2005
		 *
		 *  revisions
		 *
		 *  inputs          description                    range / units
		 *    ritrf       - position vector earth fixed    km
		 *    vitrf       - velocity vector earth fixed    km/s
		 *    aitrf       - acceleration vector earth fixedkm/s2
		 *    direct      - direction of transfer          eFrom, eTo
		 *    iau80rec    - record containing the iau80 constants rad
		 *    ttt         - julian centuries of tt         centuries
		 *    jdut1       - julian date of ut1             days from 4713 bc
		 *    lod         - excess length of day           sec
		 *    xp          - polar motion coefficient       rad
		 *    yp          - polar motion coefficient       rad
		 *    eqeterms    - terms for ast calculation      0,2
		 *    ddpsi       - delta psi correction to gcrf   rad
		 *    ddeps       - delta eps correction to gcrf   rad
		 *    nutopt      - nutation option                calc 'c', read 'r'
		 *    deltapsi    - nutation angle                 rad
		 *    deltaeps    - nutation angle                 rad
		 *
		 *  outputs       :
		 *    rgcrf       - position vector gcrf            km
		 *    vgcrf       - velocity vector gcrf            km/s
		 *    agcrf       - acceleration vector gcrf        km/s2
		 *    trans       - matrix for pef - gcrf
		 *
		 *  locals        :
		 *    trueeps     - true obliquity of the ecliptic rad
		 *    meaneps     - mean obliquity of the ecliptic rad
		 *    omega       -                                rad
		 *    prec        - matrix for mod - gcrf
		 *    nut         - matrix for tod - mod
		 *    st          - matrix for pef - tod
		 *    stdot       - matrix for pef - tod rate
		 *    pm          - matrix for itrf - pef
		 *
		 *  coupling      :
		 *   precess      - rotation for precession
		 *   nutation     - rotation for nutation
		 *   sidereal     - rotation for sidereal time
		 *   polarm       - rotation for polar motion
		 *
		 *  references    :
		 *    vallado       2007, 228
		 *
		 * --------------------------------------------------------------------------- */
		
		// Check eqeterms input
		if ((eqeterms != 0) || (eqeterms != 2)) {
			throw new IllegalArgumentException("eqeterms must be 0 or 2, not "
					+ eqeterms);
		}
		// check Matrix and VectorN argument dimensions
		if (ritrf.length != 3) {
			throw new IllegalArgumentException("ritrf length mus tbe 3, not "
					+ ritrf.length);
		}
		if (vitrf.length != 3) {
			throw new IllegalArgumentException("vitrf length mus tbe 3, not "
					+ vitrf.length);
		}
		if (aitrf.length != 3) {
			throw new IllegalArgumentException("aitrf length mus tbe 3, not "
					+ aitrf.length);
		}
		if (rgcrf.length != 3) {
			throw new IllegalArgumentException("rgcrf length mus tbe 3, not "
					+ rgcrf.length);
		}
		if (vgcrf.length != 3) {
			throw new IllegalArgumentException("vgcrf length mus tbe 3, not "
					+ vgcrf.length);
		}
		if (agcrf.length != 3) {
			throw new IllegalArgumentException("agcrf length mus tbe 3, not "
					+ agcrf.length);
		}
		// check trans, these also throw IllegalArgumentException
		trans.checkColumnDimension(3);
		trans.checkRowDimension(3);

		nutationOut no;

		Matrix st = new Matrix(3, 3);
		Matrix stdot = new Matrix(3, 3);
		Matrix pm = new Matrix(3, 3);
		double thetasa;

		// ---- find matrices
		precessOut po = IAU2000.precess(ttt, eOpt.e80);

		// if (nutopt) {
		no = IAU2000.nutation(ttt, ddpsi, ddeps, iau80rec);
		// } else {
		// // TODO What is this branch for? MATLAB doesn't have it.
		// meaneps = ((0.001813 * ttt - 0.00059) * ttt - 46.8150) * ttt
		// + 84381.448;
		// meaneps = MathUtils.mod(meaneps / 3600.0, 360.0);
		// meaneps = meaneps * MathUtils.DEG2RAD;
		// trueeps = meaneps + deltaeps;
		//
		// double cospsi = Math.cos(deltapsi);
		// double sinpsi = Math.sin(deltapsi);
		// double coseps = Math.cos(meaneps);
		// double sineps = Math.sin(meaneps);
		// double costrueeps = Math.cos(trueeps);
		// double sintrueeps = Math.sin(trueeps);
		//
		// nut.set(0, 0, cospsi);
		// nut.set(0, 1, costrueeps * sinpsi);
		// nut.set(0, 2, sintrueeps * sinpsi);
		// nut.set(1, 0, -coseps * sinpsi);
		// nut.set(1, 1, costrueeps * coseps * cospsi + sintrueeps * sineps);
		// nut.set(1, 2, sintrueeps * coseps * cospsi - sineps * costrueeps);
		// nut.set(2, 0, -sineps * sinpsi);
		// nut.set(2, 1, costrueeps * sineps * cospsi - sintrueeps * coseps);
		// nut.set(2, 2, sintrueeps * sineps * cospsi + costrueeps * coseps);
		// }

		IAU2000.sidereal(jdut1, no.deltapsi, no.meaneps, no.omega, lod,
				eqeterms, st, stdot);
		pm = IAU2000.polarm(xp, yp, ttt, eOpt.e80);

		// ---- perform transformations
		thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);

		VectorN omega_earth = new VectorN(0.0, 0.0, thetasa);

		if (itrf2gcrf) {
			// ecef2eci.m:
			// MATLAB: rpef = pm*recef;
			// MATLAB: reci = prec*nut*st*rpef;
			// CPP : matvecmult(pm, ritrf, rpef);
			// CPP : matmult( prec, nut, temp);
			// CPP : matmult( temp, st, trans);
			// CPP : matvecmult( trans, rpef, rgcrf);
			VectorN rpef = pm.times(ritrf);
			Matrix temp = po.prec.times(no.nut);
			trans = temp.times(st); // output
			rgcrf = trans.times(rpef); // output

			// MATLAB: vpef = pm*vecef;
			// MATLAB: veci = prec*nut*st*(vpef + cross(omegaearth,rpef));
			// CPP : matvecmult(pm, vitrf, vpef);
			// CPP : cross( omegaearth, rpef, omgxr);
			// CPP : addvec( 1.0, vpef, 1.0, omgxr, tempvec1);
			// CPP : matvecmult( trans, tempvec1, vgcrf);
			VectorN vpef = pm.times(vitrf);
			VectorN omgxr = omega_earth.crossProduct(rpef);
			VectorN tempvec1 = vpef.plus(omgxr);
			vgcrf = trans.times(tempvec1); // output

			// MATLAB: temp = cross(omegaearth,rpef);
			// MATLAB: aeci = prec*nut*st*( pm*aecef + cross(omegaearth,temp)
			// ...
			// MATLAB: + 2.0*cross(omegaearth,vpef) );
			// CPP : matvecmult(pm, aitrf, apef);
			// CPP : cross(omegaearth,omgxr, omgxomgxr);
			// CPP : cross( omegaearth, vpef, omgxv);
			// CPP : addvec( 1.0, apef, 1.0, omgxomgxr, tempvec);
			// CPP : addvec( 1.0, tempvec, 2.0, omgxv, tempvec1);
			// CPP : matvecmult( trans, tempvec1, agcrf);
			VectorN apef = pm.times(aitrf);
			VectorN omgxomgxr = omega_earth.crossProduct(omgxr);
			VectorN omgxv = omega_earth.crossProduct(vpef);
			VectorN tempvec = apef.plus(omgxomgxr);
			tempvec1 = tempvec.plus(omgxv);
			agcrf = trans.times(tempvec1); // output
		} else {
			// eci2ecef.m:
			// MATLAB: rpef = st'*nut'*prec'*reci;
			// MATLAB: recef = pm'*rpef;
			// CPP : mattrans(pm, pmp);
			// CPP : mattrans(st, stp);
			// CPP : mattrans(nut, nutp);
			// CPP : mattrans(prec, precp);
			// CPP : matmult( stp, nutp, temp);
			// CPP : matmult( temp, precp, trans);
			// CPP : matvecmult(trans, rgcrf, rpef);
			// CPP : matvecmult(pmp, rpef, ritrf);
			trans = st.transpose().times(
					no.nut.transpose().times(po.prec.transpose())); // output
			Matrix pmp = pm.transpose();
			VectorN rpef = trans.times(rgcrf);
			ritrf = pmp.times(rpef); // output

			// MATLAB: vpef = st'*nut'*prec'*veci - cross( omegaearth,rpef );
			// MATLAB: vecef = pm'*vpef;
			// CPP : cross( omegaearth, rpef, omgxr);
			// CPP : matvecmult(trans, vgcrf, tempvec1);
			// CPP : addvec( 1.0, tempvec1, -1.0, omgxr, vpef);
			// CPP : matvecmult( pmp, vpef, vitrf);
			VectorN omgxr = omega_earth.crossProduct(rpef);
			VectorN tempvec1 = trans.times(vgcrf);
			VectorN vpef = tempvec1.minus(omgxr);
			vitrf = pmp.times(vpef); // output

			// MATLAB: temp = cross(omegaearth,rpef);
			// MATLAB: aecef = pm'*(st'*nut'*prec'*aeci - cross(omegaearth,temp)
			// ...
			// MATLAB: - 2.0*cross(omegaearth,vpef));
			// CPP : addvec( 1.0, tempvec1, -1.0, omgxr, vpef);
			// CPP : cross( omegaearth, vpef, omgxv);
			// CPP : cross(omegaearth,omgxr, omgxomgxr);
			// CPP : matvecmult( trans, agcrf, tempvec1);
			// CPP : addvec( 1.0, tempvec1, -1.0, omgxomgxr, tempvec);
			// CPP : addvec( 1.0, tempvec, -2.0, omgxv, apef);
			// CPP : matvecmult( pmp, apef, aitrf);
			vpef = tempvec1.minus(omgxr);
			VectorN omgxv = omega_earth.crossProduct(vpef);
			VectorN omgxomgxr = omega_earth.crossProduct(omgxr);
			tempvec1 = trans.times(agcrf);
			VectorN tempvec = tempvec1.minus(omgxomgxr);
			VectorN apef = tempvec.minus(omgxv.times(2.0));
			aitrf = pmp.times(apef); // output
		}

	}

	/**
	 * This function calulates the transformation matrix that accounts for the 
	 * effects of preseccion. both the 1980 and 2000 theories are handled.  
	 * Note that the output content differs depending on the opt parameter.
	 * 
	 * @param ttt Julian centuries of TT
	 * @param opt The method option
	 * @return A precessOut class with the results.  Note that not every member in the 
	 * returned precessOut will be set.  Check its "type" member before accessing the 
	 * returned data.
	 */
	public static precessOut precess(final double ttt, final eOpt opt) {
		/*
		 * from ASTREDUC.CPP line 795.
		 * (removed from argument list because Java passes by value double&
		 * psia, double& wa, double& epsa, double& chia, double prec[3][3]
		 */
		/* -----------------------------------------------------------------------------
		 *
		 *                           function precess
		 *
		 *  this function calulates the transformation matrix that accounts for the effects
		 *    of precession. 
		 *
		 *  author        : david vallado                  719-573-2600   25 jun 2002
		 *
		 *  revisions
		 *    vallado     - conversion to c++                             21 feb 2005
		 *    vallado     - misc updates, nomenclature, etc               23 nov 2005
		 *
		 *  inputs          description                    range / units
		 *    ttt         - julian centuries of tt
		 *    opt         - method option                  e96, e80
		 *
		 *  outputs       :
		 *    prec        - transformation matrix for mod - j2000 (80 only)
		 *    psia        - cannonical precession angle    rad    (00 only)
		 *    wa          - cannonical precession angle    rad    (00 only)
		 *    epsa        - cannonical precession angle    rad    (00 only)
		 *    chia        - cannonical precession angle    rad    (00 only)
		 *    prec        - matrix converting from "mod" to gcrf
		 *
		 *  locals        :
		 *    zeta        - precession angle               rad
		 *    z           - precession angle               rad
		 *    theta       - precession angle               rad
		 *    oblo        - obliquity value at j2000 epoch "//
		 *
		 *  coupling      :
		 *    none        -
		 *
		 *  references    :
		 *    vallado       2007, 217, 228
		 * --------------------------------------------------------------------------- */
		precessOut po = new precessOut();

		double zeta = 0.0, theta = 0.0, z = 0.0, coszeta, sinzeta, costheta, sintheta, cosz, sinz, oblo = 0.0;
		/*
		 * converted to Matrices: double p1[3][3], p2[3][3], p3[3][3], p4[3][3],
		 * tr1[3][3], tr2[3][3];
		 */
		// Matrix p1 = new Matrix(3,3);
		// Matrix p2 = new Matrix(3,3);
		// Matrix p3 = new Matrix(3,3);
		// Matrix p4 = new Matrix(3,3);
		// Matrix tr1 = new Matrix(3,3);
		// Matrix tr2 = new Matrix(3,3);
		// removed for Java conversion
		// char iauhelp;
		// sethelp(iauhelp, ' ');
		double convrt = MathUtils.PI / (180.0 * 3600.0);

		switch (opt) {

		// Note, the C++ code uses this calculation for these two cases but it
		// doesn't handle the e00* cases.
		
		case e80:
		case e96:

			// The C++ and MATLAB code agree on these coefficients even though the
			// code is rearranged.
			
	        // ------------------- iau 77 precession angles --------------------
			oblo = 84381.448;
			po.psia = ((-0.001147 * ttt - 1.07259) * ttt + 5038.7784) * ttt;
			po.wa = ((-0.007726 * ttt + 0.05127) * ttt) + oblo;
			po.epsa = ((0.001813 * ttt - 0.00059) * ttt - 46.8150) * ttt + oblo;
			po.chia = ((-0.001125 * ttt - 2.38064) * ttt + 10.5526) * ttt;

			zeta = ((0.017998 * ttt + 0.30188) * ttt + 2306.2181) * ttt;
			theta = ((-0.041833 * ttt - 0.42665) * ttt + 2004.3109) * ttt;
			z = ((0.018203 * ttt + 1.09468) * ttt + 2306.2181) * ttt;
			break;
			
		// Note, the MATLAB code only handled e80 and "00".  These cases were taken
		// from the MATLAB code.  The MATLAB comments were:
		//  ------------------ iau 00 precession angles -------------------
		case e00a:
		case e00b:
	        double ttt2= ttt * ttt;
	        double ttt3= ttt2 * ttt;
            double ttt4 = ttt2 * ttt2;
            double ttt5 = ttt2 * ttt3;
            po.psia =             5038.47875*ttt - 1.07259*ttt2 - 0.001147*ttt3;
            po.wa   = 84381.448 -    0.02524*ttt + 0.05127*ttt2 - 0.007726*ttt3;
            // MATLAB had this as ea, changed to epsa:
            po.epsa   = 84381.448 -   46.84024*ttt - 0.00059*ttt2 + 0.001813*ttt3;
            // MATLAB had this as xa, changed to chia:
            po.chia   =               10.5526*ttt  - 2.38064*ttt2 - 0.001125*ttt3;
            zeta = 2.5976176 + 2306.0809506*ttt  + 0.3019015 *ttt2 + 0.0179663*ttt3 
                                - 0.0000327*ttt4 - 0.0000002*ttt5;
            theta=             2004.1917476*ttt  - 0.4269353*ttt2 - 0.0418251*ttt3 
                                - 0.0000601*ttt4 - 0.0000001*ttt5;
            z    = 2.5976176 + 2306.0803226*ttt  + 1.0947790*ttt2 + 0.0182273*ttt3 
                                + 0.0000470*ttt4 - 0.0000003*ttt5;
			break;

		default:
			throw new java.lang.IllegalArgumentException("Unexpected IAU designation in fundarg(): "
					+ opt.name());

		}

		// convert units to rad
		po.psia = po.psia * convrt;
		po.wa = po.wa * convrt;
		oblo = oblo * convrt;
		po.epsa = po.epsa * convrt;
		po.chia = po.chia * convrt;

		zeta = zeta * convrt;
		theta = theta * convrt;
		z = z * convrt;

		// Set the po.prec matrix only if using e80 or e96:
		if ((opt == eOpt.e80) | (opt == eOpt.e96)) {
			coszeta = Math.cos(zeta);
			sinzeta = Math.sin(zeta);
			costheta = Math.cos(theta);
			sintheta = Math.sin(theta);
			cosz = Math.cos(z);
			sinz = Math.sin(z);

			// ----------------- form matrix mod to gcrf ------------------
			po.prec.set(0, 0, coszeta * costheta * cosz - sinzeta * sinz);
			po.prec.set(0, 1, coszeta * costheta * sinz + sinzeta * cosz);
			po.prec.set(0, 2, coszeta * sintheta);
			po.prec.set(1, 0, -sinzeta * costheta * cosz - coszeta * sinz);
			po.prec.set(1, 1, -sinzeta * costheta * sinz + coszeta * cosz);
			po.prec.set(1, 2, -sinzeta * sintheta);
			po.prec.set(2, 0, -sintheta * cosz);
			po.prec.set(2, 1, -sintheta * sinz);
			po.prec.set(2, 2, costheta);
			
			po.type = opt; // always note which type drove the calculations

			// ----------------- do rotations instead ----------------------
			// p1 = rot3mat( z );
			// p2 = rot2mat( -theta );
			// p3 = rot3mat( zeta );
			// prec = p3*p2*p1;
		}

		return po;
	}

	/**
	 * This function calculates the transformation matrix that accounts for the
	 * effects of nutation.
	 * 
	 * @param ttt
	 *            Julian centuries of TT
	 * @param ddpsi
	 *            delta psi correction to gcrf (rad)
	 * @param ddeps
	 *            delta eps correction to gcrf (rad)
	 * @param iau80rec
	 *            record containing IAU80 constants (rad)
	 * @return
	 */
	public static nutationOut nutation(final double ttt,
			final double ddpsi, final double ddeps, final iau80data iau80rec)
	/*
	 * From ASTREDUC.CPP, line 907.
	 * Modifications: (removed C++ outputs from argument list since Java passes
	 * by value) double& deltapsi, double& deltaeps, double& trueeps, double&
	 * meaneps, double& omega, double nut[3][3] Also removed nutopt because only
	 * one option was implemented in the C++ code. There was no nutopt in the
	 * MATLAB code.
	 */

	/* -----------------------------------------------------------------------------
	 *
	 *                           function nutation
	 *
	 *  this function calulates the transformation matrix that accounts for the
	 *    effects of nutation.
	 *
	 *  author        : david vallado                  719-573-2600   27 jun 2002
	 *
	 *  revisions
	 *    vallado     - consolidate with iau 2000                     14 feb 2005
	 *    vallado     - conversion to c++                             21 feb 2005
	 *
	 *  inputs          description                    range / units
	 *    ttt         - julian centuries of tt
	 *    ddpsi       - delta psi correction to gcrf   rad
	 *    ddeps       - delta eps correction to gcrf   rad
	 *    nutopt      - nutation option                calc 'c', read 'r'
	 *    iau80rec    - record containing the iau80 constants rad
	 *
	 *  outputs       :
	 *    deltapsi    - nutation in longiotude angle   rad
	 *    trueeps     - true obliquity of the ecliptic rad
	 *    meaneps     - mean obliquity of the ecliptic rad
	 *    omega       -                                rad
	 *    nut         - transform matrix for tod -     mod
	 *
	 *  locals        :
	 *    iar80       - integers for fk5 1980
	 *    rar80       - reals for fk5 1980
	 *    l           -                                rad
	 *    ll          -                                rad
	 *    f           -                                rad
	 *    d           -                                rad
	 *    deltaeps    - change in obliquity            rad
	 *
	 *  coupling      :
	 *    fundarg     - find fundamental arguments
	 *    fmod      - modulus division
	 *
	 *  references    :
	 *    vallado       2007, 217, 228
	 * --------------------------------------------------------------------------- */
	{
		// dropped from java conversion in favor of fundargOut fo, below
		// double MathUtils.DEG2RAD, l, l1, f, d,
		// lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep,
		// precrate,
		fundargOut fo;

		// this method's Java output:
		nutationOut no = new nutationOut();

		double cospsi, sinpsi, coseps, sineps, costrueeps, sintrueeps;
		int i;
		double tempval;

		// (removed for Java conversion)
		// char iauhelp;
		// sethelp(iauhelp, ' ');

		// ---- determine coefficients for iau 1980 nutation theory ----
		no.meaneps = ((0.001813 * ttt - 0.00059) * ttt - 46.8150) * ttt
				+ 84381.448;
		no.meaneps = MathUtils.mod(no.meaneps / 3600.0, 360.0);
		no.meaneps = no.meaneps * MathUtils.DEG2RAD;

		// if ( nutopt =='c' )
		// {
		fo = IAU2000.fundarg(ttt, eOpt.e80);
		/*
		 * removed in favor of the fo instance: l, l1, f, d, omega, lonmer,
		 * lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate );
		 */

		no.deltapsi = 0.0;
		no.deltaeps = 0.0;
		for (i = 106; i >= 1; i--) {
			tempval = iau80rec.iar80[i][1] * fo.l + iau80rec.iar80[i][2]
					* fo.l1 + iau80rec.iar80[i][3] * fo.f
					+ iau80rec.iar80[i][4] * fo.d + iau80rec.iar80[i][5]
					* fo.omega;
			no.deltapsi = no.deltapsi
					+ (iau80rec.rar80[i][1] + iau80rec.rar80[i][2] * ttt)
					* Math.sin(tempval);
			no.deltaeps = no.deltaeps
					+ (iau80rec.rar80[i][3] + iau80rec.rar80[i][4] * ttt)
					* Math.cos(tempval);
		}

		// --------------- find nutation parameters --------------------
		no.deltapsi = MathUtils.mod(no.deltapsi + ddpsi / MathUtils.DEG2RAD,
				360.0)
				* MathUtils.DEG2RAD;
		no.deltaeps = MathUtils.mod(no.deltaeps + ddeps / MathUtils.DEG2RAD,
				360.0)
				* MathUtils.DEG2RAD;
		// }

		no.trueeps = no.meaneps + no.deltaeps;

		cospsi = Math.cos(no.deltapsi);
		sinpsi = Math.sin(no.deltapsi);
		coseps = Math.cos(no.meaneps);
		sineps = Math.sin(no.meaneps);
		costrueeps = Math.cos(no.trueeps);
		sintrueeps = Math.sin(no.trueeps);

		no.nut.set(0, 0, cospsi);
		no.nut.set(0, 1, costrueeps * sinpsi);
		no.nut.set(0, 2, sintrueeps * sinpsi);
		no.nut.set(1, 0, -coseps * sinpsi);
		no.nut.set(1, 1, costrueeps * coseps * cospsi + sintrueeps * sineps);
		no.nut.set(1, 2, sintrueeps * coseps * cospsi - sineps * costrueeps);
		no.nut.set(2, 0, -sineps * sinpsi);
		no.nut.set(2, 1, costrueeps * sineps * cospsi - sintrueeps * coseps);
		no.nut.set(2, 2, sintrueeps * sineps * cospsi + costrueeps * coseps);

		// n1 = rot1mat( trueeps );
		// n2 = rot3mat( deltapsi );
		// n3 = rot1mat( -meaneps );
		// nut = n3*n2*n1;

		return no;
	}

	/**
	 * This function calulates the transformation matrix that accounts for the
	 * effects of sidereal time. Notice that deltaspi should not be moded to a
	 * positive number because it is multiplied rather than used in a
	 * trigonometric argument.
	 * 
	 * @param jdut1
	 *            julian centuries of ut1 (days)
	 * @param deltapsi
	 *            nutation angle (rad)
	 * @param meaneps
	 *            mean obliquity of the ecliptic (rad)
	 * @param omega
	 *            long of asc node of moon (rad)
	 * @param lod
	 *            length of day (sec)
	 * @param eqeterms
	 *            terms for ast calculation 0,2
	 * @param st
	 *            (output) transformation matrix for pef - tod
	 * @param stdot
	 *            (output) transformation matrix for pef - tod rate
	 */
	public static void sidereal(double jdut1, double deltapsi,
			double meaneps, double omega, double lod, int eqeterms, Matrix st,
			Matrix stdot) {
		/*  (original comments from ASTREDUC.CPP)
		 * -----------------------------------------------------------------------------
		 *
		 *                           function sidereal
		 *
		 *  this function calulates the transformation matrix that accounts for the
		 *    effects of sidereal time. Notice that deltaspi should not be moded to a
		 *    positive number because it is multiplied rather than used in a
		 *    trigonometric argument.
		 *
		 *  author        : david vallado                  719-573-2600   25 jun 2002
		 *
		 *  revisions
		 *    vallado     - fix units on kinematic terms                   5 sep 2002
		 *    vallado     - add terms                                     30 sep 2002
		 *    vallado     - consolidate with iau 2000                     14 feb 2005
		 *    vallado     - conversion to c++                             21 feb 2005
		 *
		 *  inputs          description                    range / units
		 *    jdut1       - julian centuries of ut1        days
		 *    deltapsi    - nutation angle                 rad
		 *    meaneps     - mean obliquity of the ecliptic rad
		 *    omega       - long of asc node of moon       rad
		 *    lod         - length of day                  sec
		 *    eqeterms    - terms for ast calculation      0,2
		 *
		 *  outputs       :
		 *    st          - transformation matrix for pef - tod
		 *    stdot       - transformation matrix for pef - tod rate
		 *
		 *  locals        :
		 *    gmst         - mean greenwich sidereal time   0 to 2pi rad
		 *    ast         - apparent gmst                   0 to 2pi rad
		 *    hr          - hour                           hr
		 *    min         - minutes                        min
		 *    sec         - seconds                        sec
		 *    temp        - temporary vector
		 *    tempval     - temporary variable
		 *
		 *  coupling      :
		 *
		 *  references    :
		 *    vallado       2007, 217, 228
		 * --------------------------------------------------------------------------- */
		double gmst, ast, thetasa, omegaearth;
		// char iauhelp; (removed for Java conversion)

		// sethelp(iauhelp, ' '); // dropped for Java conversion

		// ------------------------ find gmst --------------------------
		gmst = gstime(jdut1);

		// ------------------------ find mean ast ----------------------
		if ((jdut1 > 2450449.5) && (eqeterms > 0)) {
			ast = gmst + deltapsi * Math.cos(meaneps) + 0.00264 * MathUtils.PI
					/ (3600 * 180) * Math.sin(omega) + 0.000063 * MathUtils.PI
					/ (3600 * 180) * Math.sin(2.0 * omega);
		} else
			ast = gmst + deltapsi * Math.cos(meaneps);

		// ast = fmod (ast,2.0*MathUtils.PI);
		ast = MathUtils.mod(ast, 2.0 * MathUtils.PI);

		thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
		omegaearth = thetasa;

		// Java conversino: directly set the elements in the given matrix
		st.set(0, 0, Math.cos(ast));
		st.set(0, 1, -Math.sin(ast));
		st.set(0, 2, 0.0);
		st.set(1, 0, Math.sin(ast));
		st.set(1, 1, Math.cos(ast));
		st.set(1, 2, 0.0);
		st.set(2, 0, 0.0);
		st.set(2, 1, 0.0);
		st.set(2, 2, 1.0);

		// compute sidereal time rate matrix
		// Java conversino: directly set the elements in the given matrix
		stdot.set(0, 0, -omegaearth * Math.sin(ast));
		stdot.set(0, 1, -omegaearth * Math.cos(ast));
		stdot.set(0, 2, 0.0);
		stdot.set(1, 0, omegaearth * Math.cos(ast));
		stdot.set(1, 1, -omegaearth * Math.sin(ast));
		stdot.set(1, 2, 0.0);
		stdot.set(2, 0, 0.0);
		stdot.set(2, 1, 0.0);
		stdot.set(2, 2, 0.0);

	}

	/**
	 * Calculates the transformation matrix that accounts for polar motion. Both
	 * the 1980 and 2000 theories are handled. Note that the rotation order is
	 * different between 1980 and 2000. 
	 * 
	 * @param xp
	 *            ploar motion coefficient (arcsec)
	 * @param yp
	 *            ploar motion coefficient (arcsec)
	 * @param ttt
	 *            julian centuries of tt (00 theory only)
	 * @param opt
	 *            method option (e00a, e00b, e96, e80)
	 * @return transformation matrix for itrf - pef
	 */
	public static Matrix polarm(final double xp, final double yp,
			final double ttt, final eOpt opt) {
		/* (original comments from ASTREDUC.CPP, line 1123)
		 * -----------------------------------------------------------------------------
		 *
		 *                           function polarm
		 *
		 *  this function calulates the transformation matrix that accounts for polar
		 *    motion. both the 1980 and 2000 theories are handled. note that the rotation
		 *    order is different between 1980 and 2000 .
		 *
		 *  author        : david vallado                  719-573-2600   25 jun 2002
		 *
		 *  revisions
		 *    vallado     - conversion to c++                             23 nov 2005
		 *
		 *  inputs          description                    range / units
		 *    xp          - polar motion coefficient       rad
		 *    yp          - polar motion coefficient       rad
		 *    ttt         - julian centuries of tt (00 theory only)
		 *    opt         - method option                  e00a, e00b, e96, e80
		 *
		 *  outputs       :
		 *    pm          - transformation matrix for itrf - pef
		 *
		 *  locals        :
		 *    convrt      - conversion from arcsec to rad
		 *    sp          - s prime value
		 *
		 *  coupling      :
		 *    none.
		 *
		 *  references    :
		 *    vallado       2007, 217, 228
		 * --------------------------------------------------------------------------- */
		double convrt, cosxp, cosyp, sinxp, sinyp, sp, cossp, sinsp; 

		/*
		 * Note, moved this output from the argument list for Java conversion
		 * since Java is pass-by-value.
		 */
		double[][] pm = new double[3][3];

		/*
		 * NOTE, the C++ header notes the xp and yp inputs as radians but the 
		 * matlab polarm.m has them as arcsec.  This conversion confirms that
		 * the inputs are arcsec.
		 */
		convrt = MathUtils.PI / (180.0 * 3600.0); // arcsec to rad

		cosxp = Math.cos(xp);
		sinxp = Math.sin(xp);
		cosyp = Math.cos(yp);
		sinyp = Math.sin(yp);

		if ((opt == IAU2000.eOpt.e80) | (opt == IAU2000.eOpt.e96)) {
			pm[0][0] = cosxp;
			pm[0][1] = 0.0;
			pm[0][2] = -sinxp;
			pm[1][0] = sinxp * sinyp;
			pm[1][1] = cosyp;
			pm[1][2] = cosxp * sinyp;
			pm[2][0] = sinxp * cosyp;
			pm[2][1] = -sinyp;
			pm[2][2] = cosxp * cosyp;

			// a1 = rot2mat(xp);
			// a2 = rot1mat(yp);
			// pm = a2*a1;
		} else {
			// approximate sp value in rad
			sp = -47.0e-6 * ttt * convrt;
			cossp = Math.cos(sp);
			sinsp = Math.sin(sp);

			// form the matrix
			pm[0][0] = cosxp * cossp;
			pm[0][1] = -cosyp * sinsp + sinyp * sinxp * cossp;
			pm[0][2] = -sinyp * sinsp - cosyp * sinxp * cossp;
			pm[1][0] = cosxp * sinsp;
			pm[1][1] = cosyp * cossp + sinyp * sinxp * sinsp;
			pm[1][2] = sinyp * cossp - cosyp * sinxp * sinsp;
			pm[2][0] = sinxp;
			pm[2][1] = -sinyp * cosxp;
			pm[2][2] = cosyp * cosxp;

			// a1 = rot1mat(yp);
			// a2 = rot2mat(xp);
			// a3 = rot3mat(-sp);
			// pm = a3*a2*a1;
		}

		return new Matrix(pm);
	}

	/**
	 * This function calulates the delaunay variables and planetary values for
	 * several theories. <br>
	 * <i>(Started using code from ASTREDUC.CPP line 695, then adapted the
	 * MATLAB version of the code, fundarg.m, because all of the systems were
	 * implemented, unlike the C++ code.)</i>
	 * 
	 * @param ttt
	 *            julian centuries of tt
	 * @param opt
	 *            method option
	 * @return fundargOut instance that contains the calculated delaunay
	 *         variables and planetary longitudes
	 */
	public static fundargOut fundarg(final double ttt,
			final IAU2000.eOpt opt)
	/*
	 * Note: Java is pass-by-value so these C++ references were moved into a
	 * returned fundargOut instance), double l, double l1, double f, double d,
	 * double omega, double lonmer, double lonven, double lonear, double lonmar,
	 * double lonjup, double lonsat, double lonurn, double lonnep, double
	 * precrate)
	 */
	{
		/* (original comments from ASTREDUC.CPP)
		 * -----------------------------------------------------------------------------
		 *
		 *                           function fundarg
		 *
		 *  this function calulates the delauany variables and planetary values for
		 *  several theories.
		 *
		 *  author        : david vallado                  719-573-2600   16 jul 2004
		 *
		 *  revisions
		 *    vallado     - conversion to c++                             23 nov 2005
		 *
		 *  inputs          description                    range / units
		 *    ttt         - julian centuries of tt
		 *    opt         - method option                  e96, e80
		 *
		 *  outputs       :
		 *    l           - delaunay element               rad
		 *    l1          - delaunay element               rad
		 *    f           - delaunay element               rad
		 *    d           - delaunay element               rad
		 *    omega       - delaunay element               rad
		 *    planetary longitudes                         rad
		 *
		 *  locals        :
		 *
		 *
		 *  coupling      :
		 *    none        -
		 *
		 *  references    :
		 *    vallado       2007, 217
		 * --------------------------------------------------------------------------- */

		// double deg2rad; // dropped for Java conversion
		// char iauhelp; // dropped for Java conversion
		// sethelp(iauhelp, ' '); // dropped for Java conversion
		// deg2rad = pi/180.0; // dropped for Java conversion
		fundargOut fo = new fundargOut();

		// (begin section from MATLAB fundarg.m)

		// ---- determine coefficients for iau 2000 nutation theory ----
		double ttt2 = ttt * ttt; // ttt squared
		double ttt3 = ttt2 * ttt; // ttt cubed
		double ttt4 = ttt2 * ttt2; // ttt^4

		switch (opt) {

		// ---- iau 2000a theory
		case e00a:

			// ------ form the delaunay fundamental arguments in deg
			fo.l = 134.96340251 + (1717915923.2178 * ttt + 31.8792 * ttt2
					+ 0.051635 * ttt3 - 0.00024470 * ttt4) / 3600.0;
			fo.l1 = 357.52910918 + (129596581.0481 * ttt - 0.5532 * ttt2
					- 0.000136 * ttt3 - 0.00001149 * ttt4) / 3600.0;
			fo.f = 93.27209062 + (1739527262.8478 * ttt - 12.7512 * ttt2
					+ 0.001037 * ttt3 + 0.00000417 * ttt4) / 3600.0;
			fo.d = 297.85019547 + (1602961601.2090 * ttt - 6.3706 * ttt2
					+ 0.006593 * ttt3 - 0.00003169 * ttt4) / 3600.0;
			fo.omega = 125.04455501 + (-6962890.5431 * ttt + 7.4722 * ttt2
					+ 0.007702 * ttt3 - 0.00005939 * ttt4) / 3600.0;

			// ------ form the planetary arguments in deg
			fo.planetlon[0] = 252.250905494 + 149472.6746358 * ttt;
			fo.planetlon[1] = 181.979800853 + 58517.8156748 * ttt;
			fo.planetlon[2] = 100.466448494 + 35999.3728521 * ttt;
			fo.planetlon[3] = 355.433274605 + 19140.299314 * ttt;
			fo.planetlon[4] = 34.351483900 + 3034.90567464 * ttt;
			fo.planetlon[5] = 50.0774713998 + 1222.11379404 * ttt;
			fo.planetlon[6] = 314.055005137 + 428.466998313 * ttt;
			fo.planetlon[7] = 304.348665499 + 218.486200208 * ttt;
			fo.precrate = 1.39697137214 * ttt + 0.0003086 * ttt2;
			break;

		// ---- iau 2000b theory
		case e00b:

			// ------ form the delaunay fundamental arguments in deg
			fo.l = 134.96340251 + (1717915923.2178 * ttt) / 3600.0;
			fo.l1 = 357.52910918 + (129596581.0481 * ttt) / 3600.0;
			fo.f = 93.27209062 + (1739527262.8478 * ttt) / 3600.0;
			fo.d = 297.85019547 + (1602961601.2090 * ttt) / 3600.0;
			fo.omega = 125.04455501 + (-6962890.5431 * ttt) / 3600.0;

			// ------ form the planetary arguments in deg
			fo.planetlon[0] = 0.0;
			fo.planetlon[1] = 0.0;
			fo.planetlon[2] = 0.0;
			fo.planetlon[3] = 0.0;
			fo.planetlon[4] = 0.0;
			fo.planetlon[5] = 0.0;
			fo.planetlon[6] = 0.0;
			fo.planetlon[7] = 0.0;
			fo.precrate = 0.0;
			break;

		// ---- iau 1996 theory
		case e96:

			fo.l = 134.96340251 + (1717915923.2178 * ttt + 31.8792 * ttt2
					+ 0.051635 * ttt3 - 0.00024470 * ttt4) / 3600.0;
			fo.l1 = 357.52910918 + (129596581.0481 * ttt - 0.5532 * ttt2
					- 0.000136 * ttt3 - 0.00001149 * ttt4) / 3600.0;
			fo.f = 93.27209062 + (1739527262.8478 * ttt - 12.7512 * ttt2
					+ 0.001037 * ttt3 + 0.00000417 * ttt4) / 3600.0;
			fo.d = 297.85019547 + (1602961601.2090 * ttt - 6.3706 * ttt2
					+ 0.006593 * ttt3 - 0.00003169 * ttt4) / 3600.0;
			fo.omega = 125.04455501 + (-6962890.2665 * ttt + 7.4722 * ttt2
					+ 0.007702 * ttt3 - 0.00005939 * ttt4) / 3600.0;
			// ------ form the planetary arguments in deg
			fo.planetlon[0] = 0.0;
			fo.planetlon[1] = 181.979800853 + 58517.8156748 * ttt; // deg
			fo.planetlon[2] = 100.466448494 + 35999.3728521 * ttt;
			fo.planetlon[3] = 355.433274605 + 19140.299314 * ttt;
			fo.planetlon[4] = 34.351483900 + 3034.90567464 * ttt;
			fo.planetlon[5] = 50.0774713998 + 1222.11379404 * ttt;
			fo.planetlon[6] = 0.0;
			fo.planetlon[7] = 0.0;
			fo.precrate = 1.39697137214 * ttt + 0.0003086 * ttt2;
			break;

		// ---- iau 1980 theory
		case e80:

			fo.l = 134.96298139 + (1717915922.6330 * ttt + 31.310 * ttt2 + 0.064 * ttt3) / 3600.0;
			fo.l1 = 357.52772333 + (129596581.2240 * ttt - 0.577 * ttt2 - 0.012 * ttt3) / 3600.0;
			fo.f = 93.27191028 + (1739527263.1370 * ttt - 13.257 * ttt2 + 0.011 * ttt3) / 3600.0;
			fo.d = 297.85036306 + (1602961601.3280 * ttt - 6.891 * ttt2 + 0.019 * ttt3) / 3600.0;
			fo.omega = 125.04452222 + (-6962890.5390 * ttt + 7.455 * ttt2 + 0.008 * ttt3) / 3600.0;
			// ------ form the planetary arguments in deg
			fo.planetlon[0] = 252.3 + 149472.0 * ttt;
			fo.planetlon[1] = 179.9 + 58517.8 * ttt; // deg
			fo.planetlon[2] = 98.4 + 35999.4 * ttt;
			fo.planetlon[3] = 353.3 + 19140.3 * ttt;
			fo.planetlon[4] = 32.3 + 3034.9 * ttt;
			fo.planetlon[5] = 48.0 + 1222.1 * ttt;
			fo.planetlon[6] = 0.0;
			fo.planetlon[7] = 0.0;
			fo.precrate = 0.0;
			break;

		default:
			throw new java.lang.IllegalArgumentException("Unexpected IAU designation in fundarg(): "
					+ opt.name());

		} // switch

		// (end section from MATLAB fundarg.m)

		// ---- convert units to rad
		fo.l = MathUtils.mod(fo.l, 360.0) * MathUtils.DEG2RAD;
		fo.l1 = MathUtils.mod(fo.l1, 360.0) * MathUtils.DEG2RAD;
		fo.f = MathUtils.mod(fo.f, 360.0) * MathUtils.DEG2RAD;
		fo.d = MathUtils.mod(fo.d, 360.0) * MathUtils.DEG2RAD;
		fo.omega = MathUtils.mod(fo.omega, 360.0) * MathUtils.DEG2RAD;

		/*
		 * For java conversion, replaced individual planet assignments with for
		 * loop, below.
		 */
		for (int i = 0; i < fundargOut.PLANETLONLENGTH; i++) {
			fo.planetlon[i] = MathUtils.mod(fo.planetlon[i], 360.0)
					* MathUtils.DEG2RAD;
		}

		fo.precrate = MathUtils.mod(fo.precrate, 360.0) * MathUtils.DEG2RAD;

		return (fo);
	}

	/**
	 * This function finds the Greenwich Mean Sidereal Time in radians (IAU-82)
	 * from a UT1 Julian Date.<br>
	 * Reference:<br>
	 * <i>Fundamentals of Astrodynamics and Applicaitons, Third Ed.</i>,
     * Vallado, David A., Microcosm Press and Springer, 2007, Section 3.5.2
     * pages 194-195 (Algorithm 15).
	 * 
	 * @param jdut1
	 *            Julian Date in UT1 (days from 4713 BC)
	 * @return Greenwich Mean Sidereal Time (0 to 2*pi rad)
	 */
	public static double gstime(final double jdut1) {
		/*
		 * (from ASTTIME.CPP ln 380)
		 * --------------------------------------------
		 * ---------------------------------
		 * 
		 * function gstime
		 * 
		 * this function finds the greenwich sidereal time (iau-82).
		 * 
		 * author : david vallado 719-573-2600 1 mar 2001
		 * 
		 * inputs description range / units jdut1 - julian date in ut1 days from
		 * 4713 bc
		 * 
		 * outputs : gstime - greenwich sidereal time 0 to 2pi rad
		 * 
		 * locals : temp - temporary variable for doubles rad tut1 - julian
		 * centuries from the jan 1, 2000 12 h epoch (ut1)
		 * 
		 * coupling : none
		 * 
		 * references : vallado 2007, 193, eq 3-43
		 * 
		 * ----------------------------------------------------------------------
		 * -----
		 */

		// (replaced by Java constructs)
		// const double twopi = 2.0 * pi;
		// const double deg2rad = pi / 180.0;
		
		double temp; // temporary variable for reals (rad)
		double tut1; // julian centuries from the jan 1, 2000 12 h epoch (ut1)

		tut1 = (jdut1 - 2451545.0) / 36525.0;
		temp = -6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1
				+ (876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841; // sec

		// convert to radians, mod for 0 to 2*pi
		// (divide by 240.0 because 1 sec = 1/240 deg)
		temp = IAU2000.mod(temp * MathUtils.DEG2RAD / 240.0, MathUtils.TWOPI);

		// ------------------------ check quadrants ---------------------
		if (temp < 0.0)
			temp += MathUtils.TWOPI;

		return temp;
	}

	/**
	 * This function reads in IAU80 nutation data from the given file.
	 * 
	 * @param filename
	 *            The complete path and name of the file containing the IAU80
	 *            data.
	 * @return A populated iau80data object if successful, NULL if all of the
	 *         constants could not be populated from the filename (premature end
	 *         of data).
	 * @throws FileNotFoundException
	 *             If the filename is not a proper file.
	 * @throws IOException
	 *             If there is a premature end to the file.
	 * @throws NumberFormatException
	 *             If the data cannot be properly parsed.
	 */
	public static iau80data iau80in(final String filename)
			throws FileNotFoundException, IOException, NumberFormatException
	/*
	 * (from ASTREDUC.CPP, line 636)
	 * --------------------------------------------
	 * ---------------------------------
	 * 
	 * function iau80in
	 * 
	 * this function initializes the nutation matricies needed for reduction
	 * calculations. the routine needs the filename of the files as input.
	 * 
	 * author : david vallado 719-573-2600 27 may 2002
	 * 
	 * revisions vallado - conversion to c++ 21 feb 2005
	 * 
	 * inputs description range / units
	 * 
	 * outputs : iau80rec - record containing the iau80 constants rad
	 * 
	 * locals : convrt - conversion factor to degrees i,j - index
	 * 
	 * coupling : none -
	 * 
	 * references :
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 * --------------------------------------------------------------------------
	 * -
	 */
	{
		iau80data iau80rec = new iau80data();

		double convrt;
		int i; // , j;

		// ------------------------ implementation -------------------
		convrt = 0.0001 / 3600.0; // 0.0001" to deg

		// infile = fopen("nut80.dat","r");
		FileReader file = new FileReader(filename);
		BufferedReader buff = new BufferedReader(file);

		// Java note: numbering intentionally starts from 1
		for (i = 1; i <= 106; i++) {

			// ret = fscanf(infile, "%d %d %d %d %d %lf %lf %lf %lf %d \n ",
			//&iau80rec.iar80[i][1],&iau80rec.iar80[i][2],&iau80rec.iar80[i][3],
			// &iau80rec.iar80[i][4],&iau80rec.iar80[i][5],
			//&iau80rec.rar80[i][1],&iau80rec.rar80[i][2],&iau80rec.rar80[i][3],
			// &iau80rec.rar80[i][4], &j);
			// for (j = 1; j <= 4; j++) {
			// iau80rec.rar80[i][j]= iau80rec.rar80[i][j] * convrt;
			// }

			String line = buff.readLine();
			if (line == null) {
				iau80rec = null; // premature end of data in a line (no line at all)
				break;				
			}
			String[] subs = line.trim().split("\\s+"); // split via whitespace

			if (subs.length < 9) {
				iau80rec = null; // premature end of data in a line
				break;
			}

			// Java note: indexing of iar80 and rar80 intentionally starts from
			// 1.
			// A bad conversion in any of these will throw a
			// NumberFormatException exception.
			iau80rec.iar80[i][1] = Integer.valueOf(subs[0]).intValue();
			iau80rec.iar80[i][2] = Integer.valueOf(subs[1]).intValue();
			iau80rec.iar80[i][3] = Integer.valueOf(subs[2]).intValue();
			iau80rec.iar80[i][4] = Integer.valueOf(subs[3]).intValue();
			iau80rec.iar80[i][5] = Integer.valueOf(subs[4]).intValue();
			iau80rec.rar80[i][1] = Double.valueOf(subs[5]).doubleValue()
					* convrt;
			iau80rec.rar80[i][2] = Double.valueOf(subs[6]).doubleValue()
					* convrt;
			iau80rec.rar80[i][3] = Double.valueOf(subs[7]).doubleValue()
					* convrt;
			iau80rec.rar80[i][4] = Double.valueOf(subs[8]).doubleValue()
					* convrt;
		}

		buff.close();
		file.close();

		return iau80rec;
	}
	
	/** 
	 * Modulus after division. <br>
	 * 
     * mod(x,y) is x - n.*y where n = floor(x./y) if y ~= 0.  If y is not an
     * integer and the quotient x./y is within roundoff error of an integer,
     * then n is that integer.  The inputs x and y must be real arrays of the
     * same size, or real scalars.
  * 
     * The statement "x and y are congruent mod m" means mod(x,m) == mod(y,m).
  * 
     * <i>modMATLAB matches MATLAB mod(x,y).</i><br>
	 * <ul>By convention:
     *    <li>mod(x,0) is x.</li>
     *    <li>mod(x,x) is 0.</li>
     *    <li>mod(x,y), for x~=y and y~=0, has the same sign as y.</li>
     *    </ul><br>
  * 
     * Note: rem(x,y), for x~=y and y~=0, has the same sign as x.
     * mod(x,y) and rem(x,y) are equal if x and y have the same sign, but
     * differ by y if x and y have different signs.
     * 
	 * @param x Numerator.
	 * @param y Denominator.
	 * @return modulus of (x/y)
	 */
	public static double mod(double x, double y) {
		
		// corner cases:
		if (y == 0.0) {
			return x;
		}
		if (x == y) {
			return 0.0;
		}
		
		double temp = Math.floor(x / y);
		double out = x - temp * y;
		
		// Convention is that out and y have the same sign:
		if (y * out < 0.0) {
			return -out;
		}
		else {
			return out;
		}
	}

	/**
	 * The remainder after division. <br>
	 * 
     * rem(x,y) is x - n.*y where n = fix(x./y) if y ~= 0.  If y is not an
     * integer and the quotient x./y is within roundoff error of an integer,
     * then n is that integer. The inputs x and y must be real arrays of the
     * same size, or real scalars.
     * 
	 * <i>rem matches MATLAB rem(x,y).</i><br>
	 * <ul>By convention:
     *  <li>rem(x,0) is NaN.</li>
     *  <li>rem(x,x), for x~=0, is 0.0</li>
     *  <li>rem(x,y), for x~=y and y~=0, has the same sign as x.</li>
     *  </ul>
	 * @param x Numerator.
	 * @param y Denominator.
	 * @return floating point remainder of (x/y)
	 */
	public static double rem(final double x, final double y) {
		
		// corner cases:
		if (y == 0.0) {
			return java.lang.Double.NaN; // div by zero
		}
		if (x == y) {
			return 0.0;
		}
		
		double temp = fix(x / y);
		double out = x - temp * y;
		
		// Convention is that out and x have the same sign:
		if (x * out < 0.0) {
			return -out;
		}
		else {
			return out;
		}
	}
	
	/**
	 * Rounds x to the nearest integers towards zero.<br>
	 * Note that this has a different rounding effect above or below zero.
	 * <i>Matches MATLAB fix(x).</i>
	 * @param x double to round.
	 * @return double rounded towards zero
	 */
	public static double fix(final double x) {
		if (x < 0.0) {
			return Math.ceil(x);
		}
		else {
			return Math.floor(x);
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
