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
import jat.matvec.data.*;
import jat.spacetime.*;
import jat.cm.*;
import jat.math.*;

/**
 * <P>
 * The GPS_SV Class models a GPS satellite which is created from the broadcast ephemeris.
 * It provides a method to generate the ECEF position of a given GPS SV at a given MJD.
 * The GPS_Constellation class reads in the broadcast ephemeris and constructs the constellation.
 * Reference: Montenbruck, Appendix A.2.2.
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
public class GPS_SV {

    private int prn = 0;

    private GPSTimeFormat toc;

    private double clockBias = 0.0;

    private double clockDrift = 0.0;

    private double clockDriftRate = 0.0;

    private double Crs = 0.0;

    private double Cuc = 0.0;

    private double Crc = 0.0;

    private double Cus = 0.0;

    private double Cic = 0.0;

    private double Cis = 0.0;

    private GPSTimeFormat toe;

    private double sqrtA = 0.0;

    private double deltaN = 0.0;

    private double ecc = 0.0;

    private double inc = 0.0;

    private double idot = 0.0;

    private double omega = 0.0;

    private double omegadot = 0.0;

    private double argp = 0.0;

    private double ma = 0.0;

    /** Creates a new instance of GPS_SV.
     * @param svid Integer containing the SV PRN.
     */
    public GPS_SV(int svid) {
        this.prn = svid;
    }

    /**
     * Set the TOC
     * @param clock GPSTimeFormat for the TOC
     */
    public void setTOC(GPSTimeFormat clock){
        this.toc = new GPSTimeFormat(clock);
    }

    /**
     * Set the TOE
     * @param clock GPSTimeFormat for the TOE
     */
    public void setTOE(GPSTimeFormat ephem){
        this.toe = new GPSTimeFormat(ephem);
    }

    /**
     * Set the clock parameters
     * @param bias clock bias
     * @param drift clock drift
     * @param driftrate clock drift rate
     */
    public void setClockParams(double bias, double drift, double driftrate){
        this.clockBias      = bias;
        this.clockDrift     = drift;
        this.clockDriftRate = driftrate;
    }

    /**
     * Set the Harmonic Corrections
     * @param crc CRC
     * @param crs CRS
     * @param cuc CUC
     * @param cus CUS
     * @param cic CIC
     * @param cis CIS
     */
    public void setHarmonicCorrections(double crc, double crs, double cuc, double cus, double cic, double cis){
        this.Crc            = crc;
        this.Crs            = crs;
        this.Cuc            = cuc;
        this.Cus            = cus;
        this.Cic            = cic;
        this.Cis            = cis;

    }

    /**
     * Set the orbit elements
     * @param sqrta square root of semi-major axis
     * @param e eccentricity
     * @param i inclination
     * @param raan right ascension of ascending node
     * @param w argument of perigee
     * @param m0 initial mean anomaly
     */
    public void setOrbit(double sqrta, double e, double i, double raan, double w, double m0){
        this.sqrtA = sqrta;
        this.ecc = e;
        this.inc = i;
        this.omega = raan;
        this.argp = w;
        this.ma = m0;
    }

    /**
     * Set the orbit corrections
     * @param dn delta N
     * @param od omega dot
     * @param id inclination dot
     */
    public void setOrbitCorrections( double dn, double od, double id){
        this.deltaN = dn;
        this.omegadot = od;
        this.idot = id;
    }

    /**
     * Return the TOE in MJD
     * @return TOE in MJD
     */
    public double getTOEmjd(){
        return this.toe.mjd();
    }
    /**
     * Return the PRN
     * @return the PRN
     */
    public int prn(){
        return this.prn;
    }

    /** Method of printing out the GPS SV data.
     * @param title Title or name for the SV.
     */
    public void print(String title){
        System.out.println("GPS SV: "+title);
        System.out.println("   PRN: "+this.prn);
        this.toc.printCalDate("   Time of Clock");
        System.out.println("toc mjd = "+this.toc.mjd());
        this.toe.printCalDate("   Time of Ephemeris");
        this.toe.print("TOE");
        System.out.println("toe mjd = "+this.toe.mjd());
        System.out.println("   clockBias: "+this.clockBias);
        System.out.println("   clockDrift: "+this.clockDrift);
        System.out.println("   clockDriftRate: "+this.clockDriftRate);
        System.out.println("   Crc: "+this.Crc);
        System.out.println("   Crs: "+this.Crs);
        System.out.println("   Cic: "+this.Cic);
        System.out.println("   Cis: "+this.Cis);
        System.out.println("   Cuc: "+this.Cuc);
        System.out.println("   Cus: "+this.Cus);
        System.out.println("   sqrtA: "+this.sqrtA);
        System.out.println("   deltaN: "+this.deltaN);
        System.out.println("   ecc: "+this.ecc);
        System.out.println("   inclination: "+this.inc);
        System.out.println("   idot: "+this.idot);
        System.out.println("   omega: "+this.omega);
        System.out.println("   omegadot: "+this.omegadot);
        System.out.println("   argument of perigee: "+this.argp);
        System.out.println("   mean anomaly: "+this.ma);
    }

    /** Compute the ECEF position vector of a GPS SV at a particular MJD.
     * @param mjd Modified Julian Date.
     * @return ECEF position vector of the GPS SV in meters.
     */
    public VectorN rECEF(double mjd){

        // compute time since Ephemeris Epoch
        double ephemTime = this.toe.mjd_utc();
        double dt = (mjd - ephemTime)*86400.0;//-biasCorrection(new GPSTimeFormat(mjd).gps_sow());
        double t0 = (double) this.toe.gps_week();
//        double dt0 = (mjd - (44244.0 + 7.0*t0))*86400.0;
        double toe_sow = this.toe.gps_sow();

        // compute mean anomaly at current time
        double sma = this.sqrtA * this.sqrtA;
        double n0 = Math.sqrt(Constants.GM_WGS84/(sma*sma*sma));  // mean motion
        double n = n0 + this.deltaN;              // apply mean motion correction
        double mt = this.ma + n*dt;

        // solve Kepler's equation
        double E = TwoBody.solveKepler(mt, this.ecc);

        // compute true anomaly
        double sinE = Math.sin(E);
        double cosE = Math.cos(E);
        double den = 1.0 - this.ecc*cosE;
        double sqrome2 = Math.sqrt(1.0 - this.ecc*this.ecc);
        double sinv = (sqrome2*sinE)/den;
        double cosv = (cosE - this.ecc)/den;

        double f = Math.atan2(sinv, cosv);

        // compute argument of latitude
        double phi = f + this.argp;
        double sin2u = Math.sin(2.0*phi);
        double cos2u = Math.cos(2.0*phi);

        // compute periodic corrections
        double dr = this.Crs*sin2u + this.Crc*cos2u;
        double du = this.Cus*sin2u + this.Cuc*cos2u;
        double di = this.Cis*sin2u + this.Cic*cos2u;

        // apply corrections
        double r = sma*(1.0 - this.ecc*cosE) + dr;
        double u = phi + du;
        double i = this.inc + this.idot*dt + di;
//        double L = this.omega + this.omegadot*dt - Constants.WE_WGS84*dt0;
        // corrected, I hope on 6/26/06 by DEG
        double L = this.omega + (this.omegadot - Constants.WE_WGS84)*dt - Constants.WE_WGS84*toe_sow;

        VectorN rvec = new VectorN(r*Math.cos(u), r*Math.sin(u), 0.0);

        RotationMatrix R = new RotationMatrix(1,-i, 3,-L);

        VectorN out = R.times(rvec);
        return out;
    }

    /** Compute the ECI position vector of a GPS SV at a particular MJD.
     * Note: this really gives the position vector in an inertial frame,
     * which is most likely not the real ECI frame since all we are doing is
     * treating the longitude of ascending node as if it were RAAN.
     * @param mjd Modified Julian Date.
     * @return ECEF position vector of the GPS SV in meters.
     */
    public VectorN rECI(double mjd){

        // compute time since Ephemeris Epoch
        double ephemTime = this.toe.mjd_utc();
        double dt = (mjd - ephemTime)*86400.0;
        double t0 = (double) this.toe.gps_week();
        double dt0 = (mjd - (44244.0 + 7.0*t0))*86400.0;

        // compute mean anomaly at current time
        double sma = this.sqrtA * this.sqrtA;
        double n0 = Math.sqrt(Constants.GM_WGS84/(sma*sma*sma));  // mean motion
        double n = n0 + this.deltaN;              // apply mean motion correction
        double mt = this.ma + n*dt;


        // solve Kepler's equation
        double E = TwoBody.solveKepler(mt, this.ecc);

        // compute true anomaly
        double sinE = Math.sin(E);
        double cosE = Math.cos(E);
        double den = 1.0 - this.ecc*cosE;
        double sqrome2 = Math.sqrt(1.0 - this.ecc*this.ecc);
        double sinv = (sqrome2*sinE)/den;
        double cosv = (cosE - this.ecc)/den;

        double f = Math.atan2(sinv, cosv);

        // compute argument of latitude
        double phi = f + this.argp;
        double sin2u = Math.sin(2.0*phi);
        double cos2u = Math.cos(2.0*phi);

        // compute periodic corrections
        double dr = this.Crs*sin2u + this.Crc*cos2u;
        double du = this.Cus*sin2u + this.Cuc*cos2u;
        double di = this.Cis*sin2u + this.Cic*cos2u;

        // apply corrections
        double r = sma*(1.0 - this.ecc*cosE) + dr;
        double u = phi + du;
        double i = this.inc + this.idot*dt + di;
        //double L = this.omega + this.omegadot*dt;
        double toe_sow = this.toe.gps_sow();
        double L = this.omega + (this.omegadot - Constants.WE_WGS84)*dt - Constants.WE_WGS84*toe_sow;

        VectorN rvec = new VectorN(r*Math.cos(u), r*Math.sin(u), 0.0);

        RotationMatrix R = new RotationMatrix(1,-i, 3,-L);

        VectorN out = R.times(rvec);
        return out;
    }
    /** Compute the WGS84 position vector of a GPS SV at a particular MJD.
     * @param mjd Modified Julian Date.
     * @return ECEF position vector of the GPS SV in meters.
     */
    public VectorN rWGS84(double mjd){

        // compute time since Ephemeris Epoch
        double ephemTime = this.toe.mjd_utc();
        double dt = (mjd - ephemTime)*86400.0;//-biasCorrection(new GPSTimeFormat(mjd).gps_sow());
        double t0 = (double) this.toe.gps_week();
        double dt0 = (mjd - (44244.0 + 7.0*t0))*86400.0;

        // compute mean anomaly at current time
        double sma = this.sqrtA * this.sqrtA;
        double n0 = Math.sqrt(Constants.GM_WGS84/(sma*sma*sma));  // mean motion
        double n = n0 + this.deltaN;              // apply mean motion correction
        double mt = this.ma + n*dt;


        // solve Kepler's equation
        double E = TwoBody.solveKepler(mt, this.ecc);

        // compute true anomaly
        double sinE = Math.sin(E);
        double cosE = Math.cos(E);
        double den = 1.0 - this.ecc*cosE;
        double sqrome2 = Math.sqrt(1.0 - this.ecc*this.ecc);
        double sinv = (sqrome2*sinE)/den;
        double cosv = (cosE - this.ecc)/den;

        double f = Math.atan2(sinv, cosv);

        // compute argument of latitude
        double phi = f + this.argp;
        double sin2u = Math.sin(2.0*phi);
        double cos2u = Math.cos(2.0*phi);

        // compute periodic corrections
        double dr = this.Crs*sin2u + this.Crc*cos2u;
        double du = this.Cus*sin2u + this.Cuc*cos2u;
        double di = this.Cis*sin2u + this.Cic*cos2u;

        // apply corrections
        double r = sma*(1.0 - this.ecc*cosE) + dr;
        double u = phi + du;
        double i = this.inc + this.idot*dt + di;
        double L = this.omega + this.omegadot*dt;

        VectorN rvec = new VectorN(r*Math.cos(u), r*Math.sin(u), 0.0);

        RotationMatrix R = new RotationMatrix(1,-i, 3,-L);

        VectorN out = R.times(rvec);
        return out;
    }
    /** Compute the ECI position vector of a GPS SV at a particular MJD.
     * Note: this really gives the position vector in an inertial frame,
     * which is most likely not the real ECI frame since all we are doing is
     * treating the longitude of ascending node as if it were RAAN.
     * @param mjd Modified Julian Date.
     * @return ECI position vector of the GPS SV in meters.
     */
    public VectorN rvECI(double mjd){

        // compute time since Ephemeris Epoch
        double ephemTime = this.toe.mjd_utc();
        double dt = (mjd - ephemTime)*86400.0;//-biasCorrection(new GPSTimeFormat(mjd).gps_sow());
        double t0 = (double) this.toe.gps_week();
        double dt0 = (mjd - (44244.0 + 7.0*t0))*86400.0;

        // compute mean anomaly at current time
        double sma = this.sqrtA * this.sqrtA;
        double n0 = Math.sqrt(Constants.GM_WGS84/(sma*sma*sma));  // mean motion
        double n = n0 + this.deltaN;              // apply mean motion correction
        double mt = this.ma + n*dt;


        // solve Kepler's equation
        double Etemp = TwoBody.solveKepler(mt, this.ecc);
        double E = MathUtils.Modulo(Etemp, 2.0*MathUtils.PI);


        // compute true anomaly
        double sinE = Math.sin(E);
        double cosE = Math.cos(E);
        double den = 1.0 - this.ecc*cosE;
        double sqrome2 = Math.sqrt(1.0 - this.ecc*this.ecc);
        double sinv = (sqrome2*sinE)/den;
        double cosv = (cosE - this.ecc)/den;

        double f = Math.atan2(sinv, cosv);

        // compute argument of latitude
        double phi = f + this.argp;
        double sin2u = Math.sin(2.0*phi);
        double cos2u = Math.cos(2.0*phi);

        // compute periodic corrections
        double dr = this.Crs*sin2u + this.Crc*cos2u;
        double du = this.Cus*sin2u + this.Cuc*cos2u;
        double di = this.Cis*sin2u + this.Cic*cos2u;

        // apply corrections
        double r = sma*(1.0 - this.ecc*cosE) + dr;
        double u = phi + du;
        double i = this.inc + this.idot*dt + di;
        double L = this.omega + this.omegadot*dt;

        double cosu = Math.cos(u);
        double sinu = Math.sin(u);

        // position in PQW
        VectorN rpqw = new VectorN(r*cosu, r*sinu, 0.0);

        // rotation from PQW to ECI
        RotationMatrix R = new RotationMatrix(1,-i, 3,-L);

        // position in ECI
        VectorN reci = R.times(rpqw);

        // compute the velocity vector
        double denom = 1.0 - ecc*Math.cos(E);
        double Edot = n / denom;
        double num = 1.0 + ecc*Math.cos(f);
        double fdot = (n * num) / (sqrome2 * denom);
        double drdot = 2.0*fdot*(this.Crs*cos2u-this.Cuc*sin2u);
        if (E < 0.0) E = E + 2.0*MathUtils.PI;
        double sqrE = Math.sqrt(E);
        double ae = sma*ecc*Math.sqrt(E)*Edot;
        double rdot =  ae + drdot;
        double dudot = 2.0*fdot*(this.Cus*cos2u-this.Cuc*sin2u);
        double didot = 2.0*fdot*(this.Cis*cos2u-this.Cic*sin2u);
        double udot = fdot + dudot;
//        System.out.println("sma = "+sma+" ecc = "+ecc+" E = "+E+" sqrE = "+sqrE+" Edot = "+Edot);

        // rdot vector in PQW
        double r1 = rdot*cosu - r*udot*sinu;
        double r2 = rdot*sinu + r*udot*cosu;
        VectorN rdotpqw = new VectorN(r1, r2, 0.0);

        // decompose the R matrix into two parts: A and B
        RotationMatrix A = new RotationMatrix(3, -L);
        RotationMatrix B = new RotationMatrix(1, -i);

        // find the derivative of A and B
        double sinL = Math.sin(L);
        double cosL = Math.cos(L);
        double sini = Math.sin(i);
        double cosi = Math.cos(i);
        Matrix A1 = new Matrix(3,3);
        A1.set(0, 0, -sinL);
        A1.set(0, 1, -cosi*cosL);
        A1.set(0, 2, cosL*sini);
        A1.set(1, 0, cosL);
        A1.set(1, 1, -cosi*sinL);
        A1.set(1, 2, sini*sinL);
        Matrix Adot = A1.times(omegadot);

        Matrix B1 = new Matrix(3,3);
        B1.set(0, 1, sini*sinL);
        B1.set(0, 2, cosi*sinL);
        B1.set(1, 1, -cosL*sini);
        B1.set(1, 2, -cosi*cosL);
        B1.set(2, 1, cosi);
        B1.set(2, 2, -sini);
        Matrix Bdot = B1.times(this.idot+didot);

        // combine to form Rdot
        Matrix Rdot1 = Adot.times(B);
        Matrix Rdot2 = A.times(Bdot);
        Matrix Rdot = Rdot1.plus(Rdot2);

        // compute the ECI velocity vector
        VectorN term1 = Rdot.times(rpqw);
        VectorN term2 = R.times(rdotpqw);
        VectorN veci = term1.plus(term2);
//        System.out.println("rdotpqw="+rdotpqw.toString());
        VectorN out = new VectorN(reci, veci);
        return out;
    }
    /** Compute the ECI position vector of a GPS SV at a particular MJD.
     * Note: this correctly treats the longitude of ascending node to provide
     * accurate GPS SV ECI position vectors, which requires the ECI to ECEF transformation
     * matrix from jat.spacetime.EarthRef
     * Added on 6/26/06 by DEG
     * @param mjd Modified Julian Date.
     * @param eci2ecef Matrix containing transformation from ECI to ECEF
     * @return ECEF position vector of the GPS SV in meters.
     */
    public VectorN rECI (double mjd, Matrix eci2ecef){
    	VectorN recef = this.rECEF(mjd);
    	VectorN out = eci2ecef.transpose().times(recef);
    	return out;
    }

    /** Compute the ECI position and velocity vectors of a GPS SV at a particular MJD.
     * Note: this correctly treats the longitude of ascending node to provide
     * accurate GPS SV ECI position and velocity vectors, which requires transformation
     * matrices from jat.spacetime.EarthRef
     * Added on 6/26/06 by DEG.
     * Reference: Tak Ebinuma's dissertation.
     * @param mjd Modified Julian Date.
     * @param poleMatrix Matrix containing Pole transformation
     * @param ghaMatrix Matrix containing GHA transformation
     * @param todMatrix Matrix containing the TOD transformation
     * @return ECEF position vector of the GPS SV in meters.
     */
    public VectorN rvECI (double mjd, Matrix poleMatrix, Matrix ghaMatrix, Matrix todMatrix){
        // compute time since Ephemeris Epoch
        double ephemTime = this.toe.mjd_utc();
        double dt = (mjd - ephemTime)*86400.0;//-biasCorrection(new GPSTimeFormat(mjd).gps_sow());
        //double t0 = (double) this.toe.gps_week();
        //double dt0 = (mjd - (44244.0 + 7.0*t0))*86400.0;
        double toe_sow = this.toe.gps_sow();

        // compute mean anomaly at current time
        double sma = this.sqrtA * this.sqrtA;
        double n0 = Math.sqrt(Constants.GM_WGS84/(sma*sma*sma));  // mean motion
        double n = n0 + this.deltaN;              // apply mean motion correction
        double mt = this.ma + n*dt;

        // solve Kepler's equation
        double E = TwoBody.solveKepler(mt, this.ecc);

        // compute true anomaly
        double sinE = Math.sin(E);
        double cosE = Math.cos(E);
        double den = 1.0 - this.ecc*cosE;
        double sqrome2 = Math.sqrt(1.0 - this.ecc*this.ecc);
        double sinv = (sqrome2*sinE)/den;
        double cosv = (cosE - this.ecc)/den;
        double f = Math.atan2(sinv, cosv);

        // compute argument of latitude
        double phi = f + this.argp;
        double sin2u = Math.sin(2.0*phi);
        double cos2u = Math.cos(2.0*phi);

        // compute periodic corrections
        double dr = this.Crs*sin2u + this.Crc*cos2u;
        double du = this.Cus*sin2u + this.Cuc*cos2u;
        double di = this.Cis*sin2u + this.Cic*cos2u;

        // apply corrections
        double r = sma*(1.0 - this.ecc*cosE) + dr;
        double u = phi + du;
        double i = this.inc + this.idot*dt + di;
//        double L = this.omega + this.omegadot*dt - Constants.WE_WGS84*dt0;
        double L = this.omega + (this.omegadot - Constants.WE_WGS84)*dt - Constants.WE_WGS84*toe_sow;

        double cosu = Math.cos(u);
        double sinu = Math.sin(u);

        // Position in orbit frame
        VectorN rpqw = new VectorN(r*cosu, r*sinu, 0.0);

        // Rotation from orbit frame to ECEF
        RotationMatrix M = new RotationMatrix(1,-i, 3,-L);
        // Compute ECEF position
        VectorN recef = M.times(rpqw);
        
        // Transform ECEF position to ECI
        Matrix temp = poleMatrix.times(ghaMatrix);
        Matrix eci2ecef = temp.times(todMatrix);
    	VectorN reci = eci2ecef.transpose().times(recef);

    	// compute the velocity vector
        double denom = 1.0 - ecc*Math.cos(E); 
        double Edot = n / denom;
        double num = 1.0 + ecc*Math.cos(f);
        double fdot = (n * num) / (sqrome2 * denom);
        double drdot = 2.0*fdot*(this.Crs*cos2u-this.Cuc*sin2u);
        if (E < 0.0) E = E + 2.0*MathUtils.PI;
        double sqrE = Math.sqrt(E);
        double ae = sma*ecc*sqrE*Edot;
        double rdot =  ae + drdot;
        double dudot = 2.0*fdot*(this.Cus*cos2u-this.Cuc*sin2u);
        double didot = 2.0*fdot*(this.Cis*cos2u-this.Cic*sin2u);
        double udot = fdot + dudot;
//        System.out.println("sma = "+sma+" ecc = "+ecc+" E = "+E+" sqrE = "+sqrE+" Edot = "+Edot);

        // rdot vector in PQW
        double r1 = rdot*cosu - r*udot*sinu;
        double r2 = rdot*sinu + r*udot*cosu;
        VectorN rdotpqw = new VectorN(r1, r2, 0.0);

        // find the derivative of M
        double sinL = Math.sin(L);
        double cosL = Math.cos(L);
        double sini = Math.sin(i);
        double cosi = Math.cos(i);
        Matrix A1 = new Matrix(3,3);
        A1.set(0, 0, -sinL);
        A1.set(0, 1, -cosi*cosL);
        A1.set(0, 2, cosL*sini);
        A1.set(1, 0, cosL);
        A1.set(1, 1, -cosi*sinL);
        A1.set(1, 2, sini*sinL);
        Matrix Adot = A1.times(omegadot - Constants.WE_WGS84);
        Matrix B1 = new Matrix(3,3);
        B1.set(0, 1, sini*sinL);
        B1.set(0, 2, cosi*sinL);
        B1.set(1, 1, -cosL*sini);
        B1.set(1, 2, -cosi*cosL);
        B1.set(2, 1, cosi);
        B1.set(2, 2, -sini);
        Matrix Bdot = B1.times(this.idot+didot);
        Matrix Mdot = Adot.plus(Bdot);

        // Get necessary transpose matrices
        Matrix Ct = todMatrix.transpose();
        Matrix Bt = poleMatrix.transpose();
        Matrix St = ghaMatrix.transpose();

        // Compute derivative of GHA Matrix (S) and its transpose
        Matrix omegaE = new Matrix(3,3);
        omegaE.set(0, 1, Constants.WE_WGS84);
        omegaE.set(1, 0, -Constants.WE_WGS84);
        Matrix Sdot = omegaE.times(ghaMatrix);
        Matrix Sdott = Sdot.transpose();

        // Form an intermediate term
        Matrix CtStBt = (Ct.times(St)).times(Bt);

        // First term of Eqn 2.23 from Tak's dissertation
        Matrix temp1 = ((Ct.times(Sdott)).times(Bt)).times(M);
        Matrix temp2 = CtStBt.times(Mdot);
        Matrix temp3 = temp1.plus(temp2);
        VectorN v1 = temp3.times(rpqw);

        // Second term of Eqn 2.23 from Tak's dissertaion
        Matrix temp4 = CtStBt.times(M);
        VectorN v2 = temp4.times(rdotpqw);

        // ECI velocity
        VectorN veci = v1.plus(v2);

        // Form output vector
        VectorN out = new VectorN(reci, veci);
        return out;
    }
    /** Compute the ECEF position and velocity vectors of a GPS SV at a particular MJD.
     * Note: this correctly treats the longitude of ascending node to provide
     * accurate GPS SV ECEF position and velocity vectors, which requires transformation
     * matrices from jat.spacetime.EarthRef
     * Added on 6/26/06 by DEG.
     * Reference: Tak Ebinuma's dissertation.
     * @param mjd Modified Julian Date.
     * @param poleMatrix Matrix containing Pole transformation
     * @param ghaMatrix Matrix containing GHA transformation
     * @param todMatrix Matrix containing the TOD transformation
     * @return ECEF position vector of the GPS SV in meters.
     */
    public VectorN rvECEF (double mjd){
        // compute time since Ephemeris Epoch
        double ephemTime = this.toe.mjd_utc();
        double dt = (mjd - ephemTime)*86400.0;//-biasCorrection(new GPSTimeFormat(mjd).gps_sow());
        double t0 = (double) this.toe.gps_week();
        double dt0 = (mjd - (44244.0 + 7.0*t0))*86400.0;
        double toe_sow = this.toe.gps_sow();

        // compute mean anomaly at current time
        double sma = this.sqrtA * this.sqrtA;
        double n0 = Math.sqrt(Constants.GM_WGS84/(sma*sma*sma));  // mean motion
        double n = n0 + this.deltaN;              // apply mean motion correction
        double mt = this.ma + n*dt;

        // solve Kepler's equation
        double E = TwoBody.solveKepler(mt, this.ecc);

        // compute true anomaly
        double sinE = Math.sin(E);
        double cosE = Math.cos(E);
        double den = 1.0 - this.ecc*cosE;
        double sqrome2 = Math.sqrt(1.0 - this.ecc*this.ecc);
        double sinv = (sqrome2*sinE)/den;
        double cosv = (cosE - this.ecc)/den;
        double f = Math.atan2(sinv, cosv);

        // compute argument of latitude
        double phi = f + this.argp;
        double sin2u = Math.sin(2.0*phi);
        double cos2u = Math.cos(2.0*phi);

        // compute periodic corrections
        double dr = this.Crs*sin2u + this.Crc*cos2u;
        double du = this.Cus*sin2u + this.Cuc*cos2u;
        double di = this.Cis*sin2u + this.Cic*cos2u;

        // apply corrections
        double r = sma*(1.0 - this.ecc*cosE) + dr;
        double u = phi + du;
        double i = this.inc + this.idot*dt + di;
//        double L = this.omega + this.omegadot*dt - Constants.WE_WGS84*dt0;
        double L = this.omega + (this.omegadot - Constants.WE_WGS84)*dt - Constants.WE_WGS84*toe_sow;

        double cosu = Math.cos(u);
        double sinu = Math.sin(u);

        // Position in orbit frame
        VectorN rpqw = new VectorN(r*cosu, r*sinu, 0.0);

        // Rotation from orbit frame to ECEF
        RotationMatrix M = new RotationMatrix(1,-i, 3,-L);

        // Compute ECEF position
        VectorN recef = M.times(rpqw);
        
        // Transform ECEF position to ECI
        //Matrix temp = poleMatrix.times(ghaMatrix);
        //Matrix eci2ecef = temp.times(todMatrix);
    	//VectorN reci = eci2ecef.transpose().times(recef);

    	// compute the velocity vector
        double denom = 1.0 - ecc*Math.cos(E);
        double Edot = n / denom;
        double num = 1.0 + ecc*Math.cos(f);
        double fdot = (n * num) / (sqrome2 * denom);
        double drdot = 2.0*fdot*(this.Crs*cos2u-this.Cuc*sin2u);
        if (E < 0.0) E = E + 2.0*MathUtils.PI;
        double sqrE = Math.sqrt(E);
        double ae = sma*ecc*sqrE*Edot;
        double rdot =  ae + drdot;
        double dudot = 2.0*fdot*(this.Cus*cos2u-this.Cuc*sin2u);
        double didot = 2.0*fdot*(this.Cis*cos2u-this.Cic*sin2u);
        double udot = fdot + dudot;
//        System.out.println("sma = "+sma+" ecc = "+ecc+" E = "+E+" sqrE = "+sqrE+" Edot = "+Edot);

        // rdot vector in PQW
        double r1 = rdot*cosu - r*udot*sinu;
        double r2 = rdot*sinu + r*udot*cosu;
        VectorN rdotpqw = new VectorN(r1, r2, 0.0);

        VectorN vecef = M.times(rdotpqw);
        
        VectorN out = new VectorN(recef, vecef);
        return out;
    }

    /**
     * Accounts for the biasCorrection term in the GPS broadcast navigation file.
     * @param ts_SOW (transmit time seconds of the GPS week
     * @return biasCorrection [sec]
     */
	public double biasCorrection(double ts_SOW) {
		double ttoc = ts_SOW - this.toc.gps_sow();
		double dtsv = this.clockBias + this.clockDrift*(ttoc)+this.clockDriftRate*(ttoc*ttoc);						
		return dtsv;
	}
}
