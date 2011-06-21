/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2009 United States Government as represented by the
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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * LeapSecondTable encapsulates the application of UTC leap second information.
 * Leap seconds are periodically added over time to account for the differences
 * between International Atomic Time (TAI) and Universal Coordinated Time (UTC)
 * as the Earth's rotation slows.  This class reads the leap second data
 * from an external data file and builds up a collection of LeapSecData.<br>
 * Leap seconds can occur at up to four different times of the year but are 
 * usually only added once or twice a year, if they are added at all.  
 * Leap seconds do not have to be added during a year and a UTC correction 
 * can actually add or take away seconds (although this has not yet occurred 
 * in practice). <br>
 * The accumulated number of leap seconds is: TAI = UTC + accumLeapSecs.
 * <br>
 * The leap second data file that this class reads has the following format:
 * <ul>
 * <li>Comment lines begin with "#",</li>
 * <li>All elements are whitespace-delimited,</li>
 * <li>The Julian Date and TAI-UTC accumulated leap seconds are all on one line, 
 * with the Julian Date first, </li>
 * <li>The Julian Date token is "=JD", followed by the Julian Date, e.g. 2441317.5,</li>
 * <li>The TAI-UTC token is "TAI-UTC=", followed by the accumulated leap seconds, e.g. 10.0,</li>
 * <li>The accumulated leap seconds are expected to be integers,</li>
 * <li>Items before and after these elements are ignored,</li>
 * <li>The Julian Dates (lines) must be in time order.</li>
 * </ul>
 * <br>
 * Example line: <br>
 * 1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0       S + (MJD - 41317.) X 0.0      S
 * @author abrown
 *
 */	 
public class LeapSecondTable {

	/**
	 * An immutable class that contains information about a leap second.  Note
	 * that technically leap seconds don't have to be positive and could be
	 * negative, although this has never been required in practice.
	 */
	public class LeapSecData {

		/**
		 * LeapSecData constructor.  Only the parent class should create these.
		 * @param accumLeapSecs The accumulated number of leap seconds.
		 * @param mjd_utc_start The MJD in UTC of when the leap second was inserted.
		 * @param mjd_utc_next The MJD in UTC of the next leap second.
		 */
		LeapSecData(int accumLeapSecs, double mjd_utc_start,
				double mjd_utc_next) {
			super();
			this.mjd_utc_start = mjd_utc_start;
			this.mjd_utc_next = mjd_utc_next;
			this.accumLeapSecs = accumLeapSecs;
		}

		/**
		 * LeapSecData constructor for the last leap second.  mjd_utc_next is set
		 * to Double.MAX_VALUE.
		 * Only the parent class should create these.
		 * @param accumLeapSecs The accumulated number of leap seconds.
		 * @param mjd_utc_start The MJD in UTC of when the leap second was inserted.
		 */
		LeapSecData(int accumLeapSecs, double mjd_utc_start) {
			super();
			this.mjd_utc_start = mjd_utc_start;
			this.mjd_utc_next = java.lang.Double.MAX_VALUE;
			this.accumLeapSecs = accumLeapSecs;
		}

		/**
		 * The MJD in UTC of when the leap second was inserted.  This is the
		 * effective time of the leap second's effects on UTC with respect to 
		 * TAI.
		 */
		public final double mjd_utc_start;

		/**
		 * The MJD in UTC of the next leap second. This is the effective
		 * time of the next leap second and is bookkept here to minimize
		 * searcing.<br>
		 * Will be set equal to java.lang.Double.MAX_VALUE if this
		 * LeapSecData instance describes the latest leap second, otherwise
		 * it is equal to another LeapSecData's mjd_utc_start.
		 */
		public final double mjd_utc_next;

		/**
		 * The accumulated number of leap seconds for the period between
		 * mjd_utc_start (inclusive) and mjd_utc_next (exclusive), where: <br>
		 * TAI = UTC + accumLeapSecs.
		 */
		public final int accumLeapSecs;
		
		/** The human-readable String representation of this object.  Formatted
		 * as:
		 * "start_date MJD UTC, next_date MJD UTC, accumleapsecs s", e.g. 
		 */
		public String toString() {
			return (this.mjd_utc_start + " MJD UTC, " 
					+ this.mjd_utc_next + " MJD UTC, "
					+ this.accumLeapSecs + " s");
		}
	}

	/**
	 * Parses the given file into an array of LeapSecData useful for storing in
	 * the LeapSecondTable. <br>
	 * This routine ignores any line beginning with "#" and calls parseDataString()
	 * to handle the parsing of each line read from the file.  Finally, this method
	 * handles the post-parsing logic to produce a proper leaps List.<br>
	 * Note, this method assumes the data in the given file is in ascending time order.
	 * @param filename The filename to read
	 * @return ArrayList of LeapSecData.
	 * @throws IOException If the file cannot be found, or if there is a 
	 * parsing error while reading the file.
	 */
	protected ArrayList<LeapSecData> parseFile(String filename) throws IOException {

		// the array used to hold the (incomplete) parsed data
		ArrayList<LeapSecData> results = new ArrayList<LeapSecData>();

		FileReader fr;
		// throws FileNotFoundException:
		fr = new FileReader(filename);

		BufferedReader buf = new BufferedReader(fr);

		try {
			/*
			 * Main parse loop
			 */
			while (buf.ready()) {
				String s = buf.readLine();

				if (s.startsWith("#")) {
					continue;
				}

				LeapSecData ld = parseDataString(s);
				if (ld == null) {
					continue;
				}
				else {
					results.add(ld);
				}

			}
		} catch (IOException e) {
			// Attempt to clean up
			try {
				buf.close();
				fr.close();
			} catch (IOException e1) {
				// ignore.
			}
			throw e; // re-throw
		}
		buf.close();
		fr.close();

		/*
		 * Post-parse: reconcile the mjd_utc_next in each LeapSecData item.
		 */
		int sz = results.size();
		
		if (sz == 0) {
			return results; // no usable data parsed
		}
		
		// The array that holds the completed data. 
		ArrayList<LeapSecData> results2 = new ArrayList<LeapSecData>();
		
		LeapSecData d1 = results.get(0);
		LeapSecData d2;
		for (int i = 0; i < sz; i++) {

			if (i == sz-1) {
				// at the last item
				results2.add(i, new LeapSecData(d1.accumLeapSecs, d1.mjd_utc_start));
			}
			else {
				// not at the last item, set with the end date from the next leap sec
				d2 = results.get(i+1); 
				results2.add(i, new LeapSecData(d1.accumLeapSecs, d1.mjd_utc_start, d2.mjd_utc_start));
				d1 = d2; // advance
			}			
		}

		return results2;
	}

	/**
	 * Parses a valid line of text for the "=JD" and "TAI-UTC=" 
	 * whitespace-separated tokens and returns a LeapSecData instance 
	 * with no mjd_utc_next set.  If the line
	 * can't be parsed or doesn't contain data then NULL is returned.
	 * @return LeapSecData without mjd_utc_next, or NULL if the line can't be parsed.
	 */
	protected LeapSecData parseDataString(String s) {
		String[] ss = s.split("\\s+"); // split the string at 1+ whitespace 

		double jd = Double.NaN;
		int tai_utc = Integer.MIN_VALUE;
		int found = 0;

		for (int i = 0; i < ss.length-1; i++) {
			if (found == 2) {
				break;
			}
			if (ss[i].compareToIgnoreCase("=JD") == 0) {
				i++;
				jd = Double.valueOf(ss[i]);
				found++;
				continue;
			}
			if (ss[i].compareToIgnoreCase("TAI-UTC=") == 0) {
				i++;
				tai_utc = Math.round(Float.valueOf(ss[i]));
				found++;
				continue;
			}
		}

		if ((jd == Double.NaN) || (tai_utc == Integer.MIN_VALUE)) {
			return null;
		}
		else {
			// create partial LeapSecData and convert to MJD from JD.
			return new LeapSecData(tai_utc, TimeUtils.JDtoMJD(jd));
		}
	}

	/**
	 * The filename that created the data in this instance.
	 */
	private final String filename;

	/**
	 * The list of LeapSecData, in (start) date order.
	 * 
	 * LeapSecondTable Design Note:
	 * The search logic depends on a properly-build LeapSecData list.
	 * The list must be in mjd_utc_start-time order, the start and next times
	 * must be the same between two LeapSecData items, and the final
	 * LeapSecData item's mjd_utc_next must be set to Double.MAX_VALUE.
	 */
	private final List<LeapSecData> leaps = new ArrayList<LeapSecData>();

	/**
	 * Create a LeapSecondTable instance by reading the given file.
	 * @param filename The full path and filename of the leap second data file.
	 * @throws IOException If there is a problem opening or parsing the given file. 
	 */
	public LeapSecondTable(final String filename) throws IOException {
		super();
		this.filename = filename;
		this.leaps.addAll(parseFile(filename));
	}

	/**
	 * Copy constructor.
	 * Performs a deep copy of the list. 
	 * @param b The leap second table to copy.
	 */
	public LeapSecondTable(final LeapSecondTable b) {
		super();
		this.filename = b.filename;
		this.leaps.addAll(b.leaps);
	}

	/**
	 * Private default constructor.  This class should always be instantiated 
	 * with a data source or a deep copy.
	 */
	@SuppressWarnings("unused")
	private LeapSecondTable() {
		super();
		this.filename = null;
	}

	/**
	 * @return the filename used as the data basis for this instance.
	 */
	public String getFilename() {
		return filename;
	}

	/**
	 * 
	 * @return True if this class has loaded leap second data, false if it has not.
	 */
	public boolean hasLeapSecs() {
		return !this.leaps.isEmpty();
	}

	/**
	 * 
	 * @return The MJD UTC date of the last leap second, or zero if none are loaded.
	 */
	public double latestLeapSec() {
		if (this.hasLeapSecs()) {
			return (this.leaps.get(this.leaps.size()-1)).mjd_utc_start;
		}
		else {
			return 0.0;
		}
	}

	/**
	 * Return the data for the leap second associated with the given UTC
	 * MJD date, if available.  Returns a LeapSecData of zero accumulated 
	 * seconds with zero for the start mdj_utc start time if mjd_utc is 
	 * before the leap second table or if the leap second table is empty.  
	 * If the given UTC MDJ date is after the last leap second then the 
	 * last LeapSecData is returned.
	 * 
	 * @param mjd_utc The MJD UTC date of interest
	 * @return The relevant LeapSecData to the mjd_utc date or a zeroed-out
	 * LeapSecData if before the table or if the table is empty.
	 */
	public LeapSecData getData(double mjd_utc) {
		
		if (this.leaps.isEmpty()) {
			/* 
			 * No table data: 
			 * zero leap secs, 
			 * start at 0.0 MJD UTC,
			 * max value for MJD UTC.
			 */
			return new LeapSecData(0, 0.0);			
		}
		
		LeapSecData retval;
		
		/*
		 * Search loop: Go backwards through the list because the time
		 * of interest is more likely to be later in the list.  With
		 * the retval.mjd_utc_next set to Double.MAX_VALUE, this should
		 * catch times after the table first.
		 */
		for (int i = this.leaps.size()-1; i > -1 ; i--) {
			retval = this.leaps.get(i);
			if ((mjd_utc < retval.mjd_utc_next) && (mjd_utc >= retval.mjd_utc_start)) {
				return retval;
			}
		}
		
		/* 
		 * Final possibility: Before table starts: 
		 * zero leap secs, 
		 * start at 0.0 MJD UTC,
		 * ends at first leap second.
		 */
		retval = this.leaps.get(0);
		return new LeapSecData(0, 0.0, retval.mjd_utc_start);
	}

	/**
	 * Returns the accumulated leap seconds for the given MJD UTC date, if
	 * available, where: <br>
	 * TAI = UTC + accumLeapSecs.<br>
	 * Returns zero if before the first leap second or if the table
	 * is empty.  Holds the value of the last specified leap second for dates
	 * after the last leap second in the table.
	 * @param mjd_utc The MJD UTC date of interest
	 * @return Returns the accumulated leap seconds for the given MJD UTC date, if
	 * available.  Returns zero if before the first leap second or if the table
	 * is empty.
	 */
	public int getAccumLeapSecs(double mjd_utc) {
		return (getData(mjd_utc)).accumLeapSecs;
	}
}
