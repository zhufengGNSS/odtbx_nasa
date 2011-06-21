/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2010 United States Government as represented by the
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

package jat.matlabInterface.unittest;

import jat.matlabInterface.BatchCalcs;
import jat.matvec.data.Matrix;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.text.DecimalFormat;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

import junit.framework.TestCase;

import org.junit.Test;

public class BatchCalcsTest extends TestCase {

	/** If true, data files will be generated which can later be used as truth files.
	 * Data files are generated in same directory as truth files and with the same name
	 * but with "_new" inserted before the extension. */
	private static boolean GENERATE_TEST_DATA = false;
	
	// Margins of error
	
	/** Error induced by floating point differences when comparing results from a Matlab and 
	 * a Java algorithm. */
	static final double DEFAULT_ERROR = 0.0001;
	
	/** Error induced by assuming precession and nutation are zero near Jan 1 2000 */
	static final double PREC_NUT_2000_ERROR = 0.0001;
	
	/** Error induced by using old precession and nutation calculations. */
	static final double PREC_NUT_OLD_ERROR = 0.00001;

	@Test
	public void testEciToEcefXform() {
		// Test against sidereal rotation near Jan 1 2000.
		String testDataFile = "jat/matlabInterface/unittest/batchcalcs_eci2ecef1.txt";
		Map<Double, Matrix> matrices = readTestFile(testDataFile);
		double[] timesArray = getTimes(matrices);
		double[][][] results = BatchCalcs.eciToEcefXform(timesArray, 0);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), PREC_NUT_2000_ERROR);
		}

		// Test over a two year period, only occasionally recomputing precession 
		// and nutation 
		testDataFile = "jat/matlabInterface/unittest/batchcalcs_eci2ecef2.txt";
		matrices = readTestFile(testDataFile);
		timesArray = getTimes(matrices);
		results = BatchCalcs.eciToEcefXform(timesArray, 10);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), PREC_NUT_OLD_ERROR);
		}
}

	@Test
	public void testEcefToEciXform() {
		// Test against sidereal rotation near Jan 1 2000.
		String testDataFile = "jat/matlabInterface/unittest/batchcalcs_ecef2eci1.txt";
		Map<Double, Matrix> matrices = readTestFile(testDataFile);
		double[] timesArray = getTimes(matrices);
		double[][][] results = BatchCalcs.ecefToEciXform(timesArray, 0);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), PREC_NUT_2000_ERROR);
		}

		// Test over a two year period, only occasionally recomputing precession 
		// and nutation 
		testDataFile = "jat/matlabInterface/unittest/batchcalcs_ecef2eci2.txt";
		matrices = readTestFile(testDataFile);
		timesArray = getTimes(matrices);
		results = BatchCalcs.ecefToEciXform(timesArray, 10);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), PREC_NUT_OLD_ERROR);
		}
	}

	@Test
	public void testJ2000ToTODXform() {
		// Test over a half day period
		String testDataFile = "jat/matlabInterface/unittest/batchcalcs_tod1.txt";
		Map<Double, Matrix> matrices = readTestFile(testDataFile);
		double[] timesArray = getTimes(matrices);
		double[][][] results = BatchCalcs.j2000ToTODXform(timesArray);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), DEFAULT_ERROR);
		}

		// Test over a two year period
		testDataFile = "jat/matlabInterface/unittest/batchcalcs_tod2.txt";
		matrices = readTestFile(testDataFile);
		timesArray = getTimes(matrices);
		results = BatchCalcs.j2000ToTODXform(timesArray);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), DEFAULT_ERROR);
		}
	}

	@Test
	public void testGetPrecession() {
		// Test over a half day period
		String testDataFile = "jat/matlabInterface/unittest/batchcalcs_precession1.txt";
		Map<Double, Matrix> matrices = readTestFile(testDataFile);
		double[] timesArray = getTimes(matrices);
		double[][][] results = BatchCalcs.getPrecession(timesArray);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), DEFAULT_ERROR);
		}

		// Test over a two year period
		testDataFile = "jat/matlabInterface/unittest/batchcalcs_precession2.txt";
		matrices = readTestFile(testDataFile);
		timesArray = getTimes(matrices);
		results = BatchCalcs.getPrecession(timesArray);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), DEFAULT_ERROR);
		}
	}

	@Test
	public void testGetNutation() {
		// Test over a half day period
		String testDataFile = "jat/matlabInterface/unittest/batchcalcs_nutation1.txt";
		Map<Double, Matrix> matrices = readTestFile(testDataFile);
		double[] timesArray = getTimes(matrices);
		double[][][] results = BatchCalcs.getNutation(timesArray);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), DEFAULT_ERROR);
		}

		// Test over a two year period
		testDataFile = "jat/matlabInterface/unittest/batchcalcs_nutation2.txt";
		matrices = readTestFile(testDataFile);
		timesArray = getTimes(matrices);
		results = BatchCalcs.getNutation(timesArray);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), DEFAULT_ERROR);
		}
	}

	@Test
	public void testGetGHAMatrix() {
		// Test over a half day period
		String testDataFile = "jat/matlabInterface/unittest/batchcalcs_gha1.txt";
		Map<Double, Matrix> matrices = readTestFile(testDataFile);
		double[] timesArray = getTimes(matrices);
		double[][][] results = BatchCalcs.getGHAMatrix(timesArray);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), DEFAULT_ERROR);
		}

		// Test over a two year period
		testDataFile = "jat/matlabInterface/unittest/batchcalcs_gha2.txt";
		matrices = readTestFile(testDataFile);
		timesArray = getTimes(matrices);
		results = BatchCalcs.getGHAMatrix(timesArray);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), DEFAULT_ERROR);
		}
	}

	@Test
	public void testGetPoleMatrix() {
		// Test over a half day period
		String testDataFile = "jat/matlabInterface/unittest/batchcalcs_pole1.txt";
		Map<Double, Matrix> matrices = readTestFile(testDataFile);
		double[] timesArray = getTimes(matrices);
		double[][][] results = BatchCalcs.getPoleMatrix(timesArray);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), DEFAULT_ERROR);
		}

		// Test over a two year period
		testDataFile = "jat/matlabInterface/unittest/batchcalcs_pole2.txt";
		matrices = readTestFile(testDataFile);
		timesArray = getTimes(matrices);
		results = BatchCalcs.getPoleMatrix(timesArray);
		if (GENERATE_TEST_DATA) {
			generateTestFile(testDataFile, timesArray, results);
		}
		for(int ctr=0; ctr<timesArray.length; ++ctr) {
			compareResults(results, ctr, matrices.get(timesArray[ctr]), DEFAULT_ERROR);
		}
}

	private Map<Double, Matrix> readTestFile(String testDataFile) {
		Map<Double, Matrix> matrices = new LinkedHashMap<Double, Matrix>();
		InputStream strm = getClass().getClassLoader().getResourceAsStream(testDataFile);
		BufferedReader rdr = new BufferedReader(new InputStreamReader(strm));
		int lineCtr = 0;
		boolean done = false;
		while (!done) {
			++lineCtr;
			try {
				String nextLine = rdr.readLine();      
				done = (nextLine == null);
				if ((nextLine != null) && !nextLine.startsWith("#") &&
						!nextLine.trim().equals(""))
				{
					String[] doubles = nextLine.trim().split("\\s+");
					assertEquals("Line " + String.valueOf(lineCtr) + " needs 10 number separated " +
							"by spaces.", 10, doubles.length);

					// Time comes in first column
					Double time = new Double(Double.parseDouble(doubles[0]));

					// Next 9 entries of entries in a 3x3 matrix (row ordered)
					Matrix m = new Matrix(3, 3);
					for(int rowctr=0; rowctr<3; ++rowctr) {
						for(int colctr=0; colctr<3; ++colctr) {
							m.set(rowctr, colctr, Double.parseDouble(doubles[3*rowctr+colctr+1]));
						}
					}
					matrices.put(time, m);
				}
			}
			catch (Exception e) {
				fail("Failure to parse doubles on line " + lineCtr + ". " +
						e.getMessage());
			}					
		}
		if (matrices.isEmpty()) {
			fail("Could not parse any entries in test file \"" + testDataFile + "\".");
		}

		return matrices;
	}

	private void generateTestFile(String originalDataFile, double[] times, double[][][] matrices) {
		// First deduce the name and location of the data file.  Should be in the same directory as the
		// truth file and use the same name with "_new" appended.
		String NEW_FILE_APPEND = "_new";
		URL truthUrl = getClass().getClassLoader().getResource(originalDataFile);
		if (truthUrl == null) {
			throw new IllegalArgumentException("Cannot create new data file as original data file, \"" +
					originalDataFile + "\" cannot be found.");
		}
		if (!truthUrl.getProtocol().equals("file")) {
			throw new IllegalArgumentException("Cannot create new data file as original data file, \"" +
					truthUrl.toString() + "\" was not found in a directory that can be written to.");
		}

		File truthFile = new File(truthUrl.getPath());
		String newName = truthFile.getName();
		int extension = newName.lastIndexOf('.');
		if (extension >= 0) {
			newName = newName.substring(0, extension) + NEW_FILE_APPEND + newName.substring(extension);
		}
		else {
			newName = newName + NEW_FILE_APPEND;
		}
		File newFile = new File(truthFile.getParentFile(), newName);


		// Now open a write and output the file.
		DecimalFormat format = new DecimalFormat("0.00000000");
		try {
			BufferedWriter wrtr = new BufferedWriter(new FileWriter(newFile));
			for(int lineCtr=0; lineCtr<times.length; ++lineCtr) {
				wrtr.write(format.format(times[lineCtr]));
				for(int rowCtr=0; rowCtr<3; ++rowCtr) {
					for(int colCtr=0; colCtr<3; ++colCtr) {
						wrtr.write(" " + format.format(matrices[rowCtr][colCtr][lineCtr]));
					}
				}
				wrtr.newLine();
			}
			wrtr.close();
		}
		catch (IOException e) {
			throw new RuntimeException("Error creating data file.", e);
		}					
	}

	/**
	 * This pulls the Doubles out of a Map<Double, Matrix> and returns it as
	 * a double[]
	 * @param matrices Map of Doubles to Matrices
	 * @return array of keys as a double[]
	 */
	private double[] getTimes(Map<Double, Matrix> matrices) {
		double[] timesArray = new double[matrices.size()];
		Iterator<Double> iter = matrices.keySet().iterator();
		for(int ctr=0; iter.hasNext(); ++ctr) {
			timesArray[ctr] = iter.next();
		}
		return timesArray;
	}


	/**
	 * Compare the matrix in the double[][][] with the target matrix
	 * @param results a 3x3xN array which is really an array of 3x3 matrices
	 * @param ctr the index of the desired matrix
	 * @param target the matrix to compare with
	 */
	private void compareResults(double[][][] results, int ctr, Matrix target, double maxError) {
		DecimalFormat fmt = new DecimalFormat("##0.0000000%");
		for(int rowctr=0; rowctr<3; ++rowctr) {
			for(int colctr=0; colctr<3; ++colctr) {
				double resultVal = results[rowctr][colctr][ctr];
				double targetVal = target.get(rowctr, colctr);
				double err = Math.abs(resultVal-targetVal);
				assertTrue("Line " + (ctr+1) + ", coordinate [" + (rowctr+1) + "," + (colctr+1) + 
						"] is off from " + target.get(rowctr, colctr) + " by " + fmt.format(err),  err <= maxError);
			}
		}
	}

}

