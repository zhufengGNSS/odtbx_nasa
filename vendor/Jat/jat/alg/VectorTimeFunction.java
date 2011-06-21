/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2002 United States Government as represented by the
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

package jat.alg;
import jat.matvec.data.*;

/**
 * The VectorTimeFunction interface provides the mechanism for passing a method
 * that evaluates a function to a solver.
 *
 * @author Tobias Berthold
 * @version 1.0
 */
public interface VectorTimeFunction
{
	/**
	 * @param x VectorN containing the required data.
	 * @param t separate data, typically time
	 * @return VectorN containing the result of the function.
	 */
	public VectorN evaluate(VectorN x, double t);
}
