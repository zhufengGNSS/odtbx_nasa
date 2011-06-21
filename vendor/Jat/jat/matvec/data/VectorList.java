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
 * 
 * File Created on Aug 25, 2003
 */
package jat.matvec.data;

import java.util.ArrayList;
 
/**
 * <P>
 * The VectorList Class provides a type-conscious ArrayList for VectorN objects.
 *
 * @author <a href="mailto:dave.gaylor@emergentspace.com">Dave Gaylor
 * @version 1.0
 */
 
public class VectorList {
	private ArrayList list = new ArrayList();
	
	/**
	 * Add a vector to the list
	 * @param x VectorN to be added
	 */
	public void add(VectorN x) {
		list.add(x);
	}
	
	/**
	 * Get a VectorN from the list
	 * @param index index of the VectorN to be returned.
	 * @return VectorN
	 */
	public VectorN get(int index) {
		return (VectorN) list.get(index);
	}
	
	/**
	 * Return the size of the list.
	 * @return the size of the list.
	 */
	public int size() {
		return list.size();
	}
	
	public boolean hasNext(int index) {
		boolean out = false;
		if (index < (this.size())) {
			out = true;
		}
		return out;
	}

}
