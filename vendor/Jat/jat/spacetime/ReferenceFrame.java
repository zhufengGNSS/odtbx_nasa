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
package jat.spacetime;

import java.io.Serializable;

/**
 * Interface used to represent any reference frame.  Subclasses must not
 * only know how to represent positions, vectors, and angles in their
 * reference frame, but also how to translate to other reference frames.
 * 
 * @author Rob Antonucci
 *
 */
public interface ReferenceFrame extends Serializable {
    
    /**
     * Creates a translate that can translate between two reference
     * frames at a given time
     * @param other another reference frame
     * @param t the time at which translation will be done
     * @return a translating object or null if the reference frame
     * does not know how to translate to the other reference frame.
     */
    public ReferenceFrameTranslater getTranslater(ReferenceFrame other, Time t);   
}
