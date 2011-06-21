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
package jat.spacecraft;

import jat.alg.integrators.*;
import jat.forces.*;
import jat.matvec.data.VectorN;

/**
 * A template class intended to be extended to implement various control
 * laws.  This class represents a unity feedback control law.
 * 
 * spacecraft[]-->(sum)-->[Dynamics]->[Control Law]-->acceleration[]
 *                  A -                             |
 *                  |_______________________________|
 * 
 * @author Richard C. Page III
 */
public class ControlLaw {
    
    public ControlLaw(){    }

    public double[] compute_control(double t, double[] x){
        return new double[x.length];
    }
    
    public double[] compute_control(double t, double[] x, double[] xrel){
        return new double[3];
    }
    
}
