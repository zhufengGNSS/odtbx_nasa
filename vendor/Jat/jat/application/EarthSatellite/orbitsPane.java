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
 */

package jat.application.EarthSatellite;

import org.eclipse.swt.SWT;
import org.eclipse.swt.custom.StackLayout;
import org.eclipse.swt.layout.FillLayout;
import org.eclipse.swt.layout.FormAttachment;
import org.eclipse.swt.layout.FormData;
import org.eclipse.swt.layout.FormLayout;
import org.eclipse.swt.layout.GridData;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.layout.RowData;
import org.eclipse.swt.layout.RowLayout;
import org.eclipse.swt.widgets.*;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.List;

public class orbitsPane extends Composite
{
  private List list;
  public orbitsPane(Composite parent)
  {
    super(parent, SWT.NONE);
  	setLayout(new GridLayout());
  	final Label orbitsLabel = new Label(this, SWT.NONE);
  	orbitsLabel.setLayoutData(new GridData(GridData.BEGINNING, GridData.CENTER, true, false));
  	orbitsLabel.setText("Orbits");

  	list = new List(this, SWT.BORDER);
  	list.setItems(new String[] {"Orbit1", "Orbit2", "Transfer Orbit"});
  	list.setData("newKey", null);
  	list.setLayoutData(new GridData(GridData.FILL, GridData.FILL, true, true));
  }
}
