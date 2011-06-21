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

package jat.test.SolSyspropagator;

import org.eclipse.swt.SWT;
import org.eclipse.swt.custom.StackLayout;
import org.eclipse.swt.events.SelectionAdapter;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.layout.FillLayout;
import org.eclipse.swt.layout.FormAttachment;
import org.eclipse.swt.layout.FormData;
import org.eclipse.swt.layout.FormLayout;
import org.eclipse.swt.layout.GridData;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.layout.RowLayout;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Dialog;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.Group;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Shell;
import org.eclipse.swt.widgets.TabFolder;
import org.eclipse.swt.widgets.TabItem;
import org.eclipse.ui.forms.widgets.ColumnLayout;

public class Model extends Dialog
{

	protected Object result;

	protected Shell shell;

	/**
	 * Create the dialog
	 * @param parent
	 * @param style
	 */
	public Model(Shell parent, int style)
	{
		super(parent, style);
	}

	/**
	 * Create the dialog
	 * @param parent
	 */
	public Model(Shell parent)
	{
		this(parent, SWT.NONE);
	}

	/**
	 * Open the dialog
	 * @return the result
	 */
	public Object open()
	{
		createContents();
		shell.open();
		shell.layout();
		Display display = getParent().getDisplay();
		while (!shell.isDisposed())
		{
			if (!display.readAndDispatch())
				display.sleep();
		}
		return result;
	}

	/**
	 * Create contents of the dialog
	 */
	protected void createContents()
	{
		shell = new Shell(getParent(), SWT.DIALOG_TRIM | SWT.APPLICATION_MODAL);
		shell.setLayout(new ColumnLayout());
		shell.setSize(500, 375);
		shell.setText("Model");

		final TabFolder tabFolder = new TabFolder(shell, SWT.NONE);

		final TabItem originTabItem = new TabItem(tabFolder, SWT.NONE);
		originTabItem.setText("Origin");

		final Group originOfTheGroup = new Group(tabFolder, SWT.NONE);
		originOfTheGroup.setText("Origin of the Reference Frame");
		originTabItem.setControl(originOfTheGroup);

		final Button sunButton = new Button(originOfTheGroup, SWT.RADIO);
		sunButton.setBounds(55, 20, 75, 20);
		sunButton.setText("Sun Center");

		final Button solarSystemBarycenterButton = new Button(originOfTheGroup, SWT.RADIO);
		solarSystemBarycenterButton.setBounds(55, 60, 140, 20);
		solarSystemBarycenterButton.setText("Solar System Barycenter");

		final Button gravitatingBodiesCenterButton = new Button(originOfTheGroup, SWT.RADIO);
		gravitatingBodiesCenterButton.setEnabled(false);
		gravitatingBodiesCenterButton.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e)
			{
			}
		});
		gravitatingBodiesCenterButton.setBounds(55, 100, 185, 20);
		gravitatingBodiesCenterButton.setText("Gravitating Bodies Center of Mass");

		final Button earthCenterButton = new Button(originOfTheGroup, SWT.RADIO);
		earthCenterButton.setBounds(55, 140, 85, 20);
		earthCenterButton.setText("Earth Center");

		final Button earthmoonBarycenterButton = new Button(originOfTheGroup, SWT.RADIO);
		earthmoonBarycenterButton.setBounds(55, 180, 135, 20);
		earthmoonBarycenterButton.setText("Earth-Moon Barycenter");

		final Button moonCenterButton = new Button(originOfTheGroup, SWT.RADIO);
		moonCenterButton.setBounds(55, 220, 85, 20);
		moonCenterButton.setText("Moon Center");

		final Button marsCenterButton = new Button(originOfTheGroup, SWT.RADIO);
		marsCenterButton.setBounds(55, 260, 80, 20);
		marsCenterButton.setText("Mars Center");

		final TabItem bodiesTabItem = new TabItem(tabFolder, SWT.NONE);
		bodiesTabItem.setText("Bodies");

		final Group gravitatingBodiesGroup = new Group(tabFolder, SWT.NONE);
		gravitatingBodiesGroup.setText("Gravitating Bodies");
		bodiesTabItem.setControl(gravitatingBodiesGroup);

		final Button sunButton_1 = new Button(gravitatingBodiesGroup, SWT.CHECK);
		sunButton_1.setText("Sun");
		sunButton_1.setBounds(45, 40, 40, 20);

		final Button earthButton = new Button(gravitatingBodiesGroup, SWT.CHECK);
		earthButton.setText("Earth");
		earthButton.setBounds(45, 100, 50, 20);

		final Button moonButton = new Button(gravitatingBodiesGroup, SWT.CHECK);
		moonButton.setText("Moon");
		moonButton.setBounds(45, 160, 50, 20);

		final Button jupiterButton = new Button(gravitatingBodiesGroup, SWT.CHECK);
		jupiterButton.setText("Jupiter");
		jupiterButton.setBounds(45, 220, 55, 20);

		final TabItem positionTabItem = new TabItem(tabFolder, SWT.NONE);
		positionTabItem.setText("Position");

		final Group group_1 = new Group(tabFolder, SWT.NONE);
		positionTabItem.setControl(group_1);

		final Button circularOrbitButton = new Button(group_1, SWT.RADIO);
		circularOrbitButton.setText("Circular Orbit");
		circularOrbitButton.setBounds(85, 55, 85, 20);

		final Button keplerElementsButton = new Button(group_1, SWT.RADIO);
		keplerElementsButton.setText("Kepler Elements");
		keplerElementsButton.setBounds(85, 130, 100, 20);

		final Button de405EphemeridesButton = new Button(group_1, SWT.RADIO);
		de405EphemeridesButton.setText("DE405 Ephemerides");
		de405EphemeridesButton.setBounds(85, 205, 120, 20);

		final TabItem atmosphereTabItem = new TabItem(tabFolder, SWT.NONE);
		atmosphereTabItem.setText("Atmosphere");

		final Group group_2 = new Group(tabFolder, SWT.NONE);
		atmosphereTabItem.setControl(group_2);

		final Button offButton = new Button(group_2, SWT.RADIO);
		offButton.setText("Off");
		offButton.setBounds(70, 73, 120, 30);

		final Button harrispriesterButton = new Button(group_2, SWT.RADIO);
		harrispriesterButton.setText("Harris-Priester");
		harrispriesterButton.setBounds(70, 176, 120, 30);

		final TabItem timeTabItem = new TabItem(tabFolder, SWT.NONE);
		timeTabItem.setText("Time");

		final TabItem forcesTabItem = new TabItem(tabFolder, SWT.NONE);
		forcesTabItem.setText("Forces");

		final Group group = new Group(tabFolder, SWT.NONE);
		forcesTabItem.setControl(group);

		final Button twobodyButton = new Button(group, SWT.RADIO);
		twobodyButton.setText("Two-Body");
		twobodyButton.setBounds(45, 47, 70, 20);

		final Button nbodyButton = new Button(group, SWT.RADIO);
		nbodyButton.setText("N-Body");
		nbodyButton.setBounds(45, 124, 60, 20);

		final Button circularlRestrictedThreeButton = new Button(group, SWT.RADIO);
		circularlRestrictedThreeButton.setText("Circularl Restricted Three Body Problem");
		circularlRestrictedThreeButton.setBounds(45, 201, 210, 20);

		final Composite composite = new Composite(shell, SWT.NONE);
		composite.setLayout(new FillLayout());

		final Button button = new Button(composite, SWT.NONE);
		button.setText("button");

		final Button button_1 = new Button(composite, SWT.NONE);
		button_1.setText("button");
		//
	}

}
