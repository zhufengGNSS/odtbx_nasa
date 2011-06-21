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

import org.eclipse.jface.dialogs.Dialog;
import org.eclipse.jface.dialogs.IDialogConstants;
import org.eclipse.swt.SWT;
import org.eclipse.swt.graphics.Point;
import org.eclipse.swt.layout.FillLayout;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Control;
import org.eclipse.swt.widgets.Group;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Shell;
import org.eclipse.swt.widgets.TabFolder;
import org.eclipse.swt.widgets.TabItem;
import org.eclipse.swt.widgets.Text;

public class Integrator extends Dialog
{

	private Text erroToleranceText;
	/**
	 * Create the dialog
	 * @param parentShell
	 */
	public Integrator(Shell parentShell)
	{
		super(parentShell);
	}

	/**
	 * Create contents of the dialog
	 * @param parent
	 */
	protected Control createDialogArea(Composite parent)
	{
		Composite container = (Composite) super.createDialogArea(parent);
		container.setLayout(new FillLayout());

		final TabFolder tabFolder = new TabFolder(container, SWT.NONE);

		final TabItem typeTabItem = new TabItem(tabFolder, SWT.NONE);
		typeTabItem.setText("Type");

		final Group group_1 = new Group(tabFolder, SWT.NONE);
		typeTabItem.setControl(group_1);

		final Button rk45Button = new Button(group_1, SWT.RADIO);
		rk45Button.setBounds(70, 74, 50, 20);
		rk45Button.setText("RK45");

		final Button rk78fehlbergButton = new Button(group_1, SWT.RADIO);
		rk78fehlbergButton.setBounds(70, 168, 95, 20);
		rk78fehlbergButton.setText("RK78 Fehlberg");

		final TabItem optionsTabItem = new TabItem(tabFolder, SWT.NONE);
		optionsTabItem.setText("Options");

		final Group group = new Group(tabFolder, SWT.NONE);
		optionsTabItem.setControl(group);

		final Label errorToleranceLabel = new Label(group, SWT.NONE);
		errorToleranceLabel.setText("Error Tolerance");
		errorToleranceLabel.setBounds(16, 50, 120, 30);

		erroToleranceText = new Text(group, SWT.BORDER);
		erroToleranceText.setText("10 E-8");
		erroToleranceText.setBounds(135, 50, 120, 30);
		//
		return container;
	}

	/**
	 * Create contents of the button bar
	 * @param parent
	 */
	protected void createButtonsForButtonBar(Composite parent)
	{
		createButton(parent, IDialogConstants.OK_ID, IDialogConstants.OK_LABEL, true);
		createButton(parent, IDialogConstants.CANCEL_ID, IDialogConstants.CANCEL_LABEL, false);
	}

	/**
	 * Return the initial size of the dialog
	 */
	protected Point getInitialSize()
	{
		return new Point(500, 375);
	}
	protected void configureShell(Shell newShell)
	{
		super.configureShell(newShell);
		newShell.setText("Integrator");
	}

}
