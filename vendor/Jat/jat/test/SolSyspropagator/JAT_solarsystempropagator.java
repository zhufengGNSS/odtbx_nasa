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
import org.eclipse.swt.events.SelectionAdapter;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.layout.FillLayout;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.Group;
import org.eclipse.swt.widgets.Menu;
import org.eclipse.swt.widgets.MenuItem;
import org.eclipse.swt.widgets.Shell;
import org.eclipse.swt.widgets.TabFolder;
import org.eclipse.swt.widgets.TabItem;

public class JAT_solarsystempropagator
{

	protected Shell shell;

	/**
	 * Launch the application
	 * @param args
	 */
	public static void main(String[] args)
	{
		try
		{
			JAT_solarsystempropagator window = new JAT_solarsystempropagator();
			window.open();
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	/**
	 * Open the window
	 */
	public void open()
	{
		final Display display = Display.getDefault();
		createContents();
		shell.open();
		shell.layout();
		while (!shell.isDisposed())
		{
			if (!display.readAndDispatch())
				display.sleep();
		}
	}

	/**
	 * Create contents of the window
	 */
	protected void createContents()
	{
		shell = new Shell();
		shell.setLayout(new FillLayout());
		shell.setSize(480, 471);
		shell.setText("JAT Solar System Propagator");

		final Menu menu = new Menu(shell, SWT.BAR);
		shell.setMenuBar(menu);

		final MenuItem fileMenuItem = new MenuItem(menu, SWT.CASCADE);
		fileMenuItem.setText("File");

		final Menu menu_3 = new Menu(fileMenuItem);
		fileMenuItem.setMenu(menu_3);

		final MenuItem newMenuItem = new MenuItem(menu_3, SWT.NONE);
		newMenuItem.setText("New");

		final MenuItem menuItem_1 = new MenuItem(menu_3, SWT.NONE);
		menuItem_1.setText("Menu item");

		final MenuItem menuItem_2 = new MenuItem(menu_3, SWT.NONE);
		menuItem_2.setText("Menu item");

		final MenuItem exitMenuItem = new MenuItem(menu_3, SWT.NONE);
		exitMenuItem.setText("Exit");

		final MenuItem optionsMenuItem = new MenuItem(menu, SWT.CASCADE);
		optionsMenuItem.setText("Options");

		final Menu menu_2 = new Menu(optionsMenuItem);
		optionsMenuItem.setMenu(menu_2);

		final MenuItem modelMenuItem = new MenuItem(menu_2, SWT.NONE);
		modelMenuItem.setText("Model");

		final MenuItem helpMenuItem = new MenuItem(menu, SWT.CASCADE);
		helpMenuItem.setText("Help");

		final Menu menu_1 = new Menu(helpMenuItem);
		helpMenuItem.setMenu(menu_1);

		final MenuItem welcomeMenuItem = new MenuItem(menu_1, SWT.NONE);
		welcomeMenuItem.setText("Welcome");

		final MenuItem aboutMenuItem = new MenuItem(menu_1, SWT.NONE);
		aboutMenuItem.setText("About");
		//
	}

}
