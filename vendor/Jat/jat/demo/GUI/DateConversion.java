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

package jat.demo.GUI;

import org.eclipse.swt.SWT;
import org.eclipse.swt.events.SelectionAdapter;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Control;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Shell;
import org.eclipse.swt.widgets.Text;
import jat.cm.*;

/**
 * SWT GUI example to convert a Calendar date to a Julian Date. The GUI was
 * designed using the SWT-Designer plugin for Eclipse. It is available at
 * http://www.instantiations.com/
 * 
 * To run in Eclipse, choose "Run As SWT Application."
 * 
 * @author Tobias Berthold
 * 
 */
public class DateConversion
{
	private Label label_1_1_1_1;
	private Label label_1_1_1;
	private Label label_1_1;
	private Label label_1;
	private Label calendarDateLabel;
	private Label label;
	private Button convertButton;
	private Text secondText;

	private Text minuteText;

	private Text hourText;

	private Text dayText;

	double cValue;

	private Text monthText;

	private Text JulianDateText;

	private Text yearText;

	long year, month, day, hour, minute, second;

	protected Shell shell;

	/**
	 * Launch the application
	 * 
	 * @param args
	 */
	public static void main(String[] args)
	{
		try
		{
			DateConversion window = new DateConversion();
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
		shell = new Shell(SWT.MIN | SWT.CLOSE | SWT.TITLE);
		shell.setSize(311, 237);
		shell.setText("Julian Date");

		yearText = new Text(shell, SWT.BORDER);
		yearText.setText("1991");
		yearText.setBounds(100, 37, 40, 20);

		convertButton = new Button(shell, SWT.NONE);
		convertButton.addSelectionListener(new SelectionAdapter()
		{
			public void widgetSelected(SelectionEvent e)
			{
				year = Long.parseLong(yearText.getText());
				month = Long.parseLong(monthText.getText());
				day = Long.parseLong(dayText.getText());
				hour = Long.parseLong(hourText.getText());
				minute = Long.parseLong(minuteText.getText());
				second = Long.parseLong(secondText.getText());
				// cValue = Double.parseDouble(yearText.getText());
				double JD = cm.juliandate(year, month, day, hour, minute, second);
				JulianDateText.setText(JD+ "");
			}
		});
		convertButton.setText("Convert");
		convertButton.setBounds(52, 83, 170, 25);

		JulianDateText = new Text(shell, SWT.BORDER);
		JulianDateText.setBounds(14, 148, 240, 25);

		monthText = new Text(shell, SWT.BORDER);
		monthText.setText("12");
		monthText.setBounds(8, 37, 25, 20);

		calendarDateLabel = new Label(shell, SWT.NONE);
		calendarDateLabel.setText("Calendar Date");
		calendarDateLabel.setBounds(31, 15, 70, 15);

		label = new Label(shell, SWT.NONE);
		label.setText("Julian Date");
		label.setBounds(23, 122, 55, 15);

		dayText = new Text(shell, SWT.BORDER);
		dayText.setBounds(54, 37, 25, 20);
		dayText.setText("31");

		hourText = new Text(shell, SWT.BORDER);
		hourText.setBounds(166, 37, 25, 20);
		hourText.setText("15");

		minuteText = new Text(shell, SWT.BORDER);
		minuteText.setBounds(212, 37, 25, 20);
		minuteText.setText("59");

		secondText = new Text(shell, SWT.BORDER);
		secondText.setBounds(258, 37, 25, 20);
		secondText.setText("59");

		label_1 = new Label(shell, SWT.NONE);
		label_1.setText("-");
		label_1.setBounds(41, 39, 5, 15);

		label_1_1 = new Label(shell, SWT.NONE);
		label_1_1.setBounds(87, 39, 5, 15);
		label_1_1.setText("-");

		label_1_1_1 = new Label(shell, SWT.NONE);
		label_1_1_1.setBounds(199, 39, 5, 15);
		label_1_1_1.setText(":");

		label_1_1_1_1 = new Label(shell, SWT.NONE);
		label_1_1_1_1.setBounds(245, 39, 5, 15);
		label_1_1_1_1.setText(":");
		shell.setTabList(new Control[] {monthText, dayText, yearText, hourText, minuteText, secondText, convertButton});
		//
	}

}
