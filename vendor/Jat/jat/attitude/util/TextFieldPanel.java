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
 */
package jat.attitude.util;

import java.awt.*;
import java.awt.event.*;


/**
 * <P>
 * The TextFieldPanel provides an efficient way of creating and organizing
 * TextField objects.
 * The main objective of this class is to create input fields with labels 
 * identifying the fields. 
 *
 * @author Noriko Takada
 * @version 1.0
 */



public class TextFieldPanel extends Panel // TextFieldPanel is itself a panel
{
	 // Create Font variables
	 Font titlef = new Font("Dialog", Font.BOLD, 16);
	 Font boldf = new Font("Dialog", Font.BOLD, 12);
	 Font italicf = new Font("Dialog", Font.ITALIC, 12);
	 Font normalf = new Font("Dialog", Font.PLAIN, 12);
	 
	 // theFrame is used in main() for demonstration purpose
	 private static Frame theFrame;
	 
	 GridBagConstraints constraint = new GridBagConstraints();		
			
	 /**
	  *	Construct a CheckboxPanel object
	  * @param	title		(String) Title of the panel
	  * @param	labels		(String[])	Labels
	  * @param	fields		(Strnig[])	Array of TextField
	  * @param	c			(Color)	Background color of the panel
	  * 
	  */
	 public TextFieldPanel(String title, String labels[], TextField fields[], Color c) 
	 {
  		// Since labels are supposed to be associated with the text fields,
  		// labels[] and fields[] are assumed to have the same length.
  		super();
  		int length = labels.length;
  		
  		setLayout(new GridBagLayout());  		
  		constraint.fill = GridBagConstraints.BOTH; // components grow in both dimensions.
		constraint.insets = new Insets(5,5,5,5);   // 5-pixel mergins on all sides
  		this.setBackground(c);
  		Label titleLabel = new Label(title);
		titleLabel.setFont(titlef);
		setConstraint(0, 0, 3, 1, 0.0, 0.0);
		add(titleLabel, constraint);
  		
  		  		
  		for(int i=0;i<length;++i)
  		{
   			Label label = new Label(labels[i], Label.RIGHT);	
      				// column, row	
      		setConstraint(0, i+1, 1, 1, 0.0, 0.0);
			add(label ,constraint);	
			
			setConstraint(1, i+1, 1, 1, 0.0, 0.0);		
   			add(fields[i], constraint);
   		}
  		
  		Label label = new Label("");	
      	setConstraint(0, length+1, 2, 1, 1.0, 1.0);
		add(label ,constraint);	
  		
 	 }// End of constructor
 
	/**
	 * Demonstrates the use of TextFieldPanel class
	 * @param	args	(String[])	Argument
	 */
	public static void main(String[] args) 
   	{
      	String title = "This is the title.";
      	String hobby[] = {"hobby1", "hobby2","hobby3"};
      	TextField hobbyfields[] = new TextField[3];
      		hobbyfields[0] = new TextField(5);
      		hobbyfields[1] = new TextField(5);
      		hobbyfields[2] = new TextField(5);
      	
      	TextFieldPanel thePanel = new TextFieldPanel(title, hobby, hobbyfields, Color.red);
      
      	      	
      	theFrame = new Frame();
      	theFrame.setTitle("TextFieldPanel Demo");
      	theFrame.addWindowListener(new WindowAdapter(){
      	public void windowClosing(WindowEvent e){
      		System.exit(0);}
      	});
      	
      	theFrame.add( thePanel, BorderLayout.CENTER );
     
      	theFrame.setSize( 200, 300 );
      	Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
      	theFrame.setLocation( (d.width - theFrame.getSize().width) / 2,
         (d.height - theFrame.getSize().height) / 2);
      	theFrame.setVisible( true );
   	} 
	
	/**
	 * Sets the necessary constraint for the GridBagLayout layout
	 * @param	gridx		(int)
	 * @param	gridy		(int)
	 * @param	gridwidth 	(int)
	 * @param	gridheight	(int)
	 * @param	weightx		(int)
	 * @param	weighty		(int)
	 */	
	void setConstraint(int gridx, int gridy,int gridwidth,
	                   int gridheight, double weightx, double weighty)
	{ 
		constraint.gridx = gridx;
		constraint.gridy = gridy;
		constraint.gridwidth = gridwidth;
		constraint.gridheight = gridheight;
		constraint.weightx = weightx;
		constraint.weighty = weighty;
	}
    
	
	
    /*
	public boolean getState(String label) 
 	{
  		TextField fields[] = (TextField[])getComponents();
  		for(int i=0;i<fields.length;++i)
   		if(label.equals(fields[i].getLabel())) 
   			return fields[i].getState();
  		return false;
 	}
 
 	public void setState(String label,boolean state) 
 	{
  		TextField fields[] = (TextField[])getComponents();
  		for(int i=0;i<fields.length;++i)
   		if(label.equals(fields[i].getLabel())) 
   			fields[i].setState(state);
 	}
 	*/

}// End of TextFieldPanel

