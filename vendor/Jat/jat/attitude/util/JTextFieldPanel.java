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

import jat.util.FileUtil;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
/**
 * <P>
 * The JTextFieldPanel provides an efficient way of creating and organizing
 * TextField objects. The class extends JPanel and is therefore a JPanel itself.
 * And all the components constructed in this class is swing components.
 * The main objective of this class is to create input fields with labels 
 * identifying the fields.
 * <P>
 * JTextFieldPanel uses BoxLayout to construct a JPanel with 1 row or column of 
 * JTextField components along with labels.
 *
 * @author 	Noriko Takada
 * @version 	Last Modified: 08/12003
 */

public class JTextFieldPanel extends JPanel // This class is the JPanel itself
{
	 // Instantiate a Font object 
     Font fancyFont = new Font("Serif", Font.BOLD | Font.ITALIC, 32);
	 Font titlef = new Font("Dialog", Font.BOLD, 16);
	 Font boldf = new Font("Dialog", Font.BOLD, 12);
	 Font italicf = new Font("Dialog", Font.ITALIC, 12);
	 Font normalf = new Font("Dialog", Font.PLAIN, 12);
	 
	 
	 // theFrame is used in main() for demonstration purpose
	 private static JFrame theFrame;
	 
	 
	 /**
	  * @param		direction	(int)1-> X_AXIS, 2-> Y_AXIS
	  * @param		title		(String) Title of the panel
	  * @param		labels		(String[]) Labels of the corresponding TextFields
	  * @param		fields		(JTextField[]) TextFields to be laid out
	  * @param		c			(Color)	Color of the panel
	  */
	 public JTextFieldPanel(int direction, String title, String labels[], 
	 						 JTextField fields[], Color c) 
	 {
  		// Since labels are supposed to be associated with the text fields,
  		// labels[] and fields[] are assumed to have the same length.
  		super();
  		int length = labels.length;
  		
  		if(direction == 1)
  		{
  			setLayout(new BoxLayout(this, BoxLayout.X_AXIS)); 
  			int panelRowWidth = length*100 + length * 5 + 5+5;
  			int panelRowHeight = 30;
  			Dimension rowSize = new Dimension (panelRowWidth, panelRowHeight);
  			this.setMaximumSize(rowSize);
  			this.setPreferredSize(rowSize);
  			this.setMinimumSize(rowSize);
  		}			
  		else
  		{
  			setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
  			int panelColumnWidth = 110;
  			int panelColumnHeight = (length+1)*25 + length * 5 + 5;
  			Dimension columnSize = new Dimension(panelColumnWidth, panelColumnHeight);
  			this.setMaximumSize(columnSize);
  			this.setPreferredSize(columnSize);
  			this.setMinimumSize(columnSize);
  		}	
  			
  		this.setBackground(c);
  		this.setBorder(BorderFactory.createTitledBorder(title));
  		
  		 		
  			/* Add JLabel and JTextField component to this JPanel */  		
  			for(int i=0;i<length;++i)
  			{
   				int fieldWidth = 50;
   				int fieldHeight = 25;
   				JPanel unitPanel = new JPanel();
   					unitPanel.setLayout(new BoxLayout(unitPanel, BoxLayout.X_AXIS));
   					Dimension unitSize = new Dimension (2*fieldWidth, fieldHeight);
   					unitPanel.setMaximumSize(unitSize);
   					unitPanel.setPreferredSize(unitSize);
   					unitPanel.setMinimumSize(unitSize);
   					unitPanel.setBackground(c);
   					if(direction ==1) //X direction
   						unitPanel.setAlignmentY(CENTER_ALIGNMENT);
   					else
   						unitPanel.setAlignmentX(LEFT_ALIGNMENT);
   						
   				Dimension fieldSize = new Dimension(fieldWidth, fieldHeight);
   				fields[i].setMaximumSize(fieldSize);
   				fields[i].setPreferredSize(fieldSize);
   				fields[i].setMinimumSize(fieldSize);
   						
   				JLabel label = new JLabel(labels[i], Label.RIGHT);	
   				unitPanel.add(Box.createHorizontalGlue());
      			unitPanel.add(label );
      			unitPanel.add(fields[i]);
      			
      			this.add(Box.createRigidArea(new Dimension(5, 5)));      			
      			this.add(unitPanel);
   			}
   			this.add(Box.createRigidArea(new Dimension(5,5)));
  			this.setAlignmentX(LEFT_ALIGNMENT);
  				
 	 }// End of constructor
 
	/*
	/**
	  * Constructor 2
	  * 
	  * @param		direction		1-> X_AXIS
	  * 							2-> Y_AXIS
	  * @param		title			Title of the panel
	  * @param		labels[]		Labels of the corresponding TextFields
	  * @param		fields[]		TextFields to be laid out
	  * @param		c				Color of the panel
	  * @param		titleImage		Icon image to be added with the title label
	  */
	 /*
	 
	 public JTextFieldPanel(int direction, String title, String labels[], 
	 						 JTextField fields[], Color c, ImageIcon titleImage) 
	 {
  		// Since labels are supposed to be associated with the text fields,
  		// labels[] and fields[] are assumed to have the same length.
  		super();
  		int length = labels.length;
  		
  		// The target of BoxLayout is 'this' JPanel. /
  		if(direction == 1)
  		{
  			setLayout(new BoxLayout(this, BoxLayout.Y_AXIS)); 
  			JPanel rowPanel = new JPanel();
  				rowPanel.setLayout(new BoxLayout(rowPanel, BoxLayout.X_AXIS));
  						  				  				
  			int panelRowWidth = length*150 + length * 5 + 5;
  			int panelRowHeight = 70;
  			Dimension rowSize = new Dimension (panelRowWidth, panelRowHeight);
  			Dimension panelSize = new Dimension (panelRowWidth, panelRowHeight + 50);
  			rowPanel.setMaximumSize(rowSize);
  			rowPanel.setPreferredSize(rowSize);
  			rowPanel.setMinimumSize(rowSize);
  			rowPanel.setAlignmentX(LEFT_ALIGNMENT);
  			this.setMaximumSize(panelSize);
  			this.setPreferredSize(panelSize);
  			this.setMinimumSize(panelSize);
  			  			
  			rowPanel.setBackground(c);
  			this.setBackground(c);
  			this.setBorder(BorderFactory.createLineBorder(Color.red));
  		
  			JLabel titleLabel = new JLabel(title);
  			// Associate the font with the label 
			titleLabel.setFont(titlef);
			titleLabel.setIcon(titleImage);
			// Align the text to the right of the Icon 
    		titleLabel.setHorizontalAlignment(JLabel.RIGHT); 
    		titleLabel.setAlignmentX(LEFT_ALIGNMENT);

  			this.add(titleLabel);
  		
  			// Add JLabel and JTextField component to this JPanel   		
  			for(int i=0;i<length;++i)
  			{
   				JLabel label = new JLabel(labels[i], Label.RIGHT);	
   				rowPanel.add(Box.createHorizontalGlue());
      			rowPanel.add(label );
      			rowPanel.add(fields[i]);
      			
   			}
   			
   			this.add(rowPanel);
  		}// End of if		
  		else
  		{
  			setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
  			int panelColumnWidth = 25;
  			int panelColumnHeight = length*70 + length * 5 + 5;
  			Dimension columnSize = new Dimension(panelColumnWidth, panelColumnHeight);
  			this.setMaximumSize(columnSize);
  			this.setPreferredSize(columnSize);
  			this.setMinimumSize(columnSize);
  			
  		
  			this.setBackground(c);
  			this.setBorder(BorderFactory.createLineBorder(Color.red));
  		
  			JLabel titleLabel = new JLabel(title);
  			// Associate the font with the label 
			titleLabel.setFont(titlef);
			titleLabel.setIcon(titleImage);
			// Align the text to the right of the Icon 
    		titleLabel.setHorizontalAlignment(JLabel.RIGHT); 

  			add(titleLabel);
  		
  			// Add JLabel and JTextField component to this JPanel   		
  			for(int i=0;i<length;++i)
  			{
   				JLabel label = new JLabel(labels[i], Label.RIGHT);
   				add(Box.createHorizontalGlue());	
      			add(label );
      			add(fields[i]);
   			}
   		
  		}//End of else
  		
  	 		
 	 }// End of constructor 2
 	 */
 	 
	/**
	 * Demonstrates the use of JTextFieldPanel class
	 * @param	args	(String[])	Argument
	 */
	public static void main(String[] args) 
   	{	
   		/* Create ImageIcon to be used with the JLabel component */
   		String images_path = FileUtil.getClassFilePath("jat.attitude.thesis", "AttitudeSimulator")+ "images/";
		ImageIcon domburi = new ImageIcon(images_path+"Domburi7.gif");
		
      	String title = "";
      	String hobby[] = {"w1", "hobby2","hobby3"};
      	JTextField hobbyfields[] = new JTextField[3];
      		hobbyfields[0] = new JTextField(5);
      		hobbyfields[1] = new JTextField(5);
      		hobbyfields[2] = new JTextField(5);
      	
      	JTextFieldPanel theJPanel = new JTextFieldPanel(2,title, hobby, hobbyfields, Color.pink);
      
      	      	
      	theFrame = new JFrame();
      	
      	theFrame.addWindowListener(new WindowAdapter(){
      	public void windowClosing(WindowEvent e){
      		System.exit(0);}
      	});
      	JPanel outerPanel = new JPanel();
      	outerPanel.setLayout(new BoxLayout(outerPanel, BoxLayout.Y_AXIS));
      	outerPanel.setSize(750, 550);
      	outerPanel.add(theJPanel);
      	//theFrame.getContentPane().setLayout(new BoxLayout(theFrame,BoxLayout.Y_AXIS));
      	theFrame.getContentPane().add( outerPanel, BorderLayout.CENTER );
      	theFrame.setTitle("JTextFieldPanel");
     	Dimension frameSize = theJPanel.getPreferredSize();
      	//theFrame.setSize((int)frameSize.getWidth(), (int)frameSize.getHeight());
      	theFrame.setSize(750, 550);
      	Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
      	theFrame.setLocation( (d.width - theFrame.getSize().width) / 2,
         (d.height - theFrame.getSize().height) / 2);
      	theFrame.setVisible(true );
   	}// End of main() 
	
	
    
}// End of file
