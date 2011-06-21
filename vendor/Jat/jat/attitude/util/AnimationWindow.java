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
 //package thesis;
package jat.attitude.util;

import java.awt.*;
import java.awt.event.*;
import java.text.*;			//Required for NumberFormat
import com.sun.j3d.utils.universe.*;
import com.sun.j3d.utils.behaviors.mouse.*;

import javax.media.j3d.*;
import javax.vecmath.*;
import javax.swing.*;			// Used for Timer & JSlider
import javax.swing.event.*;		// Required for ChangeListener

/**
 * <P>
 * AnimationWindow creates a pop-up Java3D animation window
 * @author 	Daniel Quock	
 * @version	1.3 (03/09/2004)
 * Modification since the last version
 * 		Removed: import java.applet.*;
 * 				 import com.sun.j3d.utils.geometry.*;
 */
//************************************************************************
//*
//*   AnimationWindow class
//*
//*   Written by:					Date:
//*   Daniel Quock				December 2001
//*	Original AnimationWindow
//*   Daniel Quock				July 2002
//*	Moved the AnimationWindow class into it's own file
//*
//*   Description:
//*	This routine creates the pop-up Java3D animation window
//*
//*   Global Variables:
//*	theAnimWindow (HelpWindow) = instance of itself
//*	animFrame	  (JFrame) = frame displayed in theHelpWindow
//*       i_xaxes            (float) = location of reference x-axis
//*       i_yaxes            (float) = location of reference y-axis
//*       i_zaxes            (float) = location of reference z-axis
//*       axeslength         (float) = length of graphical axis lines
//*       xcolor             (Color) = color of x-axis
//*       ycolor             (Color) = color of y-axis
//*       zcolor             (Color) = color of z-axis
//*       mtx                (float) = initial model translation x-axis
//*       mty                (float) = initial model translation y-axis
//*       mtz                (float) = initial model translation z-axis
//*       modelScale        (double) = scale of model
//*       xlength            (float) = length of satellite, x-axis
//*       ylength            (float) = length of satellite, y-axis
//*       zlength            (float) = length of satellite, z-axis
//*       time_values  (float array) = array of time values
//*       quat_values(float arrayx2) = double array of quaternions
//*       ixx                (float) = satellite principal moments of inertia
//*       iyy                (float) = satellite principal moments of inertia
//*       izz                (float) = satellite principal moments of inertia
//*       simType           (String) = string of simulation type being run
//*       curr_pt              (int) = current point in animation
//*       tot_pts              (int) = total number of points in animation
//*       start_button      (Button) = start animation button
//*       stop_button       (Button) = stop animation button
//*       reset_button      (Button) = reset view axis button
//*       timeLabel          (Label) = label of current time
//*       animSpeedLabel     (Label) = label of "Animation Speed"
//*       plus_button       (Button) = Plus Button (increase animation speed)
//*       minus_button      (Button) = Minus Button (decrease animation speed)
//*       animControl        (Timer) = timer for controlling animation
//*       nf          (NumberFormat) = NumberFormat
//*       modelTrans(TransformGroup) = Java3D transform group of model
//*       sattlTrans(TransformGroup) = Java3D transform group of satellite
//*       iner_axisTrans(TrnsfrmGrp) = Java3D transform group of inertial axis
//*       axislabelTrans(TrnsfrmGrp) = Java3D transform group of inertial axis
//*				     labels
//*	SatTrans     (Transform3D) = Java3D transformation
//*       satellite (..geometry.Box) = Java3D shape (box)
//*       inertial_axes  (LineArray) = line array of points to draw
//*				     inertial axes
//*	body_axes      (LineArray) = line array of points to draw
//*				     body axes on satellite
//*       satelliteApp  (Appearance) = appearance of satellite
//*       GraphicPanel       (Panel) = graphical panel containing 3D scene
//*       CtrlPanel          (Panel) = panel containing animation interface
//*       anim_slider      (JSlider) = slider to show time in animation
//*       delayvalue           (int) = shows delay between frames in animation
//*
//************************************************************************/   


public class AnimationWindow extends JPanel implements ActionListener, ChangeListener
{
   private static AnimationWindow theAnimWindow;
   JFrame animFrame = new JFrame(" ");
      
   // Location of Inertial Axes
   final float i_xaxes = -7.0f;
   final float i_yaxes = -7.0f;
   final float i_zaxes = 0.0f;
   final float axeslength = 2.0f;

   // Assigns the inertial and body axes' colors
   Color3f xcolor = new Color3f(1.0f, 0.0f, 0.0f);
   Color3f ycolor = new Color3f(1.0f, 1.0f, 1.0f);
   Color3f zcolor = new Color3f(0.0f, 0.0f, 1.0f);

   // Initial model translation
   final float mtx = 0f;
   final float mty = 0f;
   final float mtz = -10.0f;

   // Model Constants
   final double modelScale = 0.5d;

   // Lengths of Satellite
   float xlength, ylength, zlength;
     
   // Time values
   float time_values[];
      
   // Quaternion values
   float quat_values[][];

   // Satellite principal moments of inertia
   float ixx, iyy, izz;

   String simType;      

   int curr_pt;				// Current point in animation
   int tot_pts;				// Total number of points in anim

   Button start_button, stop_button, reset_button;
   Label timeLabel, animSpeedLabel;
   Button plus_button, minus_button;
      
   Timer animControl;

   // Set Number of digits right of decimal
  NumberFormat nf = NumberFormat.getNumberInstance();
      
   // Initialize The transformationGroups
   TransformGroup modelTrans;
   TransformGroup sattlTrans;
   TransformGroup iner_axisTrans;
   TransformGroup axislabelTrans;

   Transform3D SatTrans = new Transform3D();
      
   // Initialize Box
   com.sun.j3d.utils.geometry.Box satellite;
      
   // Initialize LineArray for inertial and body axes
   LineArray inertial_axes;
   LineArray body_axes;
      
   // Initialize Appearance attributes
   Appearance satelliteApp;

   // Panel for the buttons
   Panel GraphicPanel;
   Panel CtrlPanel;
      
   // Slider at bottom of window to display animation
   JSlider anim_slider;
      
   // Delay value for animation
   int delayValue;

//**
//*
//*   Description:<br>
//*	This routine creates the pop-up animation window<br>
//*
//*   Inputs:<br>
//*       title             (String) = title that goes in window title bar<br>
//*       ixxt               (float) = temp principal moment of inertia, x-axis<br>
//*       iyyt               (float) = temp principal moment of inertia, y-axis<br>
//*       izzt               (float) = temp princiapl moment of inertia, z-axis<br>
//*       pts                  (int) = number of animation points<br>
//*       tvars      (float arrayx2) = array of values<br>
//*       tsimType          (String) = temp type of simulation<br>
//*				     (e.g. "Gravity Gradient")<br>
//*   Outputs:<br>
//*	None<br>
//*
//*   Local Variables:<br>
//*	frameWidth	     (int) = width of frame<br>
//*	frameHeight	     (int) = height of frame<br>
//*	content	       (Container) = Container of objects<br>
//*       c               (Canvas3D) = Java3D canvas<br>
//*       scene        (BranchGroup) = Java3D branchgroup<br>
//*       u         (SimpleUniverse) = Java3D simple universe<br>
//*
//*/
  /**
   * Creates the pop-up animation window
   * @param	title 		(String) title that goes in window title bar
   * @param	ixxt 		(float) temp principal moment of inertia, x-axis
   * @param	iyyt 		(float) temp principal moment of inertia, y-axis
   * @param	izzt 		(float) temp principal moment of inertia, z-axis
   * @param	pts 		(int) number of animation points
   * @param	tvars 		(float[][])	array of values
   * @param	tsimType 	(String) temp type of simulation (e.g. "Gravity Gradient")
   *
   */
   public AnimationWindow(
             String title     , float ixxt       , float iyyt       ,
             float izzt       , int pts          , float tvars[][]  ,
             String simTypet                                         )
   {
      int frameWidth = 300;
      int frameHeight = 400;

      animFrame.setTitle(title);
      Container content = animFrame.getContentPane();
      content.setLayout(new BorderLayout());
      animFrame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
      animFrame.setBounds(50, 100, frameWidth, frameHeight);
      content.setBackground(Color.white);

      delayValue = 100;

      // Assign passed parameters into frame class
      tot_pts = pts;

      curr_pt = 0;

      nf.setMaximumFractionDigits(1);
     
      // Time Array
      time_values = new float[tot_pts+1];
         
      // Quaternions array
      quat_values = new float[4][tot_pts+1];
         
      // Assign the time and quaternion values
      for (int index = 0; index <= tot_pts; index++)
      {
         time_values[index] = tvars[0][index];
         quat_values[0][index] = tvars[1][index];	// e1
         quat_values[1][index] = tvars[2][index];	// e2
         quat_values[2][index] = tvars[3][index];	// e3
         quat_values[3][index] = tvars[4][index];	// e4
      }

      // Assign the principal moments of inertia
      ixx = ixxt;
      iyy = iyyt;
      izz = izzt;

      // Assigns the simulation type
      simType = simTypet;

      // Creates JSlider with points 0 to total points
      anim_slider = new JSlider(SwingConstants.HORIZONTAL, 0, tot_pts, 0);
      anim_slider.addChangeListener(this);
      // Fills the slider as it is finished
      anim_slider.putClientProperty("JSlider.isFilled", Boolean.TRUE);

      GraphicPanel = new Panel();
      GraphicPanel.setLayout(new BorderLayout());

      GraphicsConfiguration config = SimpleUniverse.getPreferredConfiguration();
      Canvas3D c = new Canvas3D(config);

      // Create a simple scene and attach it to the virtual universe
      BranchGroup scene = createSceneGraph();
      SimpleUniverse u = new SimpleUniverse(c);
      // This will move the ViewPlatform back
      u.getViewingPlatform().setNominalViewingTransform();
      u.addBranchGraph(scene);

      GraphicPanel.add("Center", c);
      content.add("Center", GraphicPanel);

      // Adds a Panel with all the control buttons on it
      CtrlPanel = new Panel();
         
      // *** GridBagLayout Routine ***
      // Sets up the GridBagLayout to add into Panel
      GridBagLayout anim_gbl = new GridBagLayout();
      // Sets up the GridBagConstraints
      GridBagConstraints constraints = new GridBagConstraints();
      CtrlPanel.setLayout(anim_gbl);
    
      // Sets variables to determine which line the fields appear in
      int buttonline = 0;			// Button Lines
      int scrolline = 1;			// Scrollbar Line
      int animline = 2;				// Animation Speed Line

      stop_button = new Button("Stop");
      buildConstraints(constraints, 0, buttonline, 1, 1);
      anim_gbl.setConstraints(stop_button, constraints);
      CtrlPanel.add(stop_button);
      stop_button.addActionListener(this);

      start_button = new Button("Start");
      buildConstraints(constraints, 1, buttonline, 1, 1);
      anim_gbl.setConstraints(start_button, constraints);
      CtrlPanel.add(start_button);
      start_button.addActionListener(this);

      if (simType == "Gravity Gradient")
         timeLabel = new Label("Time: " + nf.format(time_values[0]) + " orbits  ");
      else
         timeLabel = new Label("Time: " + nf.format(time_values[0]) + " secs   ");
      buildConstraints(constraints, 2, buttonline, 3, 1);
      anim_gbl.setConstraints(timeLabel, constraints);
      CtrlPanel.add(timeLabel);

      reset_button = new Button("Reset Axes");
      buildConstraints(constraints, 5, buttonline, 2, 1);
      anim_gbl.setConstraints(reset_button, constraints);
      CtrlPanel.add(reset_button);
      reset_button.addActionListener(this);

      buildConstraints(constraints, 0, scrolline, 7, 1);
      anim_gbl.setConstraints(anim_slider, constraints);
      CtrlPanel.add(anim_slider);

      minus_button = new Button("-");
      buildConstraints(constraints, 1, animline, 1, 1);
      anim_gbl.setConstraints(minus_button, constraints);
      CtrlPanel.add(minus_button);
      minus_button.addActionListener(this);
         
      animSpeedLabel = new Label("Animation Speed");
      buildConstraints(constraints, 2, animline, 3, 1);
      anim_gbl.setConstraints(animSpeedLabel, constraints);
      CtrlPanel.add(animSpeedLabel);
         
      plus_button = new Button("+");
      buildConstraints(constraints, 5, animline, 1, 1);
      anim_gbl.setConstraints(plus_button, constraints);
      CtrlPanel.add(plus_button);
      plus_button.addActionListener(this);         
      // Adds Panel to BorderLayout
      content.add("South", CtrlPanel);

      // Parameters for animation control
      // 100 = 1 sec
      animControl = new Timer(delayValue, this);
      animControl.start();
         
      animFrame.setVisible(true);
   } // End AnimationWindow constructor


//************************************************************************
//*
//*   buildConstraints
//*
//*   Description:
//*	This sets the constraints for a cell in the Grid Bag Layout
//*
//*   Inputs:
//*       gbc   (GridBagConstraints) = Grid Bag Constraints
//*       gx                   (int) = grid x-location
//*       gy                   (int) = grid y-location
//*       gw                   (int) = grid width
//*       gh                   (int) = grid height
//*
//*   Outputs:
//*	None
//*
//*   Local Variables:
//*	None
//*
//************************************************************************/
/**
 * Sets the constraints for a cell in the Grid Bag Layout
 * @param	gbc (GridBangConstraints)	Grid Bang Constraints
 * @param	gx	(int)					grid x-location
 * @param	gy	(int)					grid y-location
 * @param	gw	(int)					grid width
 * @param	gh	(int)					grid height
 */
   void buildConstraints(GridBagConstraints gbc, int gx, int gy, int gw, int gh)
   {
      // Location of field in the layout
      gbc.gridx = gx;
      gbc.gridy = gy;
      // Number of cell widths and heights field requires
      gbc.gridwidth = gw;
      gbc.gridheight = gh;
   }


//************************************************************************
//*
//*   stateChanged
//*
//*   Description:
//*	This routine handles events caused by the JSlider
//*
//*   Inputs:
//*       e            (ChangeEvent) = ChangeEvent
//*
//*   Outputs:
//*	None
//*
//*   Local Variables:
//*	new_point            (int) = new point in slider with changed loc
//*       jslider1         (JSlider) = temp slider with changed location
//*
//*
//************************************************************************/
  /**
   * Handles events caused by the JSlider
   * @param	e (ChangeEvent)	ChangeEvent)
   */
   public void stateChanged(ChangeEvent e)
   {
      int new_point;
      JSlider jslider1 = (JSlider) e.getSource();
      
      // Find the new point of animation if slider is used
      new_point = jslider1.getValue();
      // Assign the current point to the new point
      curr_pt = new_point;
         
      // Stop the animation if the current point is at the end
      if (curr_pt == tot_pts)
         animControl.stop();

      // Updates the control panel and graphic if animation is not running
      if (animControl.isRunning() == false)
      {
         animControl.start();
         SatTrans.setRotation(new
            Quat4f(quat_values[0][curr_pt],quat_values[1][curr_pt], 
            quat_values[2][curr_pt], quat_values[3][curr_pt]));
         sattlTrans.setTransform(SatTrans);
         if (simType == "Gravity Gradient")
            timeLabel.setText("Time: " + nf.format(time_values[curr_pt]) + " orbits   ");
         else
            timeLabel.setText("Time: " + nf.format(time_values[curr_pt]) + " secs   ");
         animControl.stop();
      }

   }


//************************************************************************
//*
//*   actionPerformed
//*
//*   Description:
//*	This routine handles any button presses and the animation
//*       of the satellite by drawing the next frame and updating the
//*       slider
//*
//*   Inputs:
//*       e            (ActionEvent) = ActionEvent
//*
//*   Outputs:
//*	None
//*
//*   Local Variables:
//*       modTrans     (Transform3D) = Java3D transformation
//*
//************************************************************************/
  /**
   * Handles any button presses and the animation of the satellite
   * by drawing the next frame and updating the slider
   * 
   * @param	e (ActionEvent)	ActionEvent
   */
   public void actionPerformed(ActionEvent e)
   {
      // Transform3D for satellite and model
      Transform3D modTrans = new Transform3D();
       
      if ((e.getSource() == start_button) && (curr_pt >= tot_pts))
      {
         curr_pt = 0;
         animControl.start();
      }
      else if (e.getSource() == start_button)
         animControl.start();
      else if (e.getSource() == stop_button)
         animControl.stop();
      else if (e.getSource() == reset_button)
      {
         modTrans.setTranslation(new Vector3f(mtx, mty, mtz));
         modTrans.setRotation(new Quat4f(0.0f, 0.0f, 0.0f, 1.0f));
         modTrans.setScale(new Vector3d(modelScale, modelScale, modelScale));
         modelTrans.setTransform(modTrans);
      }
      else if (e.getSource() == minus_button)
      {
         if (delayValue < 1000)
            delayValue = delayValue * 5;
         animControl.setDelay(delayValue);
      }
      else if (e.getSource() == plus_button)
      {
         if (delayValue > 5)
            delayValue = delayValue / 5;
         animControl.setDelay(delayValue);
      }

      else	// Run animation
      {
         SatTrans.setRotation(new
            Quat4f(quat_values[0][curr_pt], quat_values[1][curr_pt],
            quat_values[2][curr_pt], quat_values[3][curr_pt]));
         sattlTrans.setTransform(SatTrans);
      // Increment to next point
         curr_pt = curr_pt + 1;         
      }

      // Moves the slider as the animation is changed
      anim_slider.setValue(curr_pt);
         
      // Updates the timeLabel
      if (simType == "Gravity Gradient")
         timeLabel.setText("Time: " + nf.format(time_values[curr_pt]) + " orbits   ");
      else
         timeLabel.setText("Time: " + nf.format(time_values[curr_pt]) + " secs   ");   
      // Stop animation at end                  
      if (curr_pt >= tot_pts)
         animControl.stop();
   } // End actionPerformed



//************************************************************************
//*
//*   createSceneGraph
//*
//*   Description:
//*	This routine creates the scene graph of the universe, populates
//*       it with lights, the initial camera location, and allows
//*       mouse interaction with the scene
//*
//*   Inputs:
//*	None
//*
//*   Outputs:
//*	objRoot      (BranchGroup) = Java3D BranchGroup
//*
//*   Local Variables:
//*	objRoot      (BranchGroup) = Java3D BranchGroup
//*       bounds    (BoundingSphere) = Java3D BoundingSphere
//*       light   (DirectionalLight) = Java3D light source
//*       trans        (Transform3D) = Java3D transformation
//*
//************************************************************************/
   /**
    * Create the scene graph of the universe, populates it with lights,
    * the initial camera location, and allows mouse interaction with
    * the scene
    */
   public BranchGroup createSceneGraph()
   {
      // Create the root of the branch graph
      BranchGroup objRoot = new BranchGroup();
   
      // Create light source
      // first bound its influence
      BoundingSphere bounds = new BoundingSphere(new Point3d(0.0, 0.0, 0.0), 100.0);
      
      // now make the light itself
      DirectionalLight light = new DirectionalLight(new Color3f(1.0f, 1.0f, 1.0f),
         new Vector3f(0f, 0f, -1f));
         
      light.setInfluencingBounds(bounds);
      objRoot.addChild(light);
         
      // TransformGroup that moves the whole model
      modelTrans = new TransformGroup();
         
      // Transform3D
      Transform3D trans = new Transform3D();
      // Initially move it away from camera
      trans.setTranslation(new Vector3f(mtx, mty, mtz));
      // and shrink everything so it fits in field of view
      trans.setScale(new Vector3d(modelScale, modelScale, modelScale));
      modelTrans.setTransform(trans);
         
      // Allows mouse to alter scene
      modelTrans.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
      modelTrans.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
      
      // Create Mouse Rotate
      MouseRotate modelRot = new MouseRotate();
      modelRot.setTransformGroup(modelTrans);
      modelTrans.addChild(modelRot);
      modelRot.setSchedulingBounds(bounds);
         
      // Create the Zoom
      MouseZoom modelZoom = new MouseZoom();
      modelZoom.setTransformGroup(modelTrans);
      modelTrans.addChild(modelZoom);
      modelZoom.setSchedulingBounds(bounds);
         
      // Create the translate
      MouseTranslate modelTranslate = new MouseTranslate();
      modelTranslate.setTransformGroup(modelTrans);
      modelTrans.addChild(modelTranslate);
      modelTranslate.setSchedulingBounds(bounds);
         
      modelTrans.addChild(generateAxes());
      modelTrans.addChild(generateAxesLabels());
      modelTrans.addChild(generateSatellite());
         
      objRoot.addChild(modelTrans);
      return objRoot;
   } // End createSceneGraph
      

//************************************************************************
//*
//*   generateAxes
//*
//*   Description:
//*	This routine creates the inertial reference axis in the
//*       animation
//*
//*   Inputs:
//*	None
//*
//*   Outputs:
//*	iner_axisTrans  (TransformGroup) = Java3D TransformGroup
//*
//*   Local Variables:
//*	iner_axisTrans  (TransformGroup) = Java3D TransformGroup
//*       origin           (Point3f) = origin point
//*       trans        (Transform3D) = Java3D transformation
//*
//************************************************************************/
  /**
   * Creates the inertial reference axis in the animation
   */
   public TransformGroup generateAxes()
   {
      Point3f origin = new Point3f(0.0f, 0.0f, 0.0f);

      iner_axisTrans = new TransformGroup();
      
      Transform3D trans = new Transform3D();
      trans.setTranslation(new Vector3f(i_xaxes, i_yaxes, i_zaxes));
      iner_axisTrans.setTransform(trans);
      
      inertial_axes = new LineArray(6, LineArray.COORDINATES | LineArray.COLOR_3);
      
      // X-Axis
      inertial_axes.setCoordinate(0, origin);
      inertial_axes.setCoordinate(1, new Point3f(axeslength, 0.0f, 0.0f));
      inertial_axes.setColor(0, xcolor);
      inertial_axes.setColor(1, xcolor);

      // Y-Axis
      inertial_axes.setCoordinate(2, origin);
      inertial_axes.setCoordinate(3, new Point3f(0.0f, axeslength, 0.0f));
      inertial_axes.setColor(2, ycolor);
      inertial_axes.setColor(3, ycolor);

      // Z-Axis
      inertial_axes.setCoordinate(4, origin);
      inertial_axes.setCoordinate(5, new Point3f(0.0f, 0.0f, axeslength));
      inertial_axes.setColor(4, zcolor);
      inertial_axes.setColor(5, zcolor);

      iner_axisTrans.addChild(new Shape3D(inertial_axes));
      return iner_axisTrans;
   } // End generateAxes


//************************************************************************
//*
//*   generateAxesLabels
//*
//*   Description:
//*	This routine creates the labels for the inertial reference
//*       axis in the animation (The white "X", "Y", and "Z")
//*
//*   Inputs:
//*	None
//*
//*   Outputs:
//*	axislabelTrans  (TransformGroup) = Java3D TransformGroup
//*
//*   Local Variables:
//*	axislabelTrans  (TransformGroup) = Java3D TransformGroup
//*       font3d            (Font3D) = Java3D Font3D
//*       xfont             (Text3D) = 3D "X"
//*       xshape           (Shape3D) = Java3D Shape3D of xfont
//*       yfont             (Text3D) = 3D "Y"
//*       yshape           (Shape3D) = Java3D Shape3D of yfont
//*       zfont             (Text3D) = 3D "Z"
//*       zshape           (Shape3D) = Java3D Shape3D of zfont
//*
//************************************************************************/
  /**
   * Creates the labels for the inertial reference axis in the animation
   * (The white "X", "Y"m and "Z")
   */
   public TransformGroup generateAxesLabels()
   {
      axislabelTrans = new TransformGroup();
    
      Font3D font3d = new Font3D(new Font("Display", Font.PLAIN, 1),
         new FontExtrusion());
         
      // X-Axis Label
      Text3D xfont = new Text3D(font3d, new String("X"),
         new Point3f((i_xaxes + axeslength + 0.5f), (i_yaxes - 0.25f), 0.0f));
      Shape3D xshape = new Shape3D(xfont);
         
      // Y-Axis Label
      Text3D yfont = new Text3D(font3d, new String("Y"),
         new Point3f((i_xaxes - 0.25f), (i_yaxes + axeslength + 0.5f), 0.0f));
      Shape3D yshape = new Shape3D(yfont);
         
      // Z-Axis Label
      Text3D zfont = new Text3D(font3d, new String("Z"),
         new Point3f((i_xaxes-0.5f), (i_yaxes-0.5f), (axeslength+0.5f)));
      Shape3D zshape = new Shape3D(zfont);
         
      axislabelTrans.addChild(xshape);
      axislabelTrans.addChild(yshape);
      axislabelTrans.addChild(zshape);
         
      return axislabelTrans;
   } // End generateAxesLabels


//************************************************************************
//*
//*   generateSatellite
//*
//*   Description:
//*	This routine creates the box shape of the satellite using the
//*	principal moments of inertia.
//*
//*   Inputs:
//*	None
//*
//*   Outputs:
//*	sattlTrans  (TransformGroup) = Java3D TransformGroup
//*
//*   Local Variables:
//*	sattlTrans  (TransformGroup) = Java3D TransformGroup
//*       black            (Color3f) = color black
//*       origin           (Point3f) = origin point
//*       trans        (Transform3D) = Java3D Transform3D
//*
//************************************************************************/ 
  /**
   * Create the box shape of the satellite using the principal 
   * moments of inertia
   */
   public TransformGroup generateSatellite()
   {
      // Assigns colors
      Color3f black = new Color3f(0.0f, 0.0f, 0.0f);         

      Point3f origin = new Point3f(0.0f, 0.0f, 0.0f);

      Quat4f init_quats = new Quat4f(
         quat_values[0][0], quat_values[1][0],
         quat_values[2][0], quat_values[3][0]);
      sattlTrans = new TransformGroup();
      sattlTrans.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
    
      Transform3D trans = new Transform3D();
      trans.setTranslation(new Vector3d(0.0d, 0.0d, 0.0d));
      // Set initial Quaternion Values
      trans.setRotation(init_quats);
      sattlTrans.setTransform(trans);

      satelliteApp = new Appearance();
      satelliteApp.setCapability(Appearance.ALLOW_MATERIAL_WRITE);
      satelliteApp.setMaterial(new Material(black, black, new Color3f(1.0f, 0.0f, 0.0f), black, 0f));
            
      // Find the satellite dimensions
      findDimensions();
         
      satellite = new com.sun.j3d.utils.geometry.Box(xlength, ylength, zlength, satelliteApp);
      

      // Places body axes on satellite
      body_axes = new LineArray(6, LineArray.COORDINATES | LineArray.COLOR_3);
    
      // X-Axis
      body_axes.setCoordinate(0, origin);
      body_axes.setCoordinate(1, new Point3f(xlength + 1.5f, 0.0f, 0.0f));
      body_axes.setColor(0, xcolor);
      body_axes.setColor(1, xcolor);

      // Y-Axis
      body_axes.setCoordinate(2, origin);
      body_axes.setCoordinate(3, new Point3f(0.0f, ylength + 1.5f, 0.0f));
      body_axes.setColor(2, ycolor);
      body_axes.setColor(3, ycolor);

      // Z-Axis
      body_axes.setCoordinate(4, origin);
      body_axes.setCoordinate(5, new Point3f(0.0f, 0.0f, zlength+1.5f));
      body_axes.setColor(4, zcolor);
      body_axes.setColor(5, zcolor);

      sattlTrans.addChild(new Shape3D(body_axes));
      sattlTrans.addChild(satellite);
      return sattlTrans;

   } // End generateSatellite
      



//************************************************************************
//*
//*   findDimensions
//*
//*   Description:
//*	This routine finds the normalized lengths of the satellite
//*       based on the principal moments of inertia
//*
//*   Inputs:
//*	None
//*
//*   Outputs:
//*	None
//*
//*   Local Variables:
//*       x_temp             (float) = temp length along x-axis
//*       y_temp             (float) = temp length along y-axis
//*       y_temp             (float) = temp length along y-axis
//*       max_length         (float) = maximum length of temp lengths
//*       largest_length     (float) = length of largest side
//*
//************************************************************************/ 
  /**
   * Finds the normalized lengths of the satellite based on the 
   * principal moments of inertia
   */
   public void findDimensions()
   {
      float x_temp, y_temp, z_temp;
      float max_length;
      final float largest_length;

      largest_length = 3f;

      // Assuming that the satellite is rectangular-shaped
      // and homogeneous, then the equations for inertias are:
      // Ixx = (m/12)*(b^2 + c^2)
      // Iyy = (m/12)*(a^2 + c^2)
      // Izz = (m/12)*(a^2 + b^2)
      // drop the (m/12) and we can solve for a, b, and c:
      // a^2 = Iyy - 0.5*(Ixx + Iyy - Izz)
      // b^2 = Izz - Iyy + 0.5*(Ixx + Iyy - Izz)
      // c^2 = 0.5*(Izz + Iyy - Izz)

      x_temp = (float)Math.sqrt((double)(iyy - 0.5f*(ixx + iyy - izz)));
      y_temp = (float)Math.sqrt((double)(izz - iyy + 0.5f*(ixx + iyy - izz)));
      z_temp = (float)Math.sqrt((double)(0.5f*(izz + iyy - izz)));

      // Find the maximum between a_temp, b_temp, and c_temp
      max_length = Math.max(x_temp, y_temp);
      max_length = Math.max(max_length, z_temp);

      // Set max_length so that you will get 
      max_length = max_length / largest_length;

      // Divide each length by max_length for the largest length
      xlength = x_temp / max_length;
      ylength = y_temp / max_length;
      zlength = z_temp / max_length;
   } // End findDimensions

} // End AnimationWindow