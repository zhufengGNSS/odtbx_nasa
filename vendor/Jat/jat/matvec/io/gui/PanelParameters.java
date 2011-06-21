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

package jat.matvec.io.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

/**
 * <p>Titre : JAva MAtrix TOols</p>
 * <p>Description : builds a JPanel containing fields for setting parameters.</p>
 * @author Yann RICHET
 */

public class PanelParameters extends JPanel implements FocusListener {

  private Dimension defaultSize ;

  private String[] paramLabels;
  private String[] paramValues;
  private String[][] paramChoices;

  private JLabel[] labels;
  private JComboBox[] fields;

  public PanelParameters(String[] lab, String[] val) {
    paramLabels = lab;
    paramValues = val;
    paramChoices = new String[paramLabels.length][1];
    for (int i = 0; i < paramLabels.length; i++) {
      paramChoices[i][0] = paramValues[i];
    }

    setComponents();
    setAppearence();
    draw();
  }

  public PanelParameters(String[] lab, String[][] ch) {
    paramLabels = lab;
    paramValues = new String[paramLabels.length];
    paramChoices = ch;
    for (int i = 0; i < paramLabels.length; i++) {
      paramValues[i] = paramChoices[i][0];
    }

    setComponents();
    setAppearence();
    draw();
  }

  public PanelParameters(String[] lab) {
    paramLabels = lab;
    paramValues = new String[paramLabels.length];
    paramChoices = new String[paramLabels.length][1];

    setComponents();
    setAppearence();
    draw();
  }

  private void setComponents() {
    labels = new JLabel[paramLabels.length];
    fields = new JComboBox[paramLabels.length];
    for (int i = 0 ; i < paramLabels.length; i++) {
      labels[i] = new JLabel(paramLabels[i],JLabel.RIGHT);
      fields[i] = new JComboBox(paramChoices[i]);
      fields[i].setEditable(true);
    }
    defaultSize = new Dimension(400,paramLabels.length*30);
  }

  private void setAppearence() {
    setPreferredSize(defaultSize);
    setSize(defaultSize);
  }

  private void updateValues() {
    for (int i=0;i<paramLabels.length;i++) {
      paramValues[i] = (String)(fields[i].getSelectedItem());
    }
  }

  public void focusLost(FocusEvent e) {
    updateValues();
  }

  public void focusGained(FocusEvent e) {}

  public String[] getValues() {
    updateValues();
    return paramValues;
  }

  private void buildConstraints(GridBagConstraints gbc,int gx, int gy, int gw, int gh, int wx, int wy) {
    gbc.gridx = gx;
    gbc.gridy = gy;
    gbc.gridwidth = gw;
    gbc.gridheight = gh;
    gbc.weightx = wx;
    gbc.weighty = wy;
  }

  private void draw() {
    JPanel panel = new JPanel();

    GridBagLayout gbl = new GridBagLayout();
    GridBagConstraints c = new GridBagConstraints();
    panel.setLayout(gbl);

    for (int i = 0; i < paramLabels.length; i++) {
      fields[i].addFocusListener(this);

      // Ajout du panel de la chaine
      buildConstraints(c,0,i,1,1,50,20);
      c.anchor = GridBagConstraints.EAST;
      gbl.setConstraints(labels[i],c);
      panel.add(labels[i]);

      // Ajout du panel de la chaine
      buildConstraints(c,1,i,1,1,50,20);
      c.fill = GridBagConstraints.HORIZONTAL;
      gbl.setConstraints(fields[i],c);
      panel.add(fields[i]);
    }

    JScrollPane scrollPane = new JScrollPane(panel);

    scrollPane.setPreferredSize(getSize());
    scrollPane.setSize(getSize());

    this.setLayout(new BorderLayout());
    this.add(scrollPane,BorderLayout.CENTER);

  }

}
