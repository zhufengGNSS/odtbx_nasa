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

import jat.matvec.data.Matrix;
import java.lang.Double;

import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.JScrollPane;

//import javax.swing.*;
import java.awt.*;
//import java.awt.event.*;

import javax.swing.table.*;
//import javax.swing.event.TableModelListener;
//import javax.swing.event.TableModelEvent;

public class MatrixTable extends JPanel {

  private JTable table;
  private DoubleModel model;
  private boolean modificationEnabled = false;

  private Dimension defaultSize = new Dimension(400,400);


  public MatrixTable(Matrix m) {
    setModel(m);
    setAppearence();
    toWindow();
  }

  public MatrixTable(double[][] d) {
    setModel(new Matrix(d));
    setAppearence();
    toWindow();
  }

  private void setAppearence() {
    setPreferredSize(defaultSize);
    setSize(defaultSize);
  }


  private void setModel(Matrix m) {
    model = new DoubleModel(m,modificationEnabled);
  }

  private void toWindow() {
    table = new JTable(model);

    table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
    table.setCellSelectionEnabled(true);
    //table.setRowSelectionAllowed(false);
    //table.setColumnSelectionAllowed(false);

    JPanel panel = new JPanel();
    panel.add(table);

    JScrollPane scrollPane = new JScrollPane(panel);

    scrollPane.setPreferredSize(getSize());
    scrollPane.setSize(getSize());

    this.setLayout(new BorderLayout());
    this.add(scrollPane,BorderLayout.CENTER);
  }

  public void update(Matrix m) {
    model.setMatrix(m);
  }

  public Matrix getMarix() {
    return model.getMatrix();
  }

  public void setModificationEnabled() {
    modificationEnabled = true;
  }

  public void setModificationDisabled() {
    modificationEnabled = false;
  }

  private class DoubleModel extends AbstractTableModel {

    private Matrix M;
    private boolean modificationEnabled;

    public int getRowCount() {
      return M.getRowDimension();
    }

    public int getColumnCount() {
      return M.getColumnDimension();
    }

    public DoubleModel(Matrix m,boolean mE) {
      this.setMatrix(m);
      modificationEnabled = mE;
    }

    public void setValueAt(int i, int j,double v) {
        M.set(i,j,v);
    }

    public void setValueAt(Object o,int i, int j) {
        M.set(i,j,(Double.parseDouble((String)o)));
        fireTableCellUpdated(i, j);
    }

     public void setMatrix(Matrix m) {
        M = m;
    }

    public Matrix getMatrix() {
        return M;
    }

    public Object getValueAt(int i, int j) {
        return new Double(M.get(i,j));
    }

    public boolean isCellEditable(int i, int j) {
      return modificationEnabled;
    }

  }
}
