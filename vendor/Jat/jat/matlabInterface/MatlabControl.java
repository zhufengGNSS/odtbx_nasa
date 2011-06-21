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
*/



/**
* This class interfaces with the current Matlab session, allowing you
* to call matlab commands from Java objects
*
* @author <a href="mailto:kamin@cs.berkeley.edu">Kamin Whitehouse</a>
*/

//package net.tinyos.matlab;

package jat.matlabInterface;

import com.mathworks.jmi.*;

public class MatlabControl {
   Matlab matlab = null; //this is the com.mathworks.jmi.Matlab class,which has functionality allowing one to interact with the matlab session.
   boolean useCb=false;
   Object returnVal;
   String callbackFunction; 

   /***************CONSTRUCTORS****************/
   /*** usually, the default constructor with no arguments is fine.
    * Sometimes, a callback function is useful.  A callback function
    * is allows the user to pass the command and its arguments to a
    * matlab function, which will figure out the best way to executue
    * it, instead of trying to execute it directly from java.  This
    * is often useful for handling/printing errors as well as dealing
    * with native matlab data types, like cell arrays.  The cell
    * array is not really converted to a java type, only the handle
    * of the cell array is converted to a java type, so often a
    * matlab callback function is useful for dealing specially with
    * cell arrays, etc.
    ***/

   public MatlabControl() {
       this(false);
   }

   public MatlabControl(boolean useCallback) {
	this(useCallback,new String("matlabControlcb"));
   }

   public MatlabControl(boolean useCallback, String CallBackFunction) {
       try {
           if (matlab == null)
               matlab = new Matlab();//this command links to the current matlab session
       } catch (Exception e) {
	    System.out.println(e.toString());
       }
       returnVal = new String("noReturnValYet");
       this.useCb=useCallback;
	callbackFunction=CallBackFunction;
   }

   /***************USER-LEVEL FUNCTIONS****************/
   /***call these from any java thread (as long as the thread
    * originated from within matlab) 
    ***/

   /**Evaluate a string, Matlab script, or Matlab function**/
   public void eval(String Command) {
       Matlab.whenMatlabReady(new MatlabEvalCommand(Command,useCb));
   }

   /**Evaluate a Matlab function that requires arguments.  Each element of
    the "args" vector is an argument to the function "Command"**/
   public void feval(String Command, Object[] args) {
       Matlab.whenMatlabReady(new MatlabFevalCommand(Command, args,useCb));
   }

   /**Evaluate a Matlab function that requires arguments and provide return arg.
    * Each element of the "args" vector is an argument to the function "Command"**/
   public Object blockingFeval(String Command, Object[] args) throws InterruptedException {
       returnVal = new String("noReturnValYet");
       Matlab.whenMatlabReady(new MatlabBlockingFevalCommand(Command, args, useCb, this));
       if (returnVal.equals(new String("noReturnValYet"))) {
           synchronized(returnVal){
                returnVal.wait();
           }
       }
       return returnVal;
   }

   /**Echoing the eval statement is useful if you want to see in
    * matlab each time that a java function tries to execute a matlab
    * command **/
   public void setEchoEval(boolean echo){
	Matlab.setEchoEval(echo);
   }

   /**********TEST FUNCTIONS***********************/
   /***call these functions from within Matlab itself.  These are examples of the general execution order:
	1. instantiate java object from matlab
	2. spawn a new Java thread
	3. call matlab functions from new java thread  

	EXAMPLE (from matlab prompt):

	>> mc=MatlabControl;
	>> mc.testEval('x = 5')
	x =
	       5
	> mc.testFeval('help',{'sqrt'})
	  SQRT    Square root.
	     SQRT(X) is the square root of the elements of X. Complex
	     results are produced if X is not positive.
	
	     See also SQRTM.
	
	     Overloaded methods
	     help sym/sqrt.m
	
	>> mc.testBlockingFeval('sqrt',{x})
	       2.2361
   ****/

   public void testEval(final String Command) {
   class Caller extends Thread{
        public void run(){
             try{
                  eval(Command);
             } catch(Exception e){
		  //                System.out.println(e.toString());
	      }
        }
   }
   Caller c = new Caller();
   c.start();
   }

   public void testFeval(final String Command, final Object[] args) {
   class Caller extends Thread{
        public void run(){
             try{
                  feval(Command, args);
             } catch(Exception e){
		  //                System.out.println(e.toString());
	      }
        }
   }
   Caller c = new Caller();
   c.start();
   }

   public void testBlockingFeval(final String Command, final Object[] args) {
   class Caller extends Thread{
        public void run(){
             try{
                  Object rets[] = new Object[1];
                  rets[0] = blockingFeval(Command, args);
                  feval(new String("disp"),rets );
             } catch(Exception e){
		  //                System.out.println(e.toString());
	      }
        }
   }
   Caller c = new Caller();
   c.start();
   }



   /********   INTERNAL FUNCTIONS AND CLASSES  *******/

   public void setReturnVal(Object val) {
       synchronized(returnVal){
           Object oldVal = returnVal;
           returnVal = val;
           oldVal.notifyAll();
       }
   }

   /** This class is used to execute a string in Matlab **/
   protected class MatlabEvalCommand implements Runnable {
       String command;
       boolean useCallback,eval;
       Object[] args;

       public MatlabEvalCommand(String Command, boolean useCallback) {
           command = Command;
           this.useCallback = useCallback;
	    eval=true;
           args=null;
       }

       protected Object useMatlabCommandCallback(String command, Object[] args){
           int numArgs = (args==null)? 0 : args.length;
           Object newArgs[] = new Object[numArgs+1] ;
           newArgs[0]=command;
           for(int i=0;i<numArgs;i++){
               newArgs[i+1] = args[i];
           }
           try{
               return matlab.mtFevalConsoleOutput(new String("matlabControlcb"), newArgs, 0);
           }
           catch(Exception e){
		//                System.out.println(e.toString());
		return null;
	    }
       }

       public void run() {
           try {
               if(useCallback){
                   useMatlabCommandCallback(command, args);
               }
               else if(eval){
                   matlab.evalConsoleOutput(command);
               }
               else{
                   matlab.fevalConsoleOutput(command, args);
                   // For Matlab R13 use the line below instead:
                   // matlab.fevalConsoleOutput(command, args, 0, null);
                   }
           } catch (Exception e) {
               System.out.println(e.toString());
           }
       }
   }

   /** This class is used to execute a function in matlab and pass paramseters**/
   protected class MatlabFevalCommand extends MatlabEvalCommand {

       public MatlabFevalCommand(String Command, Object[] Args, boolean useCallback) {
           super(Command, useCallback);
           args = Args;
	    eval=false;
       }

   }

   /** This class is used to execute a function in matlab and pass paramseters
    *  and it also return arguments**/
   protected class MatlabBlockingFevalCommand extends MatlabFevalCommand {
       MatlabControl parent;

       public MatlabBlockingFevalCommand(String Command, Object[] Args, boolean useCallback, MatlabControl parent) {
           super(Command, Args, useCallback);
           this.parent=parent;
       }

       public void run() {
           try {
               if(useCallback){
                   parent.setReturnVal(useMatlabCommandCallback(command, args));
               }
               else{
                   parent.setReturnVal(matlab.mtFevalConsoleOutput(command, args, 0));
               }
           } catch (Exception e) {
		//                System.out.println(e.toString());
           }
       }
   }
}
