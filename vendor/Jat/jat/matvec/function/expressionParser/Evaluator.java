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

package jat.matvec.function.expressionParser;

import java.util.HashMap;
import jat.matvec.data.Matrix;
import jat.matvec.data.RandomMatrix;
import java.lang.Double;

/************************************************************************
 * <i>Mathematic expression evaluator.</i> Supports the following functions:
 * +, -, *, /, ^, %, cos, sin, tan, acos, asin, atan, sqrt, sqr, log, min, max, ceil, floor, abs, neg, rndr.<br>
 * When the getValue() is called, a Double object is returned. If it returns null, an error occured.<p>
 * <pre>
 * Sample:
 * MathEvaluator m = new MathEvaluator("-5-6/(-2) + sqr(15+x)");
 * m.addVariable("x", 15.1d);
 * System.out.println( m.getValue() );
 * </pre>
 * @version 1.1
 * @author 	The-Son LAI, <a href="mailto:Lts@writeme.com">Lts@writeme.com</a>
 * @date     April 2001
 ************************************************************************/
public class Evaluator
{
   	protected static 	Operator[] 	operators 	= null;
	private 			Node 		node       	= null;
	private 			String  	expression 	= null;
    private 			HashMap	 	variables  	= new HashMap();

	/***
     * creates an empty MathEvaluator. You need to use setExpression(String s) to assign a math expression string to it.
     */
	public Evaluator()
	{
    	init();
	}

    /***
     * creates a MathEvaluator and assign the math expression string.
     */
	public Evaluator(String s)
	{
    	init();
		setExpression(s);
	}

    private void init()
    {
       	if ( operators == null ) initializeOperators();
    }

    /***
     * adds a variable and its value in the MathEvaluator
     */
    public void addVariable(String v, Matrix val)
    {
		variables.put(v, val);
    }

    /***
     * adds a variable and its value in the MathEvaluator
     */
    public void addVariable(String v, double val)
    {
		variables.put(v, new Double(val));
    }

    /***
     * adds a variable and its value in the MathEvaluator
     */
    public void addVariable(String v, Double val)
    {
		variables.put(v, val);
    }

    /***
     * sets the expression
     */
    public void setExpression(String s)
    {
	    expression = s;
    }

    /***
     * resets the evaluator
     */
    public void reset()
    {
    	node 		= null;
        expression 	= null;
        variables 	= new HashMap();
    }

    /***
     * trace the binary tree for debug
     */
	public void trace()
    {
    	try
        {
            node = new Node(expression);
            node.trace();
        }
        catch (Exception e)
        {
        	e.printStackTrace();
        }
    }

    /***
     * evaluates and returns the value of the expression
     */
    public Object getValue()
    {
       	if (expression == null) return null;

    	try
        {
            node = new Node(expression);
            return (evaluate(node));
        }
        catch (Exception e)
        {
        	e.printStackTrace();
            return null;
        }
    }

    private static Object evaluate(Node n)
    {
        if ( n.hasOperator() && n.hasChild() )
        {
            if ( n.getOperator().getType() == 1 )
                n.setValue ( evaluateExpression1( n.getOperator(), evaluate( n.getLeft() )) );
            else if ( n.getOperator().getType() == 2 )
                n.setValue( evaluateExpression2( n.getOperator(), evaluate( n.getLeft() ), evaluate( n.getRight() ) ) );
        }
        return (n.getValue());
    }

    private static Object evaluateExpression1(Operator o, Object o1)
    {

        String op 	= o.getOperator();
        Object res 	= null;

        boolean o1isMatrix;
        boolean o2isMatrix;

        try {
          Matrix m1 = (Matrix)o1;
          o1isMatrix = true;
        } catch (Exception e) {
          o1isMatrix = false;
        }

        if (o1isMatrix) {
          Matrix f1 = (Matrix)o1;

            if       ( "t".equals(op) )  	    res = (Object)(f1.transpose()) ;
            else if  ( "inv".equals(op) )  	  res = (Object)(f1.inverse()) ;
            else if  ( "diag".equals(op) )  	res = (Object)(f1.diag()) ;
            else if  ( "det".equals(op) )  	  res = (Object)(new Double((double)(f1.det()))) ;
            else if  ( "trace".equals(op) )  	res = (Object)(new Double((double)(f1.trace()))) ;
            else if  ( "rank".equals(op) )  	res = (Object)(new Double((double)(f1.rank()))) ;
            else if  ( "sum".equals(op) )  	  res = (Object)(f1.sum()) ;
            else if  ( "prod".equals(op) )  	res = (Object)(f1.prod()) ;
            else if  ( "min".equals(op) )  	  res = (Object)(f1.min()) ;
            else if  ( "max".equals(op) )  	  res = (Object)(f1.max()) ;
            else if  ( "mean".equals(op) )  	res = (Object)(((RandomMatrix)(f1)).mean()) ;
            else if  ( "cov".equals(op) )  	  res = (Object)(((RandomMatrix)(f1)).cov()) ;
            else if  ( "var".equals(op) )  	  res = (Object)(((RandomMatrix)(f1)).var()) ;
            else if  ( "cor".equals(op) )  	  res = (Object)(((RandomMatrix)(f1)).cor()) ;

        } else {
          double f1 = ((Double)o1).doubleValue();

            if       ( "cos".equals(op) )  	  res = (Object)(new Double(Math.cos(f1))) ;
            else if  ( "sin".equals(op) )  	  res = (Object)(new Double(Math.sin(f1))) ;
            else if  ( "tan".equals(op) )  	  res = (Object)(new Double(Math.tan(f1))) ;
            else if  ( "acos".equals(op) )    res = (Object)(new Double(Math.acos(f1))) ;
            else if  ( "asin".equals(op) )    res = (Object)(new Double(Math.asin(f1))) ;
            else if  ( "atan".equals(op) )    res = (Object)(new Double(Math.atan(f1))) ;
            else if  ( "sqrt".equals(op) )    res = (Object)(new Double(Math.sqrt(f1))) ;
            else if  ( "log".equals(op) )  	  res = (Object)(new Double(Math.log(f1))) ;
            else if  ( "exp".equals(op) )  	  res = (Object)(new Double(Math.exp(f1))) ;
            else if  ( "floor".equals(op) )   res = (Object)(new Double(Math.floor(f1))) ;
            else if  ( "ceil".equals(op) )    res = (Object)(new Double(Math.ceil(f1))) ;
            else if  ( "abs".equals(op) )  	  res = (Object)(new Double(Math.abs(f1))) ;
        }

        return res;
    }

    private static Object evaluateExpression2(Operator o, Object o1, Object o2)
    {

        String op 	= o.getOperator();
        Object res 	= null;

        boolean o1isMatrix;
        boolean o2isMatrix;

        try {
          Matrix m1 = (Matrix)o1;
          o1isMatrix = true;
        } catch (Exception e) {
          o1isMatrix = false;
        }

        try {
          Matrix m2 = (Matrix)o2;
          o2isMatrix = true;
        } catch (Exception e) {
          o2isMatrix = false;
        }

        if (o1isMatrix) {
          Matrix f1 = (Matrix)o1;
          if (o2isMatrix) {
            Matrix f2 = (Matrix)o2;

            if       ( "+".equals(op) ) 	    res = (Object)(f1.plus(f2)) ;
            else if  ( "-".equals(op) ) 	    res = (Object)(f1.minus(f2)) ;
            else if  ( "*".equals(op) ) 	    res = (Object)(f1.times(f2)) ;
            else if  ( "/".equals(op) )  	    res = (Object)(f1.divide(f2)) ;

          } else {
            double f2 = ((Double)o2).doubleValue();

            if       ( "*".equals(op) ) 	    res = (Object)(f1.times(f2)) ;
            else if  ( "/".equals(op) )  	    res = (Object)(f1.divide(f2)) ;

            else if  ( "sort".equals(op) )  	res = (Object)(f1.sort(Math.round((float)f2))) ;
            else if  ( "find".equals(op) )  	res = (Object)(f1.find(f2)) ;
          }

        } else {
          double f1 = ((Double)o1).doubleValue();
          if (o2isMatrix) {
            Matrix f2 = (Matrix)o2;

            if       ( "*".equals(op) ) 	    res = (Object)(f2.times(f1)) ;

          } else {
            double f2 = ((Double)o2).doubleValue();

            if     	 ( "+".equals(op) ) 	    res = (Object)(new Double(f1 + f2)) ;
            else if  ( "-".equals(op) ) 	    res = (Object)(new Double(f1 - f2)) ;
            else if  ( "*".equals(op) ) 	    res = (Object)(new Double(f1 * f2)) ;
            else if  ( "/".equals(op) )  	    res = (Object)(new Double(f1 / f2)) ;

            else if  ( "^".equals(op) )  	    res = (Object)(new Double(Math.pow(f1,f2))) ;
            else if  ( "%".equals(op) )  	    res = (Object)(new Double(f1 % f2)) ;
            else if  ( "min".equals(op) )  	  res = (Object)(new Double(Math.min(f1,f2))) ;
            else if  ( "max".equals(op) )  	  res = (Object)(new Double(Math.max(f1,f2))) ;
            else if  ( "&".equals(op) )  	    res = (Object)(new Double(((f1>0)&&(f2>0)) ? (1) : (0))) ;
            else if  ( "|".equals(op) )  	    res = (Object)(new Double(((f1>0)||(f2>0)) ? (1) : (0))) ;
            else if  ( "<".equals(op) )  	    res = (Object)(new Double((f1<f2) ? (1) : (0))) ;
            else if  ( ">".equals(op) )  	    res = (Object)(new Double((f1>f2) ? (1) : (0))) ;
            else if  ( "=".equals(op) )  	  res = (Object)(new Double((f1==f2) ? (1) : (0))) ;

          }
        }

        return res;
    }

    private void initializeOperators()
    {
        operators     = new Operator[39];

        operators[0]  = new Operator("+"	, 2, 0);
        operators[1]  = new Operator("-"	, 2, 0);
        operators[2]  = new Operator("*"	, 2, 10);
        operators[3]  = new Operator("/"	, 2, 10);
        operators[4]  = new Operator("t"	, 1, 20);
        operators[5]  = new Operator("inv"	, 1, 20);
        operators[6]  = new Operator("diag"	, 1, 20);
        operators[7]  = new Operator("det"	, 1, 20);
        operators[8]  = new Operator("trace"	, 1, 20);
        operators[9]  = new Operator("rank"	, 1, 20);
        operators[10]  = new Operator("sum"	, 1, 20);
        operators[11]  = new Operator("prod"	, 1, 20);
        operators[12]  = new Operator("min"	, 1, 20);
        operators[13]  = new Operator("max"	, 1, 20);
        operators[14]  = new Operator("mean"	, 1, 20);
        operators[15]  = new Operator("cov"	, 1, 20);
        operators[16]  = new Operator("var"	, 1, 20);
        operators[17]  = new Operator("cor"	, 1, 20);
        operators[18]  = new Operator("sort"	, 2, 20);
        operators[19]  = new Operator("find"	, 2, 20);
        operators[20]  = new Operator("%"	, 2, 10);
        operators[21]  = new Operator("^"	, 2, 10);
        operators[22]  = new Operator("cos"	, 1, 20);
        operators[23]  = new Operator("sin"	, 1, 20);
        operators[24]  = new Operator("tan"	, 1, 20);
        operators[25]  = new Operator("acos"	, 1, 20);
        operators[26]  = new Operator("asin"	, 1, 20);
        operators[27]  = new Operator("atan"	, 1, 20);
        operators[28]  = new Operator("sqrt"	, 1, 20);
        operators[29]  = new Operator("exp"	, 1, 20);
        operators[30]  = new Operator("log"	, 1, 20);
        operators[31]  = new Operator("floor"	, 1, 20);
        operators[32]  = new Operator("ceil"	, 1, 20);
        operators[33]  = new Operator("abs"	, 1, 20);
        operators[34]  = new Operator("&"	, 2, 0);
        operators[35]  = new Operator("|"	, 2, 0);
        operators[36]  = new Operator("<"	, 2, 0);
        operators[37]  = new Operator(">"	, 2, 0);
        operators[38]  = new Operator("="	, 2, 0);

    }

    /***
     * gets the variable's value that was assigned previously
     */
    public Object getVariable(String s)
    {
    	return variables.get(s);
    }

    private Object getObject(String s)
    {
    	if ( s == null ) return null;

        Object res = null;
        try
        {
            res = new Double(Double.parseDouble(s));
        }
        catch(Exception e)
        {
        	return getVariable(s);
        }

        return res;
    }

    protected Operator[] getOperators()
    {
    	return operators;
    }

    protected class Operator
    {
    	private String op;
        private int type;
        private int priority;

        public Operator(String o, int t, int p)
        {
        	op = o;
            type = t;
            priority = p;
    	}

        public String getOperator()
        {
        	return op;
        }

        public void setOperator(String o)
        {
        	op = o;
        }

        public int getType()
        {
			return type;
        }

        public int getPriority()
        {
			return priority;
        }
    }

    protected class Node
    {
        public String 	nString		= null;
    	public Operator nOperator 	= null;
        public Node 	nLeft		= null;
        public Node 	nRight		= null;
        public Node 	nParent		= null;
        public int		nLevel		= 0;
        public Object  	nValue		= null;

    	public Node(String s) throws Exception
        {
        	init(null, s, 0);
        }

    	public Node(Node parent, String s, int level) throws Exception
        {
        	init(parent, s, level);
        }

        private void init(Node parent, String s, int level) throws Exception
        {
            s = removeIllegalCharacters(s);
            s = removeBrackets(s);
            s = addZero(s);
        	if ( checkBrackets(s) != 0 ) throw new Exception("Wrong number of brackets in [" + s + "]");

            nParent				= parent;
			nString 			= s;
            nValue				= getObject(s);
            nLevel 				= level;
            int sLength  		= s.length();
            int inBrackets		= 0;
            int startOperator   = 0;

            for (int i=0; i<sLength; i++)
            {
        		if ( s.charAt(i) == '(' )
                	inBrackets++;
                else if ( s.charAt(i) == ')' )
                	inBrackets--;
				else
                {
                    // the expression must be at "root" level
                    if ( inBrackets == 0 )
                    {
                        Operator o = getOperator(nString,i);
                        if ( o != null )
                        {
                        	// if first operator or lower priority operator
							if ( nOperator == null || nOperator.getPriority() >= o.getPriority() )
                            {
                            	nOperator 		= o;
                            	startOperator 	= i;
                            }
                        }
                    }
                }
            }

            if ( nOperator != null )
            {
                // one operand, should always be at the beginning
                if ( startOperator==0 && nOperator.getType() == 1 )
                {
                    // the brackets must be ok
                    if ( checkBrackets( s.substring( nOperator.getOperator().length() ) ) == 0 )
                    {
                        nLeft  = new Node( this, s.substring( nOperator.getOperator().length() ) , nLevel + 1);
                        nRight = null;
                        return;
                    }
                    else
                    	throw new Exception("Error during parsing... missing brackets in [" + s + "]");
                }
                // two operands
                else if ( startOperator > 0 && nOperator.getType() == 2 )
                {
//                    nOperator = nOperator;
                    nLeft 	= new Node( this, s.substring(0,  startOperator), nLevel + 1 );
                    nRight 	= new Node( this, s.substring(startOperator + nOperator.getOperator().length()), nLevel + 1);
                }
            }
        }

        private Operator getOperator(String s, int start)
        {
        	Operator[] operators = getOperators();
        	String temp = s.substring(start);
            temp = getNextWord(temp);

            for (int i=0; i<operators.length; i++)
            {
                if ( temp.startsWith(operators[i].getOperator()) ) {
                    return operators[i];}

            }
            return null;
        }

        private String getNextWord(String s)
        {
        	int sLength = s.length();
            for (int i=1; i<sLength; i++)
            {
            	char c = s.charAt(i);
                if ( (c > 'z' || c < 'a') && (c > '9' || c < '0') )
                	return s.substring(0, i);
            }
            return s;
    	}

        /***
         * checks if there is any missing brackets
         * @return true if s is valid
         */
        protected int checkBrackets(String s)
        {
            int sLength  	= s.length();
            int inBracket   = 0;

            for (int i=0; i<sLength; i++)
            {
                if      ( s.charAt(i) == '(' && inBracket >= 0 )
                    inBracket++;
                else if ( s.charAt(i) == ')' )
                    inBracket--;
            }

            return inBracket;
        }

        /***
         * returns a string that doesnt start with a + or a -
         */
        protected String addZero(String s)
        {
        	if ( s.startsWith("+") || s.startsWith("-") )
            {
            	int sLength  	= s.length();
                for (int i=0; i<sLength; i++)
                {
                	if ( getOperator(s, i) != null )
                    	return "0" + s;
                }
            }

            return s;
        }

        /***
         * displays the tree of the expression
         */
		public void trace()
        {
            String op = getOperator() == null ? " " : getOperator().getOperator() ;
            _D( op + " : " + getString() );
			if ( this.hasChild() )
            {
                if ( hasLeft() )
                	getLeft().trace();
                if ( hasRight() )
                	getRight().trace();
            }
        }

        protected boolean hasChild()
        {
        	return ( nLeft != null || nRight != null );
        }

        protected boolean hasOperator()
        {
        	return ( nOperator != null );
        }

        protected boolean hasLeft()
        {
        	return ( nLeft != null );
        }

        protected Node getLeft()
        {
        	return nLeft;
        }

        protected boolean hasRight()
        {
        	return ( nRight != null );
        }

        protected Node getRight()
        {
        	return nRight;
        }

        protected Operator getOperator()
        {
        	return nOperator;
        }

        protected int getLevel()
        {
        	return nLevel;
        }

        protected Object getValue()
        {
        	return nValue;
        }

        protected void setValue(Object f)
        {
        	nValue = f;
        }

        protected String getString()
        {
        	return nString;
        }

        /***
         * Removes spaces, tabs and brackets at the begining
         */
        public String removeBrackets(String s)
        {
            String res = s;
            if ( s.length() > 2 && res.startsWith("(") && res.endsWith(")") && checkBrackets(s.substring(1,s.length()-1)) == 0 )
            {
                res = res.substring(1, res.length()-1 );
            }
            if ( res != s )
                return removeBrackets(res);
			else
         	   return res;
        }

        /***
         * Removes illegal characters
         */
        public String removeIllegalCharacters(String s)
        {
        	char[] illegalCharacters = { ' ' };
			String res = s;

            for ( int j=0; j<illegalCharacters.length; j++)
            {
                int i = res.lastIndexOf(illegalCharacters[j], res.length());
                while ( i != -1 )
                {
                    String temp = res;
                    res = temp.substring(0,i);
                    res += temp.substring(i + 1);
                    i = res.lastIndexOf(illegalCharacters[j], s.length());
                }
            }
            return res;
        }

        protected void _D(String s)
        {
        	String nbSpaces = "";
            for (int i=0; i<nLevel; i++) nbSpaces += "  ";
            System.out.println(nbSpaces + "|" + s);
        }
    }

    protected static void _D(String s)
    {
    	System.err.println(s);
	}
}
