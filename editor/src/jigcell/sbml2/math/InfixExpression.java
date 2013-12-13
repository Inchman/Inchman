package jigcell.sbml2.math;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

/**
 * This code takes a parse tree composed of Nodes and handles all things related
 * to generating an infix expression.
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See
 * LICENSE for more details.
 *
 * @author Jason Zwolak
 * @author Nicholas Allen
 */

public class InfixExpression {

    private boolean lookupSymbols = true;
    private ArrayList expression;
    private SymbolTable symbolTable;
    private SymbolTable localSymbols;
    private Node parseTree;
    private boolean useHtml;
                          // if true then the infix expression will be generated
                          // as an html string.  Errors found in the Nodes will
                          // be colored in red text.
    // map provides a map from SBML constants to Fortran constants (e.g., true)
    private static Hashtable map = null;
    // order provides integers for operator precedence.
    private static Hashtable order = null;
    // associative is true if this operator level is associative for all
    // operators at this level.
    private static Hashtable associative = null;

    /*
     * Public Constants
     *
     * for mode
     */
    public static final int DONT_LOOKUP_SYMBOLS = 1;
    public static final int NORMAL              = 0;

    public InfixExpression(
        SymbolTable ls,
        SymbolTable st,
        Node parseTree,
        int mode
    )
        throws Exception
    {
        this.useHtml = false;
        if ( map == null || order == null || associative == null )
            initMaps();
        this.symbolTable  = st;
        this.localSymbols = ls;
        this.parseTree    = parseTree;
        if ( mode == DONT_LOOKUP_SYMBOLS ) lookupSymbols = false;
        expression = null;
    }

    public InfixExpression(
        SymbolTable ls,
        SymbolTable st,
        Node parseTree
    )
        throws Exception
    {
        this(ls,st,parseTree,NORMAL);
    }

    public InfixExpression(
        SymbolTable ls,
        SymbolTable st,
        MathMLExpression mathML
    )
        throws Exception
    {
        this(ls,st,mathML.getExpression(),NORMAL);
    }

    public InfixExpression( SymbolTable st, MathMLExpression mathML)
        throws Exception
    {
        this(new SymbolTable(), st, mathML);
    }

    public InfixExpression( SymbolTable st, Node parseTree )
        throws Exception
    {
        this(new SymbolTable(), st, parseTree, NORMAL );
    }

    public InfixExpression( MathMLExpression mathML )
        throws Exception
    {
        this(new SymbolTable(), new SymbolTable(), mathML);
    }

    public InfixExpression( Node parseTree )
        throws Exception
    {
        this(new SymbolTable(), new SymbolTable(), parseTree, NORMAL );
    }

    public InfixExpression( Node parseTree, int mode )
        throws Exception
    {
        this(new SymbolTable(), new SymbolTable(), parseTree, mode );
    }

    static private void initMaps() {

        map = new Hashtable(100);
        map.put( "true",       "1>0"  );
        map.put( "false",      "0>1" );
        map.put( "eq",         "=="    );
        map.put( "neq",        "!="    );
        map.put( "gt",         ">"    );
        map.put( "lt",         "<"    );
        map.put( "geq",        ">="    );
        map.put( "leq",        "<="    );
        map.put( "divide",     "/"       );
        map.put( "power",      "^"      );
        map.put( "and",        "&"   );
        map.put( "or",         "|"    );
        map.put( "times",      "*"       );
        map.put( "plus",       "+"       );
        map.put( "minus",      "-"       );
        map.put( "not",        "NOT"   );
        map.put( "abs",        "abs"     );
        map.put( "exp",        "exp"     );
        map.put( "ln",         "log"     );
        map.put( "floor",      "flr"   );
        map.put( "ceiling",    "flr" );
        map.put( "sin",        "sin"     );
        map.put( "cos",        "cos"     );
        map.put( "tan",        "tan"     );
        map.put( "csc",        "1/sin"     );
        map.put( "sec",        "1/cos"     );
        map.put( "cot",        "1/tan"     );
        map.put( "csch",        "1/sinh"     );
        map.put( "sech",        "1/cosh"     );
        map.put( "coth",        "1/tanh"     );
        map.put( "sinh",       "sinh"    );
        map.put( "cosh",       "cosh"    );
        map.put( "tanh",       "tanh"    );
        map.put( "arcsin",     "asin"    );
        map.put( "arccos",     "acos"    );
        map.put( "arctan",     "atan"    );
        map.put( "arccsc",     "asin"    );
        map.put( "arcsec",     "acos"    );
        map.put( "arccot",     "atan"    );
        map.put( "arcsinh",     "log"    );
        map.put( "arccosh",     "log"    );
        map.put( "arctanh",     "0.5*log"    );
        map.put( "arccsch",     "log"    );
        map.put( "arcsech",     "log"    );
        map.put( "arccoth",     "0.5*log"    );
        map.put( "exponentiale", "exp(1)" );
        map.put( "pi",           "pi");

        order = new Hashtable(30);
        order.put( "power",  new Integer(1) );
        order.put( "times",  new Integer(2) );
        order.put( "divide", new Integer(2) );
        order.put( "plus",   new Integer(3) );
        order.put( "minus",  new Integer(3) );
        order.put( "eq",     new Integer(4) );
        order.put( "neq",    new Integer(4) );
        order.put( "gt",     new Integer(4) );
        order.put( "lt",     new Integer(4) );
        order.put( "geq",    new Integer(4) );
        order.put( "leq",    new Integer(4) );
        order.put( "not",    new Integer(5) );
        order.put( "and",    new Integer(6) );
        order.put( "or",     new Integer(7) );
        order.put( "xor",     new Integer(6) );

        associative = new Hashtable(12);
        associative.put( new Integer(1), new Boolean(true)  );
        associative.put( new Integer(2), new Boolean(false)  );
        associative.put( new Integer(3), new Boolean(false)  );
        associative.put( new Integer(4), new Boolean(false) );
        associative.put( new Integer(5), new Boolean(true)  );
        associative.put( new Integer(6), new Boolean(true)  );
        associative.put( new Integer(7), new Boolean(true)  );

    }

    /**
     * Changes the function map to map <i>function</i> to <i>targetFunction</i>.
     * For example: if you have the <i>ln</i> function its default map is to the
     * <i>log</i> function.  You can call this function to change its map to be
     * something else like the <i>ln</i> function.  The call looks like:
     * modifyMap("ln","ln").  To change it back call modifyMap("ln","log").
     * NOTE: map is static so the changes affect all InfixExpressions.
     */
    static public void modifyMap( String function, String targetFunction )
        throws Exception
    {
        if ( map == null ) initMaps();
        if ( !map.containsKey(function) ) {
            throw new Exception("function "+function+" does not exist in map.");
        }
        map.put(function,targetFunction);
    }

    /**
     * Convert a parse tree into Fortran 90 code.  The Fortran 90 code is
     * returned as a list of strings.  The list can then be merged and line
     * breaks added to ensure the Fortran 90 line length limit is not exceeded.
     * <p>
     * When identifiers are reached they are looked up in the symbol table.  If
     * a Fortran identifier is not found then an Exception is thrown.
     * <p>
     * case 1
     * cn, ci, true, false, exponentiale
     * These are terminal nodes and will return a literal, identifier, or
     * constant.
     * <p>
     * case 2
     * eq, neq, gt, lt, geq, leq, divide, power, and, or, times, plus
     * These are nodes with 2 operands.  Recursion will be performed on the
     * operands and the first will be placed left of the operator and the second
     * will be placed right of the operator.  Parenthesis will be added around a
     * child if it is an operator of lower precidence than the current one.
     * <p>
     * case 3
     * minus
     * These are nodes with 1 or 2 operands.  If there is only one then the
     * operand is placed to the right of the operator.  Otherwise this node is
     * treated like a binary node.
     * <p>
     * case 4
     * root, log, xor
     * Special functions that are builtin.  Sometimes these are implemented
     * inline because they have no Fortran equivalent.
     * <p>
     * case 5
     * abs, exp, ln, floor, ceiling, not, and all Fortran supported
     * trig functions
     * Generate the appropriate Fortran function and pass as an argument the
     * child.
     * <p>
     * case 6
     * function
     * A user defined function with arbitrary arguments.  Generate
     * "USERFUNC( ARGS )".
     * <p>
     * default case
     * csymbol, nan, inf, piecewise, piece, otherwise, diff, pi, factorial
     * These nodes are unrecognized, unhandlable, or not yet implemented and
     * will cause an Exception.
     */
    private ArrayList genExpression( Node current ) throws Exception {
        ArrayList retArray = new ArrayList();
        if ( current == null ) return retArray;
        String name = current.getSimpleName().toLowerCase();
        if ( current.hasError() && useHtml ) retArray.add("<span style=\"color:red;\">");
        if ( name.equals("cn") ) {
            retArray.add( current.getValue());

        } else if ( name.equals("ci") ) {
            String value = (String)current.getValue();
            String sym = lookupSymbols ?
                localSymbols.lookupSymbolById(value):
                value
            ;
            if ( sym == null && !lookupSymbols )
                throw new Exception("Symbol not looked up and null.");
            if ( sym == null ) { sym = symbolTable.lookupSymbolById(value); }
            if ( sym == null )
                throw new Exception("Symbol not found: '"+value+"'");
            retArray.add(sym);

        } else if ( name.matches("true|false|exponentiale|pi") ) {
            retArray.add(map.get(name));

        } else if (name.equals ("xor")) {
            ArrayList leftChild, rightChild;
            leftChild  = genExpression(current.getChild(0));
            rightChild = genExpression(current.getChild(1));
            int opPrec = ((Integer)order.get(name)).intValue();
            Integer temp = (Integer)order.get(
                current.getChild(0).getSimpleName().toLowerCase()
            );
            if ( temp != null && temp.intValue() > opPrec) {
                leftChild.add(0,"(");
                leftChild.add(")");
            }
            temp = (Integer)order.get(
                current.getChild(1).getSimpleName().toLowerCase()
            );
            if ( temp != null &&
                    ( temp.intValue() > opPrec ||
                      temp.intValue() == opPrec &&
                      !((Boolean)associative.get(temp)).booleanValue() )
            ) {
                rightChild.add(0,"(");
                rightChild.add(")");
            }
            retArray.add ("(");
            retArray.addAll(leftChild);
            retArray.add (") & (NOT (");
            retArray.addAll (rightChild);
            retArray.add(")) | (");
            retArray.addAll (rightChild);
            retArray.add (") & (NOT (");
            retArray.addAll(leftChild);
            retArray.add ("))");
        } else if ( name.matches(
            "eq|neq|gt|lt|geq|leq|divide|power|and|or|times|plus"
        ) || (
            name.equals("minus") && current.getNumChildren() == 2
        ) ) {
            /*
            - If either child is an operator of lower precedence then enclose
              that child in parenthesis (lower precedence means a higher number
              in the map "order").
            - Call genExpression on each child and assemble the string
              leftchild+operator+rightchild.
            - Note here that >= opPrec is used instead of > opPrec.  This is
              because some of these binary operators are not associative with
              other operators of the same precedence.
            */
            String op;
            ArrayList leftChild, rightChild;
            op = (String)map.get(name);
            leftChild  = genExpression(current.getChild(0));
            rightChild = genExpression(current.getChild(1));
            int opPrec = ((Integer)order.get(name)).intValue();
            Integer temp = (Integer)order.get(
                current.getChild(0).getSimpleName().toLowerCase()
            );
            if ( temp != null && temp.intValue() > opPrec) {
                leftChild.add(0,"(");
                leftChild.add(")");
            }
            temp = (Integer)order.get(
                current.getChild(1).getSimpleName().toLowerCase()
            );
            if ( temp != null &&
                    ( temp.intValue() > opPrec ||
                      temp.intValue() == opPrec &&
                      !((Boolean)associative.get(temp)).booleanValue() )
            ) {
                rightChild.add(0,"(");
                rightChild.add(")");
            }
            retArray.addAll(leftChild);
            retArray.add(" " + op + " ");
            retArray.addAll(rightChild);

        } else if ( name.equals("minus") ) {
            // Only the unary minus operator is treated here.  Detection of the
            // binary minus operator is above.
            int opPrec = ((Integer)order.get(name)).intValue();
            Integer temp = (Integer)order.get(
                current.getChild(0).getSimpleName().toLowerCase()
            );
            ArrayList child = genExpression(current.getChild(0));
            if ( temp != null && temp.intValue() >= opPrec ) {
                child.add(0,"(");
                child.add(")");
            }
            retArray.add(map.get(name));
            retArray.addAll(child);

        } else if ( name.equals("root") ) {
            // only works for square roots
            double d = Double.parseDouble(current.getChild(1).getValue());
            if ( d == 2.0 ) {
               retArray.add("sqrt(");
               retArray.addAll(genExpression(current.getChild(0)));
               retArray.add(")");
            } else {
               retArray.add("(");
               retArray.addAll(genExpression(current.getChild(0)));
               retArray.add(")^(1/" + d + ")");
            }

        } else if ( name.equals("ln") ){
             retArray.add("("+map.get("ln")+"(");
             retArray.addAll(genExpression(current.getChild(0)));
             retArray.add("))");
        }
        else if ( name.equals("log") ) {
            // special code is constructed for logs of bases other than exp and
            // 10
            int n = Integer.parseInt(current.getChild(0).getValue());
            if ( n == 10 ) {
                retArray.add("log10(");
                retArray.addAll(genExpression(current.getChild(1)));
                retArray.add(")");
            } else {
                retArray.add("(log(");
                retArray.addAll(genExpression(current.getChild(1)));
                retArray.add(")/log(");
                retArray.add(String.valueOf(n));
                retArray.add("))");
            }

        } else if ( name.equals("not") ) {
            // not is treated different than abs, exp, ln, etc. because it isn't
            // a function.
            String op = (String)map.get(name);
            int opPrec = ((Integer)order.get(name)).intValue();
            Integer temp = (Integer)order.get(
                current.getChild(0).getSimpleName().toLowerCase()
            );
            ArrayList child = genExpression(current.getChild(0));
            if ( temp != null && temp.intValue() > opPrec ) {
                child.add(0,"(");
                child.add(")");
            }
            retArray.add(op);
            retArray.addAll(child);

        } else if ( name.matches(
            "abs|exp|floor|sin|cos|tan|sinh|cosh|tanh|csc|sec|cot|csch|sech|coth|arcsin|arccos|arctan"
        ) ) {
            retArray.add(map.get(name));
            retArray.add("(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add(")");
        } else if (name.equals ("ceiling")) {
            retArray.add(map.get(name));
            retArray.add("(0.9999999+");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add(")");
        } else if ( name.matches(
            "arccsc|arcsec|arccot"
        ) ) {
            retArray.add(map.get(name));
            retArray.add("(1/(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add("))");

        } else if (
            name.equals("function")
        ) {
            // look up user function name and
            // return name+"("+args+")";
            String value = (String)current.getValue();
            String sym = lookupSymbols ?
                localSymbols.lookupSymbolById(value):
                value
            ;
            if ( sym == null && !lookupSymbols )
                throw new Exception("Symbol not looked up and null.");
            sym = (sym==null)?symbolTable.lookupSymbolById(value):sym;
            if ( sym == null )
                throw new Exception("Symbol not found: "+value);
            ArrayList args = new ArrayList();
            int numChildren = current.getNumChildren();
            for ( int i = 0; i < numChildren; i++ ) {
                if ( i != 0 ) args.add(",");
                args.addAll(genExpression(current.getChild(i)));
            }
            retArray.add(sym);
            retArray.add("(");
            retArray.addAll(args);
            retArray.add(")");

        } else if(name.equals("piecewise")){
            if(current.getNumChildren() > 2)
                throw new Exception("Only if then else is supported");
            retArray.add("if (");
            retArray.addAll(genExpression(current.getChild(0).getChild(1)));
            retArray.add(") then (");
            retArray.addAll(genExpression(current.getChild(0).getChild(0)));
            retArray.add(") else (");
            retArray.addAll(genExpression(current.getChild(1).getChild(0)));
            retArray.add(")");
        } else if (name.equals ("arcsinh")) {
            retArray.add(map.get(name));
            retArray.add("((");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add(")+((");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add(")^2+1)^0.5)");
        } else if (name.equals ("arccosh")) {
            retArray.add(map.get(name));
            retArray.add("((");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add(")+((");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add(")^2-1)^0.5)");
        } else if (name.equals ("arctanh")) {
            retArray.add(map.get(name));
            retArray.add("((1+(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add("))/(1-(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add(")))");
        } else if (name.equals ("arcsech")) {
            retArray.add(map.get(name));
            retArray.add("((1+(1-(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add(")^2)^0.5)/(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add("))");
        } else if (name.equals ("arccsch")) {
            retArray.add(map.get(name));
            retArray.add("((1+(1+(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add(")^2)^0.5)/(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add("))");
        } else if (name.equals ("arccoth")) {
            retArray.add(map.get(name));
            retArray.add("((1+(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add("))/((");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add(")-1))");
        } else if (name.equals ("factorial")) {
            retArray.add ("((2*(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add (")+1/3)*pi)^0.5*(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add (")^(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add (")*exp(-(");
            retArray.addAll(genExpression(current.getChild(0)));
            retArray.add ("))");
        } else if (name.equals("csymbol")) {
           String symbol = current.getAttributes ().getValue ("definitionURL");
           if (symbol.equals ("http://www.sbml.org/sbml/symbols/time")) {
              retArray.add ("t");
           } else {
              throw new Exception ("Unsupported csymbol: " + symbol);
           }
        } else if (name.matches("nan|inf")) {
           throw new Exception(name+" not supported.");
        } else {
           throw new Exception("should never throw this, but did for "+name);
        }

        if ( current.hasError() && useHtml ) retArray.add("</span>");
        return retArray;

    }

    public Integer getOrder( Node node ) throws Exception {
        String simpleName = node.getSimpleName().toLowerCase();
        if ( !simpleName.equals("function") ) {
            return ((Integer)order.get(simpleName));
        }
        String value = (String)node.getValue();
        String sym = lookupSymbols ?
            localSymbols.lookupSymbolById(value):
            value
        ;
        if ( sym == null && !lookupSymbols )
            throw new Exception("Symbol not looked up and null.");
        sym = (sym==null)?symbolTable.lookupSymbolById(value):sym;
        if ( sym == null )
            throw new Exception("Symbol not found: "+value);
        return ((Integer)order.get(sym));
    }

    /**
     * Merges strings in expression and addes line breaks to maintain a width no
     * greater than width.  This method is specific to Fortran code as it adds
     * ampersands at the end of lines.  Indent will be added at the begining of
     * each new line (but not the first line).
     */
    public static String mergeExpression(
        ArrayList expression,
        String indent,
        int width,
        String prefix,
        boolean useHtml
    ) {
        StringBuffer curLine = new StringBuffer("");
        if (useHtml) curLine.append("<html>");
        if ( prefix != null ) curLine.append(prefix);
        for ( Iterator i = expression.iterator(); i.hasNext(); ) {
           String curToken = (String)i.next();
           curLine.append(curToken);
        }
        if (useHtml) curLine.append("</html>");
        return curLine.toString();
    }
    public static String mergeExpression(
        ArrayList expression,
        String indent,
        int width,
        String prefix
    ) {
        return mergeExpression( expression, indent, width, prefix, false );
    }

    public String getExpression()
        throws Exception
    {
        return getExpression("",80,"");
    }
    public String getExpression( String indent, int width )
        throws Exception
    {
        return getExpression( indent, width, "" );
    }
    public String getExpression( String indent, int width, String prefix )
        throws Exception
    {
        if ( expression == null ) expression = genExpression(parseTree);
        return mergeExpression( expression, indent, width, prefix, useHtml );
    }

    public void setUseHtml( boolean useHtml ) {
       this.useHtml = useHtml;
    }

}
