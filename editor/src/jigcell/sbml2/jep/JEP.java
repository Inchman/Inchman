/**
 * JEP - Java Expression Parser JEP is a Java package for parsing and evaluating mathematical expressions. It currently supports user defined
 * variables, constant, and functions. A number of common mathematical functions and constants are included. Author: Nathan Funk Copyright (C)
 * 2000 Nathan Funk JEP is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. JEP is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with
 * JEP; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

package jigcell.sbml2.jep;

import java.io.Reader;
import java.io.StringReader;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Vector;
import jigcell.sbml2.jep.function.PostfixMathCommand;
import jigcell.sbml2.jep.function.PostfixMathCommandI;
import jigcell.sbml2.jep.type.DoubleNumberFactory;
import jigcell.sbml2.jep.type.NumberFactory;

/**
 * The JEP class is the main interface with which the user should interact. It contains all neccessary methods to parse and evaluate
 * expressions.
 *
 * <p>
 * The most important methods are parseExpression(String), for parsing the mathematical expression, and getValue() for obtaining the value of
 * the expression.
 * </p>
 *
 * <p>
 * Visit <a href="http://jep.sourceforge.net/">http://jep.sourceforge.net/</a> for the newest version of JEP, and complete documentation. This
 * is an edited and enhanced version not aligned with the sourceforge project.
 * </p>
 *
 * @author Marc Vass
 */

public class JEP {

   /**
    * Debug flag for extra command line output
    */

   private final static boolean debug = false;

   /**
    * MathML functions
    */

   private final static HashMap mathMLFunctions;

   /**
    * Implicit multiplication option
    */

   protected boolean implicitMul;

   /**
    * Function Table
    */

   protected FunctionTable funTab;

   /**
    * Symbol Table
    */

   protected Hashtable symTab;

   /**
    * Error List
    */

   protected Vector errorList;

   /**
    * Traverse option
    */

   private boolean traverse;

   /**
    * Evaluator
    */

   private EvaluatorVisitor ev;

   /**
    * Node at the top of the parse tree
    */

   private Node topNode;

   /**
    * Number factory
    */

   private NumberFactory numberFactory;

   /**
    * The parser object
    */

   private Parser parser;

   static {
      mathMLFunctions = new HashMap ();
      mathMLFunctions.put ("sin", new PostfixMathCommand (1));
      mathMLFunctions.put ("cos", new PostfixMathCommand (1));
      mathMLFunctions.put ("tan", new PostfixMathCommand (1));
      mathMLFunctions.put ("sec", new PostfixMathCommand (1));
      mathMLFunctions.put ("csc", new PostfixMathCommand (1));
      mathMLFunctions.put ("cot", new PostfixMathCommand (1));
      mathMLFunctions.put ("sech", new PostfixMathCommand (1));
      mathMLFunctions.put ("csch", new PostfixMathCommand (1));
      mathMLFunctions.put ("coth", new PostfixMathCommand (1));
      mathMLFunctions.put ("arcsech", new PostfixMathCommand (1));
      mathMLFunctions.put ("arcsec", new PostfixMathCommand (1));
      mathMLFunctions.put ("arccsc", new PostfixMathCommand (1));
      mathMLFunctions.put ("arccsch", new PostfixMathCommand (1));
      mathMLFunctions.put ("arccot", new PostfixMathCommand (1));
      mathMLFunctions.put ("arccoth", new PostfixMathCommand (1));
      mathMLFunctions.put ("arcsin", new PostfixMathCommand (1));
      mathMLFunctions.put ("arccos", new PostfixMathCommand (1));
      mathMLFunctions.put ("arctan", new PostfixMathCommand (1));
      mathMLFunctions.put ("sinh", new PostfixMathCommand (1));
      mathMLFunctions.put ("cosh", new PostfixMathCommand (1));
      mathMLFunctions.put ("tanh", new PostfixMathCommand (1));
      mathMLFunctions.put ("arcsinh", new PostfixMathCommand (1));
      mathMLFunctions.put ("arccosh", new PostfixMathCommand (1));
      mathMLFunctions.put ("arctanh", new PostfixMathCommand (1));
      mathMLFunctions.put ("log", new PostfixMathCommand (2));
      mathMLFunctions.put ("ln", new PostfixMathCommand (1));
      mathMLFunctions.put ("floor", new PostfixMathCommand (1));
      mathMLFunctions.put ("ceiling", new PostfixMathCommand (1));
      mathMLFunctions.put ("factorial", new PostfixMathCommand (1));
      mathMLFunctions.put ("exp", new PostfixMathCommand (1));
      mathMLFunctions.put ("root", new PostfixMathCommand (2));
      mathMLFunctions.put ("abs", new PostfixMathCommand (1));
      mathMLFunctions.put ("plus", new PostfixMathCommand (2));
      mathMLFunctions.put ("minus", new PostfixMathCommand (2));
      mathMLFunctions.put ("times", new PostfixMathCommand (2));
      mathMLFunctions.put ("divide", new PostfixMathCommand (2));
      mathMLFunctions.put ("eq", new PostfixMathCommand (2));
      mathMLFunctions.put ("neq", new PostfixMathCommand (2));
      mathMLFunctions.put ("gt", new PostfixMathCommand (2));
      mathMLFunctions.put ("lt", new PostfixMathCommand (2));
      mathMLFunctions.put ("geq", new PostfixMathCommand (2));
      mathMLFunctions.put ("leq", new PostfixMathCommand (2));
      mathMLFunctions.put ("and", new PostfixMathCommand (2));
      mathMLFunctions.put ("or", new PostfixMathCommand (2));
      mathMLFunctions.put ("xor", new PostfixMathCommand (2));
      mathMLFunctions.put ("not", new PostfixMathCommand (1));
      mathMLFunctions.put ("power", new PostfixMathCommand (1));
   }

   public static boolean isMathMLFunction (String s) {
      return mathMLFunctions.containsKey (s);
   }

   /**
    * Creates a new JEP instance with the default settings.
    *
    * <p>
    * Traverse = false<br> Implicit multiplication = false<br> Number Factory = DoubleNumberFactory
    * </p>
    */

   public JEP () {
      topNode = null;
      traverse = false;
      implicitMul = false;
      numberFactory = new DoubleNumberFactory ();
      initFunTab ();
      symTab = new Hashtable ();
      errorList = new Vector ();
      ev = new EvaluatorVisitor ();
      parser = new Parser (new StringReader (""));

      //Ensure errors are reported for the initial expression
      //e.g. No expression entered
      parseExpression ("");
   }

   /**
    * Creates a new JEP instance with custom settings. If the numberFactory_in is null, the default number factory is used.
    *
    * @param traverse_in The traverse option.
    * @param implicitMul_in The implicit multiplication option.
    * @param numberFactory_in The number factory to be used.
    */

   public JEP (boolean traverse_in, boolean allowUndeclared_in, boolean implicitMul_in, NumberFactory numberFactory_in) {
      topNode = null;
      traverse = traverse_in;
      implicitMul = implicitMul_in;
      if (numberFactory_in == null)
         numberFactory = new DoubleNumberFactory ();
      else
         numberFactory = numberFactory_in;
      initFunTab ();
      symTab = new Hashtable ();
      errorList = new Vector ();
      ev = new EvaluatorVisitor ();
      parser = new Parser (new StringReader (""));

      //Ensure errors are reported for the initial expression
      //e.g. No expression entered
      parseExpression ("");
   }

   /**
    * Adds a new function to the parser. This must be done before parsing an expression so the parser is aware that the new function may be
    * contained in the expression.
    *
    * @param functionName The name of the function
    * @param function The function object that is used for evaluating the function
    */

   public void addFunction (String functionName, PostfixMathCommandI function) {
      funTab.put (functionName, function);
   }

   /**
    * Adds the standard functions to the parser. If this function is not called before parsing an expression, functions such as sin() or cos()
    * would produce an "Unrecognized function..." error. In most cases, this method should be called immediately after the JEP object is
    * created.
    */

   public void addStandardFunctions () {

      //add functions to Function Table
      funTab.put ("sin", new PostfixMathCommand (1));
      funTab.put ("cos", new PostfixMathCommand (1));
      funTab.put ("tan", new PostfixMathCommand (1));
      funTab.put ("sec", new PostfixMathCommand (1));
      funTab.put ("csc", new PostfixMathCommand (1));
      funTab.put ("cot", new PostfixMathCommand (1));
      funTab.put ("sech", new PostfixMathCommand (1));
      funTab.put ("csch", new PostfixMathCommand (1));
      funTab.put ("coth", new PostfixMathCommand (1));
      funTab.put ("arcsech", new PostfixMathCommand (1));
      funTab.put ("arcsec", new PostfixMathCommand (1));
      funTab.put ("arccsc", new PostfixMathCommand (1));
      funTab.put ("arccsch", new PostfixMathCommand (1));
      funTab.put ("arccot", new PostfixMathCommand (1));
      funTab.put ("arccoth", new PostfixMathCommand (1));
      funTab.put ("arcsin", new PostfixMathCommand (1));
      funTab.put ("arccos", new PostfixMathCommand (1));
      funTab.put ("arctan", new PostfixMathCommand (1));
      funTab.put ("sinh", new PostfixMathCommand (1));
      funTab.put ("cosh", new PostfixMathCommand (1));
      funTab.put ("tanh", new PostfixMathCommand (1));
      funTab.put ("arcsinh", new PostfixMathCommand (1));
      funTab.put ("arccosh", new PostfixMathCommand (1));
      funTab.put ("arctanh", new PostfixMathCommand (1));
      funTab.put ("log", new PostfixMathCommand (2));
      funTab.put ("log10", new PostfixMathCommand (1));
      funTab.put ("ln", new PostfixMathCommand (1));
      funTab.put ("floor", new PostfixMathCommand (1));
      funTab.put ("ceiling", new PostfixMathCommand (1));
      funTab.put ("factorial", new PostfixMathCommand (1));
      funTab.put ("exp", new PostfixMathCommand (1));
      funTab.put ("sqrt", new PostfixMathCommand (1));
      funTab.put ("root", new PostfixMathCommand (2));
      funTab.put ("abs", new PostfixMathCommand (1));
      funTab.put ("delay", new PostfixMathCommand (2));
      funTab.put ("sqrt", new PostfixMathCommand (1));
   }

   /**
    * Reports information on the errors that occured during the most recent action.
    *
    * @return A string containing information on the errors, each separated by a newline character; null if no error has occured
    */

   public String getErrorInfo () {
      if (hasError ()) {
         String str = "";

         // iterate through all errors and add them to the return string
         for (int i = 0; i < errorList.size (); i++)
            str += errorList.elementAt (i) + "\n";

         return str;
      } else
         return null;
   }

   /**
    * Returns the number factory.
    */

   public NumberFactory getNumberFactory () {
      return numberFactory;
   }

   // get the mathml2.0/sbml2 string for expression_in
   public String getSBMLString (String expression_in) throws ParseException {
      Reader reader = new StringReader (expression_in);
      symTab = new Hashtable (); // reset the table of variables
      try {

         // try parsing
         errorList.removeAllElements ();
         topNode = parser.parseStream (reader, this);
      } catch (Throwable e) {

         // an exception was thrown, so there is no parse tree
         topNode = null;

         // check the type of error
         if (e instanceof ParseException)

            // the ParseException object contains additional error
            // information
            errorList.addElement (((ParseException) e).getErrorInfo ());
         else {

            // if the exception was not a ParseException, it was most
            // likely a syntax error
            if (debug) {
               System.out.println (e.getMessage ());
               e.printStackTrace ();
            }
            errorList.addElement ("Syntax error");
         }
      }

      // If traversing is enabled, print a dump of the tree to
      // standard output
      if (!hasError ()) {
         ParserVisitor v = new ParserDumpVisitor ();
         topNode.jjtAccept (v, null);
         return (((ParserDumpVisitor) v).getSBMLString ());
      } else
         throw new ParseException(getErrorInfo());
   }

   /**
    * Returns the top node of the expression tree. Because all nodes are pointed to either directly or indirectly, the entire expression tree
    * can be accessed through this node. It may be used to manipulate the expression, and subsequently evaluate it manually.
    *
    * @return The top node of the expression tree
    */

   public Node getTopNode () {
      return topNode;
   }

   // get the variables in the symbol table
   public Enumeration getVariables () {
      return symTab.elements ();
   }

   /**
    * Gets the variables in the given equation.
    *
    * @param eq The equation.
    *
    * @return An Enumeration of the variables.
    */

   public Enumeration getVariablesInExpression (String eq) {
      parseExpression (eq);
      return getVariables ();
   }

   /**
    * Returns true if an error occured during the most recent action (parsing or evaluation).
    *
    * @return Returns <code>true</code> if an error occured during the most recent action (parsing or evaluation).
    */

   public boolean hasError () {
      return !errorList.isEmpty ();
   }

   /**
    * Creates a new FunctionTable object as funTab.
    */

   public void initFunTab () {

      //Init FunctionTable
      funTab = new FunctionTable ();
   }

   // is the string a function
   public boolean isFunction (String s) {
      return funTab.containsKey (s);
   }

   /**
    * Parses the expression. If there are errors in the expression, they are added to the <code>errorList</code> member.
    *
    * @param expression_in The input expression string
    */

   public void parseExpression (String expression_in) {
      Reader reader = new StringReader (expression_in);
      symTab = new Hashtable (); // reset the table of variables
      try {

         // try parsing
         errorList.removeAllElements ();
         topNode = parser.parseStream (reader, this);
      } catch (Throwable e) {

         // an exception was thrown, so there is no parse tree
         topNode = null;

         // check the type of error
         if (e instanceof ParseException)

            // the ParseException object contains additional error
            // information
            errorList.addElement (((ParseException) e).getErrorInfo ());
         else {

            // if the exception was not a ParseException, it was most
            // likely a syntax error
            if (debug) {
               System.out.println (e.getMessage ());
               e.printStackTrace ();
            }
            errorList.addElement ("Syntax error");
         }
      }

      // If traversing is enabled, print a dump of the tree to
      // standard output
      if (traverse && !hasError ()) {
         ParserVisitor v = new ParserDumpVisitor ();
         topNode.jjtAccept (v, null);
         System.out.println (((ParserDumpVisitor) v).getSBMLString ());
      }
   }

   /**
    * Removes a function from the parser.
    *
    * @return If the function was added earlier, the function class instance is returned. If the function was not present, <code>null</code> is
    *         returned.
    */

   public Object removeFunction (String name) {
      return funTab.remove (name);
   }

   /**
    * Sets the value of the implicit multiplication option. If this option is set to true before parsing, implicit multiplication will be
    * allowed. That means that an expression such as
    * <pre>"1 2"</pre>
    * is valid and is interpreted as
    * <pre>"1*2"</pre>
    * .
    *
    * <p>
    * The default value is false.
    * </p>
    *
    * @param value The boolean implicit multiplication option.
    */

   public void setImplicitMul (boolean value) {
      implicitMul = value;
   }

   /**
    * Sets the value of the traverse option. setTraverse is useful for debugging purposes. When traverse is set to true, the parse-tree will be
    * dumped to the standard ouput device.
    *
    * <p>
    * The default value is false.
    * </p>
    *
    * @param value The boolean traversal option.
    */

   public void setTraverse (boolean value) {
      traverse = value;
   }
}
