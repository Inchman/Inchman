package jigcell.sbml2.math;

import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.util.StringTokenizer;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import org.xml.sax.Attributes;
import org.xml.sax.helpers.AttributesImpl;
import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

/**
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 *
 * @author Jason Zwolak
 * @author Marc Vass
 * @author Nicholas Allen
 * @author Ranjit Randhawa
 */

public class MathMLExpression {

   /*
    * tabs - used when printing the parse tree.  tabs contains the current
    *        level of indentation.
    * mathMLParseTree - the MathML parse tree as constructed by MathMLHandler
    * lambdaFunction - is true if this mathml represents a lambda function as
    *                  defined in the MathML spec.  The
    *                  jigcell.modelbuilder.sbml parser returns the contents of
    *                  the lambda function only.  Therefore it is possible to
    *                  have a function with no arguments that is
    *                  indistinguishable from an expression.  The programmer
    *                  must therefore tell MathMLExpression if it is a lambda
    *                  function.
    * simpleTree - Represents a simple mathematical parse tree.  All
    *              annotations removed and all ambiguities removed.  The
    *              struction is a binary tree (except for some elements that
    *              are not handled: etc.)
    * functionArgs - The function argument names.
    * lookForFunctionArgs - Specifies if the genSimpleTree method should
    *                       look for function args or not.  See the
    *                       genSimpleTree method.
    */
   private static StringBuffer tabs = new StringBuffer ();
   private ArrayList functionArgs = new ArrayList ();
   private boolean lambdaFunction;
   private boolean lookForFunctionArgs = true;
   private Node mathMLParseTree;
   private Node simpleTree;
   private StringBuffer mathML;

   /**
    * Used as a return value in solve.  left is the left hand side of an equation and right is the right hand side.
    */

   protected class Equation {
      Node left;
      Node right;
   }

   /**
    * Receives events from SAX and builds a parse tree.  The parse tree exactly matches the MathML file.  For every element there is a node.
    * Each node has a child for each element it encloses.  The MathML <code>&lt;apply&gt;&lt;times/&gt;&lt;cn&gt;5&lt;cn&gt;&lt;cn&gt;2
    * &lt;cn&gt;&lt;/apply&gt;</code> has a root node with "apply" and 3 children with "times", "5", and "2" (the order is preserved).
    */

   protected class MathMLHandler extends DefaultHandler {
      protected Node root;
      private Node current;

      public void characters (char ch [], int start, int length) {
         current.appendToValue (new String (ch, start, length));
      }

      public void endDocument () {}

      public void endElement (String uri, String localName, String qName) {

         // Collapse white space according to the MathML spec.
         current.setValue (current.getValue ().trim ().replaceAll ("[ \t\n\r]+", " "));
         current = current.getParent ();
      }

      public void startDocument () {
         root = new Node ("", "", "TREEROOT", null);
         current = root;
      }

      public void startElement (String uri, String localName, String qName, Attributes attributes) {
         Node newNode = new Node (current, uri, localName, qName, new AttributesImpl (attributes));
         current.addChild (newNode);
         current = newNode;
      }
   }

   /**
    * Gets all identifiers in the parse tree cur and returns them.
    */

   public String [] getIdentifiers () throws Exception {
      return getIdentifiers (getExpression());
   }

   public static String [] getIdentifiers (Node cur) {
      return (String []) getIdentifiersAsList (cur).toArray (new String [0]);
   }

   public static void printParseTree (Node cur) {
      System.out.println (tabs + "node: " + cur.getQName () + ", value: " + cur.getValue ().trim ());
      tabs.append ("    ");
      for (int i = 0; i < cur.getNumChildren (); i++)
         printParseTree (cur.getChild (i));
      tabs.delete (tabs.length () - 4, tabs.length ());
   }

   public ArrayList getIdentifiersAsList () throws Exception {
      return getIdentifiersAsList (getExpression());
   }

   public static ArrayList getIdentifiersAsList (Node cur) {
      String type = cur.getSimpleName ().toLowerCase ();
      ArrayList arrayList = new ArrayList ();
      if (type.equals ("ci"))
         arrayList.add (cur.getValue ());
      else
         for (int i = 0; i < cur.getNumChildren (); i++)
            arrayList.addAll (getIdentifiersAsList (cur.getChild (i)));
      return arrayList;
   }

   public MathMLExpression (String mathML) throws SAXException, IOException, ParserConfigurationException {
      parse (new InputSource (new StringReader (mathML)));
   }

   public MathMLExpression (Node mathMLParseTree){
	      this.mathMLParseTree = mathMLParseTree;
   }
   
   /**
    * Converts the MathML parse tree into a traditional math parse tree and returns the result.  If this MathMLExpression is a function then the
    * function arguments are removed from the tree and just the function definition is returned.
    */

   public Node getExpression () throws Exception {
      if (simpleTree == null)
         simpleTree = genSimpleTree (mathMLParseTree);
      return simpleTree;
   }

   /**
    * Return an array of strings of the function arguments.  If this MathMLExpression is not a function then return null.
    */

   public String [] getFunctionArguments () throws Exception {
      if (simpleTree == null)
         simpleTree = genSimpleTree (mathMLParseTree);
      return (String []) functionArgs.toArray (new String [0]);
   }

   public Node getMathMLParseTree () {
      return mathMLParseTree;
   }

   public boolean isLambdaFunction () {
      return lambdaFunction;
   }

   /**
    * Tell this MathMLExpression if it represents a lambda function or not.  If it represents a lambda function it will look for function
    * arguments first then a mathematical expression.  If this MathMLExpression doesn't represent a lambda function then it will give an error
    * if it finds function arguments.  Note that a lambda function may have zero arguments (making it indistinguishable from an expression).
    */

   public void setLambdaFunction (boolean lambdaFunction) {
      this.lambdaFunction = lambdaFunction;
   }

   /**
    * Solves for var and returns a parse tree.  The parse tree contains the expression for the value of var assuming the expression represented
    * by this MathMLExpression is equal to zero.  This routine only works on
    * expressions containing only -, +, and *.  An exception is thrown if
    * other operations are found.  Furthermore, the tree must be a sum of products.  The returned tree is a simpleTree.
    */

   public Node solve (String var) throws Exception {
      if (simpleTree == null)
         simpleTree = genSimpleTree (mathMLParseTree);
      Equation equation = solve_ (simpleTree, var);
      if (equation.right == null)
         throw new Exception ("Variable doesn't occur on the right side "
               + "of its expression: " + var);
      if (equation.left == null)
         throw new Exception ("Variable occurs alone on the right side of its " 
               + "expression: " + var);

      // construct -left/right from equation left == - right * var
      Node minus = new Node ("minus", "");
      boolean isOne = false;
      try { 
         if ( Double.parseDouble(equation.right.getValue()) == 1.0 )
            isOne = true;
      } catch ( NumberFormatException e ) {};
      if ( isOne ) {
         minus.addChild(equation.left);
      } else {
         Node divide = new Node ("divide", "");
         minus.addChild (divide);
         divide.addChild (equation.left);
         divide.addChild (equation.right);
      }
      return minus;
   }

   /**
    * Create a simple mathematical parse tree (roughly a binary tree).  If simpleTree is already defined then return it.  Otherwise traverse
    * mathMLParseTree recursively.
    *
    * <p>
    * For all nodes be sure children are not null after a call to genSimpleTree if they were not null before a call.  Set lookForFunctionArgs to
    * false on most of the nodes (if they aren't a bvar or a skipable).  Be sure to make all children point to the current node as they by
    * default point to the parent of the child they were cloned from.
    * </p>
    *
    * <p>
    * case 1 For nodes of cn, ci, csymbol, true, false, nan, pi, inf, exp: return a copy of the node.  Exceptions: when there are children.
    * </p>
    *
    * <p>
    * case 2 For nodes of piecewise, piece, otherwise: return a copy of the node.  No exceptions.
    * </p>
    *
    * <p>
    * case 3 For node of apply: copy first child and operate as if the remaining nodes of apply are children of the first child.  Exception:
    * apply has no children.
    * </p>
    *
    * <p>
    * case 4 For nodes of eq, neq, gt, lt, geq, leq, divide, power, root, and, or, xor: recurse on each child, copy and use the resulting trees
    * from the recursion as the new children, return the copy.  Exceptions: children not equal to two.
    * </p>
    *
    * <p>
    * case 5 For nodes of times and plus: if no children return identity (0 for plus, 1 for times), if one child, recurse on child and return
    * child. If two or more children copy node and remove first child.  Recurse on copy and recurse on removed child.  Add the two recursions
    * as children of a fresh copy and return.  No exceptions.
    * </p>
    *
    * <p>
    * case 6 For nodes of minus: recurse on children and return copy with trees from children as new children.  Exceptions: more than 2
    * children, zero children.
    * </p>
    *
    * <p>
    * case 7 For node log: if two children then recurse on children and make copy with new children.  If one child then recurse on child and
    * make copy with &lt;cn&gt;10&lt;/cn&gt; as the first child and the other child from the recursion.  Exceptions: zero children, more than 2
    * children, if 2 children and first child is not logbase.
    * </p>
    *
    * <p>
    * case 8 For nodes logbase and degree: recurse on child and return child. Exception: not exactly one child.
    * </p>
    *
    * <p>
    * case 9 For node bvar: If lookForFunctionArgs and lambdaFunction then add to functionArgs.  Exceptions: not exactly one argument, argument
    * isn't ci, didn't expect bvar here. Note: the first non-skipable element that isn't a bvar sets lookForFunctionArgs = false.
    * </p>
    *
    * <p>
    * case 10 For nodes abs, exp, ln, floor, ceiling, factorial, not, and all trig functions: recurse on child and return copy with new child.
    * Exception: number of children isn't exactly one.
    * </p>
    *
    * <p>
    * case 11 (obsolete? the stripComments routine takes care of this) For nodes semantics and math: recurse on children and return one child.
    * If all children return null then return null. Exception: more than one child evaluated to non-null, no children. Note: these nodes are
    * skippable.
    * </p>
    *
    * <p>
    * case 12 (obsolete? the stripComments routine takes care of this) For node annotation: return null.  No Exceptions. Note: this is a
    * skippable node.
    * </p>
    *
    * <p>
    * case 13 (obsolete? the stripComments routine takes care of this) For node annotation-xml: If it is a content annotation recurse on its
   * child and return the new child.  If it is another type of annotation return null.  Exceptions: not exactly one child (for content only).
      * Note: this is skippable if it isn't content.
      * </p>
      *
      * <p>
      * case 14 For the root node call stripComments and then process the tree as normal. This case is also for lambda nodes because lambda used
      * to not be present in the MathML from jigcell.modelbuilder.sbml.  Lambda could be handled better by recognizing it as its own case and
      * collecting the bvars then building the expression.
      * </p>
      *
      * <p>
      * case 15 Sometimes apply takes a ci as its first argument.  In this case ci is a function name.  The ci is changed to function to be
      * distinguishable from other ci's.
      * </p>
      */

   private Node genSimpleTree (Node current) throws Exception {
      String name = current.getSimpleName ().toLowerCase ();
      if (name.equals ("ci") || name.equals ("csymbol") || name.equals ("true") || name.equals ("false") || name.equals ("notanumber") ||
            name.equals ("pi") || name.equals ("infinity") || name.equals ("exponentiale")) {
         if (current.getNumChildren () != 0)
            throw new Exception ("Expected 0 children for " + name + " and there are " + current.getNumChildren () + " children.");
         return new Node (current);
      } else if (name.equals ("cn")) {
         if (current.getNumChildren () == 0)
            return new Node (current);
         if (current.getNumChildren () != 1)
            throw new Exception ("Node cn must have 0 or 1 children.");
         Node newNode = new Node (current);
         String type = newNode.getAttributes ().getValue ("type");
         StringTokenizer values = new StringTokenizer (newNode.getValue ());
         if (type.equals ("e-notation"))
            newNode.setValue (String.valueOf (Double.parseDouble (values.nextToken () + "E" + values.nextToken ())));
         else if (type.equals ("rational"))
            newNode.setValue (String.valueOf (Double.parseDouble (values.nextToken ()) / Double.parseDouble (values.nextToken ())));
         else
            throw new Exception ("Unknown cn type: " + type + ".");
         newNode.removeAllChildren ();
         return newNode;
      } else if (name.equals ("piecewise") || name.equals ("piece") || name.equals ("otherwise")) {
         int numChildren = current.getNumChildren ();
         Node newNode = new Node (current);
         newNode.removeAllChildren ();
         if (name.equals ("piecewise")) {
            for (int i = 0; i < numChildren; i++)
               newNode.addChild (new Node (genSimpleTree (current.getChild (i))));
         } else if (name.equals ("piece")) {
            newNode.addChild (new Node (genSimpleTree (current.getChild (0))));
            newNode.addChild (new Node (genSimpleTree (current.getChild (1))));
         } else if (name.equals ("otherwise"))
            newNode.addChild (new Node (genSimpleTree (current.getChild (0))));
         return newNode;
      } else if (name.equals ("apply")) {
         int numChildren = current.getNumChildren ();
         if (numChildren == 0)
            throw new Exception ("Node apply must have at least one child.");
         Node newNode = new Node (current.getChild (0));
         // change the name from ci to function, note the value of this node
         // contains the function name.
         if (newNode.getSimpleName ().equalsIgnoreCase ("ci"))
            newNode.setQName ("function");
         newNode.children = new ArrayList (current.children);
         newNode.children.remove (0); // No longer need the function name as the
                                      // first child, the name is in the parent.
         return genSimpleTree (newNode);
      } else if (name.equals ("eq") || name.equals ("neq") || name.equals ("gt") || name.equals ("lt") || name.equals ("geq") ||
            name.equals ("leq") || name.equals ("divide") || name.equals ("power")) {
         if (current.getNumChildren () != 2)
            throw new Exception ("Number of children must be equal to 2 for binary " + "operator " + name + ".  There are " +
                  current.getNumChildren () + " children.");
         Node newNode = new Node (current);
         newNode.removeAllChildren ();
         newNode.addChild (genSimpleTree (current.getChild (0)));
         newNode.addChild (genSimpleTree (current.getChild (1)));
         return newNode;
      } else if (name.equals ("and") || name.equals ("or") || name.equals ("xor")) {
         if (current.getNumChildren () < 2)
            throw new Exception ("Operator " + name + " must have at least two children.");
         Node newNode = new Node (current);
         newNode.removeAllChildren ();
         for (int i = 0; i < current.getNumChildren (); i++)
            newNode.addChild (genSimpleTree (current.getChild (i)));
         while (newNode.getNumChildren () > 2) {
            Node rightNode = newNode.removeChild (newNode.getNumChildren () - 1);
            Node leftNode = newNode.removeChild (newNode.getNumChildren () - 1);
            Node middleNode = new Node ("math:" + name, "");
            middleNode.addChild (leftNode);
            middleNode.addChild (rightNode);
            newNode.addChild (middleNode);
         }
         return newNode;
      } else if (name.equals ("root")) {
         if (current.getNumChildren () == 1) {
            if (current.getChild (0).getSimpleName ().equalsIgnoreCase ("degree"))
               throw new Exception ("Operator root is missing its radicand.");
            Node degree = new Node ("math:degree", "");
            degree.addChild (new Node ("math:cn", "2"));
            current.addChild (degree);
         }
         if (current.getNumChildren () != 2)
            throw new Exception ("Number of children must be equal to 2 for binary " + "operator " + name + ".  There are " +
                  current.getNumChildren () + " children.");
         Node newNode = new Node (current);
         newNode.removeAllChildren ();
         if (current.getChild (0).getSimpleName ().equalsIgnoreCase ("degree")) {
            newNode.addChild (genSimpleTree (current.getChild (1)));
            newNode.addChild (genSimpleTree (current.getChild (0)));
         } else {
            newNode.addChild (genSimpleTree (current.getChild (0)));
            newNode.addChild (genSimpleTree (current.getChild (1)));
         }
         return newNode;

         // case 5
      } else if (name.equals ("times") || name.equals ("plus")) {
         int numChildren = current.getNumChildren ();
         switch (numChildren) {
            case 0 :if (name.equals ("times"))
                       return new Node ("cn", "1");
                    else
                       return new Node ("cn", "0");
            case 1 :return genSimpleTree (current.getChild (0));
            default :Node leftNode = genSimpleTree (current.getChild (0));
                     Node rightNode = new Node (current);
                     Node newNode = new Node (current);
                     rightNode.removeChild (0);
                     rightNode = genSimpleTree (rightNode);
                     newNode.removeAllChildren ();
                     newNode.addChild (leftNode);
                     newNode.addChild (rightNode);
                     return newNode;
         }

         // case 6
      } else if (name.equals ("minus")) {
         int numChildren = current.getNumChildren ();
         if (numChildren > 2 || numChildren < 1)
            throw new Exception ("Minus expects 1 or 2 children and has " + numChildren + ".");
         Node newNode = new Node (current);
         newNode.removeAllChildren ();
         newNode.addChild (genSimpleTree (current.getChild (0)));
         if (numChildren == 2)
            newNode.addChild (genSimpleTree (current.getChild (1)));
         return newNode;

         // case 7
      } else if (name.equals ("log")) {
         int numChildren = current.getNumChildren ();
         if (numChildren > 2 || numChildren < 1)
            throw new Exception ("Log expects 1 or 2 children and has " + numChildren + ".");
         Node newNode = new Node (current);
         newNode.removeAllChildren ();
         if (numChildren == 2) {
            if (!current.getChild (0).getSimpleName ().equalsIgnoreCase ("logbase"))
               throw new Exception ("Log expects first child to be logbase when there are 2 " + "children.  First child is " +
                     current.getChild (0).getSimpleName () + ".");
            newNode.addChild (genSimpleTree (current.getChild (0)));
            newNode.addChild (genSimpleTree (current.getChild (1)));
         } else {
            newNode.addChild (new Node ("cn", "10"));
            newNode.addChild (genSimpleTree (current.getChild (0)));
         }
         return newNode;

         // case 8
      } else if (name.equals ("logbase") || name.equals ("degree")) {
         if (current.getNumChildren () != 1)
            throw new Exception ("Logbase and degree expect exactly 1 child and have " + current.getNumChildren ());
         return genSimpleTree (current.getChild (0));

         // case 9
      } else if (name.equals ("bvar")) {
         if (!lambdaFunction)
            throw new Exception ("Bvar found in non-lambdaFunction.");
         if (!lookForFunctionArgs)
            throw new Exception ("Bvar found out of place in a lambdaFunction.");
         if (current.getNumChildren () != 1)
            throw new Exception ("Bvar expects exactly 1 child and has " + current.getNumChildren ());
         Node ciNode = current.getChild (0);
         if (!ciNode.getSimpleName ().equalsIgnoreCase ("ci"))
            throw new Exception ("Bvar expects its child to be ci and its child is " + ciNode.getSimpleName () + ".");
         functionArgs.add (ciNode.getValue ());
         return null;

         // case 10
      } else if (name.equals ("abs") || name.equals ("exp") || name.equals ("ln") || name.equals ("floor") || name.equals ("ceiling") ||
            name.equals ("factorial") || name.equals ("not") || name.equals ("sin") || name.equals ("cos") || name.equals ("tan") ||
            name.equals ("sec") || name.equals ("csc") || name.equals ("cot") || name.equals ("sinh") || name.equals ("cosh") ||
            name.equals ("tanh") || name.equals ("sech") || name.equals ("csch") || name.equals ("coth") || name.equals ("arcsin") ||
            name.equals ("arccos") || name.equals ("arctan") || name.equals ("arcsec") || name.equals ("arccsc") || name.equals ("arccot") ||
            name.equals ("arcsinh") || name.equals ("arccosh") || name.equals ("arctanh") || name.equals ("arcsech") || name.equals ("arccsch") ||
            name.equals ("arccoth")) {
         if (current.getNumChildren () != 1)
            throw new Exception ("The function " + name + " epects one argument and has " + current.getNumChildren () + ".");
         Node newNode = new Node (current);
         newNode.removeAllChildren ();
         newNode.addChild (genSimpleTree (current.getChild (0)));
         return newNode;

         // case 11 - obsolete (?) by stripComments
      } else if (name.equals ("semantics") || name.equals ("math")) {
         if (current.getNumChildren () != 1)
            throw new Exception ("The element " + name + " can only have one child with " + "this parser and has " + current.getNumChildren () +
                  ".");
         return genSimpleTree (current.getChild (0));

         // case 14
      } else if (name.equals ("treeroot") || name.equals ("lambda")) {
         Node[] nodes = stripComments (current);
         if ( nodes == null ) return null;
         current = nodes[0];
         int numChildren = current.getNumChildren ();
         int i = 0;
         lambdaFunction = true;
         lookForFunctionArgs = true;
         Node newNode = genSimpleTree (current.getChild (i));
         Node leftOverNode = null;
         while (newNode == null && i < numChildren - 1) {
            i++;
            newNode = genSimpleTree (current.getChild (i));
         }
         lookForFunctionArgs = false;
         while (leftOverNode == null && i < numChildren - 1) {
            i++;
            leftOverNode = genSimpleTree (current.getChild (i));
         }
         if (leftOverNode != null)
            throw new Exception ("TREEROOT has more than one expression.");
         lambdaFunction = false;
         return newNode;

         // case 15
      } else if (name.equals ("function")) {
         int numChildren = current.getNumChildren ();
         Node newNode = new Node (current);
         newNode.removeAllChildren ();
         for (int i = 0; i < numChildren; i++)
            newNode.addChild (genSimpleTree (current.getChild (i)));
         return newNode;

      }

      throw new Exception ("Unknown node with name " + name + ".");

   }

   /**
    * Splits a term into two components: coefficient and variable.  The variable that this term is to be split from may be
    * specified in <i>var</i>.  If <i>var</i> is specified then the value of the return result in <i>coefficient</i> is
    * not meaningful.  If <i>var</i> is not specified then it is assumed that only one id (or variable) is present in this
    * term and it will be split from the term and returned in <i>id</i>, the rest of the expression will be evaluated to
    * a Double and returned in <i>coefficient</i>.
    */
   public static SplitTerm splitTerm( Node current, String var ) throws IllegalArgumentException {
      String name = current.getSimpleName().toLowerCase();
      SplitTerm splitTerm = new SplitTerm();
      if ( name.equals("ci") ) {
         if ( var == null ) {
            splitTerm.id = current.getValue();
         } else {
            if ( var.equals(current.getValue()) ) {
               splitTerm.id = var;
            } else {
               splitTerm.termWithoutVar = current;
            }
         }
      } else if ( name.equals("cn") ) {
         splitTerm.termWithoutVar = current;
         try {
            splitTerm.coefficient = new Double(current.getValue());
         } catch (NumberFormatException e) {
            current.setError(true);
            throw new IllegalArgumentException("Invalid double '"+current.getValue()+"'",e);
         }
      } else if ( name.equals("times") ) {
         SplitTerm first  = splitTerm(current.getChild(0),var);
         SplitTerm second = splitTerm(current.getChild(1),var);
         if ( first.id != null && second.id != null ) {
            current.setError(true);
            throw new IllegalArgumentException("Only one species id per term is supported, two were found.");
         }
         splitTerm.id = first.id==null?second.id:first.id;
         if ( first.coefficient != null && second.coefficient != null ) {
            splitTerm.coefficient = new Double(first.coefficient.doubleValue()*second.coefficient.doubleValue());
         } else {
            splitTerm.coefficient = first.coefficient==null?second.coefficient:first.coefficient;
         }
         Node left       = first.termWithoutVar;
         Node right      = second.termWithoutVar;
         if ( left != null && right != null ) {
            splitTerm.termWithoutVar = new Node("times","",left,right);
         } else {
            splitTerm.termWithoutVar = left==null?right:left;
            if ( splitTerm.termWithoutVar == null ) {
               splitTerm.termWithoutVar = new Node("cn","1.0");
            }
         }
      } else if ( name.equals("plus") ) {
         SplitTerm first  = splitTerm(current.getChild(0),var);
         SplitTerm second = splitTerm(current.getChild(1),var);
         if ( first.id != null || second.id != null ) {
            current.setError(true);
            throw new IllegalArgumentException("Species ids may not appear in a sum within a product.");
         }
         if ( first.coefficient != null && second.coefficient != null ) {
            splitTerm.coefficient = new Double(first.coefficient.doubleValue()+second.coefficient.doubleValue());
         } else {
            splitTerm.coefficient = first.coefficient==null?second.coefficient:first.coefficient;
         }
         splitTerm.termWithoutVar = current;
      } else if ( name.equals("minus") && current.getNumChildren() == 1 ) {
         SplitTerm first  = splitTerm(current.getChild(0),var);
         if ( first.coefficient == null ) splitTerm.coefficient = new Double(-1.0);
         else splitTerm.coefficient = new Double(-first.coefficient.doubleValue());
         splitTerm.id = first.id;
         if ( splitTerm.id == null ) {
            splitTerm.termWithoutVar = current;
         } else {
            splitTerm.termWithoutVar = new Node("cn","-1");
         }
      } else if ( name.equals("minus") && current.getNumChildren() == 2 ) {
         SplitTerm first  = splitTerm(current.getChild(0),var);
         SplitTerm second = splitTerm(current.getChild(1),var);
         if ( first.id != null || second.id != null ) {
            current.setError(true);
            throw new IllegalArgumentException("Species ids may not appear in a sum within a product.");
         }
         if ( first.coefficient != null && second.coefficient != null ) {
            splitTerm.coefficient = new Double(first.coefficient.doubleValue()-second.coefficient.doubleValue());
         } else {
            splitTerm.coefficient = first.coefficient==null?new Double(-second.coefficient.doubleValue()):first.coefficient;
         }
         splitTerm.termWithoutVar = current;
      } else {
         current.setError(true);
         throw new IllegalArgumentException("Unsupported operation: '"+name+"'");
      }
      return splitTerm;
   }

   private void parse (InputSource is) throws SAXException, IOException, ParserConfigurationException {
      SAXParser parser = SAXParserFactory.newInstance ().newSAXParser ();
      MathMLHandler handler = new MathMLHandler ();
      parser.parse (is, handler);
      mathMLParseTree = handler.root;
   }

   /**
    * Remove the current node from it's parent and adjust the parse tree to
    * make it valid MathML parse tree.  This operates on a MathML parse tree as
    * returned by MathMLHandler.  Returns the new root node if the root changes
    * as a result of this call.
    *
    * This method is under construction.  It currently only works if the Math
    * only contains minus, plus, and times operators.
    */
   public static Node removeNode( Node current ) {
      Node root = null;
      Node parent = current.getParent();
      if ( parent == null ) return null;
      root = parent;
      int index = parent.removeChild(current);
      if ( parent.getSimpleName().equalsIgnoreCase("apply") ) {
         if ( parent.getNumChildren() == 1 ) {
            root = removeNode(parent);
         } else if ( parent.getNumChildren() == 2 ) {
            if ( parent.getChild(0).getSimpleName().equalsIgnoreCase("minus") && index == 2 ) {
               Node grandParent = parent.getParent();
               if ( grandParent == null ) {
                  root = parent.getChild(1);
                  root.setParent(null);
               } else {
                  int parentIndex = grandParent.removeChild(parent);
                  grandParent.addChild(parentIndex,parent.getChild(1));
                  root = grandParent;
               }
            }
         }
      }
      while ( root.getParent() != null ) root = root.getParent();
      return root;
   }

   /**
    * Internal routine for separating left and right side of equations.  Terms with var are put on the right and var is factored out.  Terms
    * without var are put on the left.
    */
   private Equation solve_ (Node cur, String var) throws Exception {
      String type = cur.getSimpleName ().toLowerCase ();
      Equation equation = new Equation ();
      equation.left = equation.right = null;
      if (!type.matches ("plus|minus|times|ci|cn"))
         throw new Exception ("Unsupported operation: " + type);
      else if (type.equals ("plus")) {
         if (cur.getNumChildren () != 2)
            throw new Exception ("The plus operator must have exactly 2 children" + ", it has " + cur.getNumChildren ());
         Equation child0 = solve_ (cur.getChild (0), var);
         Equation child1 = solve_ (cur.getChild (1), var);
         equation.left = child0.left;
         equation.right = child0.right;
         if (equation.left == null)
            equation.left = child1.left;
         else if (child1.left != null) {
            equation.left = new Node ("plus", "");
            equation.left.addChild (child0.left);
            equation.left.addChild (child1.left);
         }
         if (equation.right == null)
            equation.right = child1.right;
         else if (child1.right != null) {
            equation.right = new Node ("plus", "");
            equation.right.addChild (child0.right);
            equation.right.addChild (child1.right);
         }
      } else if (type.equals ("minus")) {
         int numChildren = cur.getNumChildren ();
         if (numChildren != 2 && numChildren != 1)
            throw new Exception ("The minus operator must have 1 or 2 children" + ", it has " + numChildren);
         if (numChildren == 1) {
            Equation child0 = solve_ (cur.getChild (0), var);
            if (child0.left != null) {
               equation.left = new Node ("minus", "");
               equation.left.addChild (child0.left);
            }
            if (child0.right != null) {
               equation.right = new Node ("minus", "");
               equation.right.addChild (child0.right);
            }
         } else {
            Equation child0 = solve_ (cur.getChild (0), var);
            Equation child1 = solve_ (cur.getChild (1), var);
            if (child1.left != null) {
               equation.left = new Node ("minus", "");
               if (child0.left != null)
                  equation.left.addChild (child0.left);
               equation.left.addChild (child1.left);
            } else
               equation.left = child0.left;
            if (child1.right != null) {
               equation.right = new Node ("minus", "");
               if (child0.right != null)
                  equation.right.addChild (child0.right);
               equation.right.addChild (child1.right);
            } else
               equation.right = child0.right;
         }
      } else if (type.equals ("times")) {
         SplitTerm splitTerm = splitTerm(cur,var);
         if ( splitTerm.id != null ) {
            equation.right = splitTerm.termWithoutVar;
         } else {
            equation.left = splitTerm.termWithoutVar;
         }
      } else if (type.equals ("ci")) {
         if (cur.getValue ().equals (var))
            equation.right = new Node ("cn", "1");
         else
            equation.left = cur;
      } else if (type.equals ("cn"))
         equation.left = cur;
      return equation;
   }

   /**
    * Return a parse tree with only elements that can be translated into code.
    *
    * <ul>
    * <li>
    * Remove all annotation elements and their children.
    * </li>
    * <li>
    * Ignore semantic elements and process their children.
    * </li>
    * <li>
    * Ignore math elements and process their children.
    * </li>
    * <li>
    * Ignore annotation-xml and its children if it doesn't have the encoding  "MathML-Content".  Otherwise, ignore and process its children.
    * </li>
    * </ul>
    */

   private Node [] stripComments (Node current) throws Exception {
      String name = current.getSimpleName ().toLowerCase ();

      if (name.equals ("annotation"))
         return null;
      else if (name.equals ("annotation-xml")) {
         if (!current.getAttributes ().getValue ("encoding").equalsIgnoreCase ("MathML-Content"))
            return null;
         ArrayList newNodes = new ArrayList ();
         int numChildren = current.getNumChildren ();
         int i;
         int j;
         for (i = 0; i < numChildren; i++) {
            Node tempNodes [] = stripComments (current.getChild (i));
            if ( tempNodes == null ) continue;
            for (j = 0; j < tempNodes.length; j++) {
               if (tempNodes [j] != null)
                  newNodes.add (tempNodes [j]);
            }
         }
         if (newNodes.size () == 0)
            return null;
         return (Node []) newNodes.toArray (new Node [0]);

      } else if (name.equals ("treeroot")) {
         int numChildren;
         Node newNode = new Node (current);
         newNode.removeAllChildren ();
         numChildren = current.getNumChildren ();
         int i;
         int j;
         for (i = 0; i < numChildren; i++) {
            Node tempNodes [] = stripComments (current.getChild (i));
            if ( tempNodes == null ) continue;
            for (j = 0; j < tempNodes.length; j++) {
               if (tempNodes [j] != null)
                  newNode.addChild (tempNodes [j]);
            }
         }
         if (newNode.getNumChildren () == 0)
            return null;
         Node retNode [] = {newNode};
         return retNode;

      } else if (name.matches ("semantics|math")) {
         ArrayList newNodes = new ArrayList ();
         int numChildren = current.getNumChildren ();
         int i;
         int j;
         for (i = 0; i < numChildren; i++) {
            Node tempNodes [] = stripComments (current.getChild (i));
            if ( tempNodes == null ) continue;
            for (j = 0; j < tempNodes.length; j++) {
               if (tempNodes [j] != null)
                  newNodes.add (tempNodes [j]);
            }
         }
         if (newNodes.size () == 0)
            return null;
         return (Node []) newNodes.toArray (new Node [0]);

      } else {
         Node newNode = new Node (current);
         newNode.removeAllChildren ();
         int numChildren = current.getNumChildren ();
         int i;
         int j;
         for (i = 0; i < numChildren; i++) {
            Node tempNodes [] = stripComments (current.getChild (i));
            if ( tempNodes == null ) continue;
            for (j = 0; j < tempNodes.length; j++) {
               if (tempNodes [j] != null)
                  newNode.addChild (tempNodes [j]);
            }
         }
         Node retNode [] = {newNode};
         return retNode;
      }
   }

   public void searchAndRecordIds (Set setOfUsedIds) {
      Set bvarSet = new HashSet ();
      // invalidate simpleTree so that it gets regenerated next time around
      simpleTree = null; 

      // ignore root -- root acts as a container
      Node node = mathMLParseTree.getChild (0);

      while (node != null) {
         String nodeType = node.getQName ();
         // Populate the bvarSet if you encounter a lambda node
         if (nodeType.endsWith ("lambda")) {
            int numberOfBvars = node.getNumChildren () - 1;
            for (int i = 0; i < numberOfBvars; i++) {
               Node bvarNode = node.getChild (i);
               if (bvarNode.getQName ().endsWith ("bvar")) 
                  bvarSet.add (bvarNode.getChild (0).getValue ());
            }
         }

         // Update the id of the node 
         if (nodeType != null && nodeType.endsWith ("ci")) {
            if (!bvarSet.contains (node.getValue ())) 
               setOfUsedIds.add (node.getValue ());
         }

         if (node.getNumChildren () > 0) 
            node = node.getChild (0);
         else {   
            // find the parent level
            Node parent = node.getParent ();
            List children = parent.getListOfChildren ();
            int siblingIndex = children.indexOf (node) + 1;
            while (siblingIndex == children.size () && node != mathMLParseTree) {
               // use child-parent link to get to the parent level
               node = node.getParent ();
               parent = node.getParent ();
               if (parent == mathMLParseTree)
                  return;
               children = parent.getListOfChildren ();
               siblingIndex = children.indexOf (node) + 1;
            }
            node = parent.getChild (siblingIndex);
         }
      }
      bvarSet.clear ();
   }

   public void searchAndReplaceIds (Map oldToNewIds) {
      Set bvarSet = new HashSet ();
      // invalidate simpleTree so that it gets regenerated next time around
      simpleTree = null; 

      // ignore root -- root acts as a container
      Node node = mathMLParseTree.getChild (0);

      while (node != null) {
         String nodeType = node.getQName ();
         // Populate the bvarSet if you encounter a lambda node
         if (nodeType.endsWith ("lambda")) {
            int numberOfBvars = node.getNumChildren () - 1;
            for (int i = 0; i < numberOfBvars; i++) {
               Node bvarNode = node.getChild (i);
               if (bvarNode.getQName ().endsWith ("bvar")) {
                  //This previous version did not get the actual ci values for the bvar
                  //parameters in the functionDefinitions
                  //bvarSet.add (node.getValue ());
                  bvarSet.add (bvarNode.getChild (0).getValue ());
               }
            }
         }

         // Update the id of the node 
         if (nodeType != null && nodeType.endsWith ("ci")) {
            String nodeId = (String)oldToNewIds.get (node.getValue ());
            if (nodeId != null && !bvarSet.contains (node.getValue ())) 
               node.setValue (nodeId);
         }

         if (node.getNumChildren () > 0) 
            node = node.getChild (0);
         else {   
            // find the parent level
            Node parent = node.getParent ();
            List children = parent.getListOfChildren ();
            int siblingIndex = children.indexOf (node) + 1;
            while (siblingIndex == children.size () && node != mathMLParseTree) {
               // use child-parent link to get to the parent level
               node = node.getParent ();
               parent = node.getParent ();
               if (parent == mathMLParseTree)
                  return;
               children = parent.getListOfChildren ();
               siblingIndex = children.indexOf (node) + 1;
            }
            node = parent.getChild (siblingIndex);
         }
      }
      bvarSet.clear ();
   }

   public String toString () {
      mathML = new StringBuffer ();
      try {
         generateMathML (mathMLParseTree.getChild (0));
      } catch (Exception e) {
         e.printStackTrace ();
      }
      // Remove the last newline/carriage return
      mathML.delete (mathML.length () - 1, mathML.length ());
      return mathML.toString ();
   }

   private void generateMathML (Node cur) throws Exception{
      String nodeType = cur.getQName ();
      if (nodeType == null)
         throw new Exception ("Node type was null");

      String nodeValue = cur.getValue ().trim ();
      if (nodeValue == null)
         nodeValue = "";

      // Special format for bvars, all on same line.
      // Every other node needs indenting
      if (!cur.getParent ().getQName ().endsWith ("bvar"))
         mathML.append (tabs);

      mathML.append ("<" + nodeType);

      // Loop through node attributes
      //int numberOfAttributes = cur.getAttributes ().getLength ();
      //for (int i = 0; i < numberOfAttributes; i++) {
      //   mathML.append (" " + cur.getAttributes ().getQName (i) + " " + cur.getAttributes ().getType (i) + " " + cur.getAttributes ().getValue (i));
      //}

      // Node has no children and no value, close the tag off and return
      if (cur.getNumChildren () == 0 && nodeValue.equals ("")) {
         mathML.append ("/>\n");
         return;
      } 

      // CI, CN case. Print out value and then close tag on one line
      if (cur.getNumChildren () == 0) {
         mathML.append (">" + nodeValue + "</" + nodeType + ">\n");
         return;
      }

      // All other cases print out children and node value.
      // Note the bvars will have no newline but instead
      // will be printed out on one line
      mathML.append (">" + nodeValue);
      if (!nodeType.endsWith ("bvar"))
         mathML.append ("\n");

      tabs.append ("    ");
      for (int i = 0; i < cur.getNumChildren (); i++) 
         generateMathML (cur.getChild (i));
      tabs.delete (tabs.length () - 4, tabs.length ());

      // Ensure that bvars are formatted all on one line
      if (nodeType.endsWith ("bvar")) {
         mathML.delete (mathML.length () - 1, mathML.length ());
         mathML.append ("</" + cur.getQName () + ">\n");
      } 
      // For all non-bvar cases 
      else  
         mathML.append (tabs + "</" + cur.getQName () + ">\n");
   }
}
