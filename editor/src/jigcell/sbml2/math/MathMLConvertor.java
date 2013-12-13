package jigcell.sbml2.math;

import java.util.Hashtable;

/**
 * A class for converting MathML to normal math. This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more
 * details.
 *
 * @author Marc Vass
 * @author Nicholas Allen
 */

public class MathMLConvertor {
   private static MathMLConvertorSAX saxConvertor;
   static boolean firstPiece = true;
   private Hashtable associative = new Hashtable ();
   private Hashtable order = new Hashtable ();

   static {
      saxConvertor = new MathMLConvertorSAX ();
   }

   public static void insertFunction (String function) {
      saxConvertor.insertFunction (function);
   }

   public static void reset () {
      saxConvertor.reset ();
   }

   /**
    * Gets the normal math representation of the string.
    *
    * @param mathML The MathML string
    * @param isLamda Is the MathML String a lambda string?
    *
    * @return The normal representation
    */

   public static String toNormalMath (String mathML, boolean isLambda, boolean isRatelaw) throws Exception {
      firstPiece = true;
      saxConvertor = new MathMLConvertorSAX (mathML, isLambda, isRatelaw);

      // run down the tree
      StringBuffer normalMath = new StringBuffer ();
      SBMLNode runningNode = saxConvertor.getRoot ();
      if (saxConvertor.getRoot () == null)
         return "";

      normalMath.append (treeOutput (runningNode));
      if (normalMath.toString ().trim ().length () == 0)
         return "Error!";
      return normalMath.toString ();
   }

   private static String treeOutput (SBMLNode rootNode) {
      boolean operator = false;
      if (rootNode != null && rootNode.element != null && rootNode.element.toString () != null && rootNode.element.toString ().startsWith ("@"))
         operator = true;
      if (operator) {
         String temp = rootNode.element.toString ().substring (1, rootNode.element.toString ().length ());
         if (temp.equals ("piecewise")) {
            StringBuffer b = new StringBuffer ("(");
            for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
               if (((SBMLNode) rootNode.children.get (i)).element.toString ().equals ("piece") && i != 0)
                  b.append ("else");
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i + 1 != ((SBMLNode) rootNode).children.size ())
                  b.append ("\n");
            }
            b.append (")");
            return b.toString ();
         } else if (temp.equals ("piece")) {
            if (firstPiece) {
               firstPiece = false;
               return "if (" + treeOutput ((SBMLNode) rootNode.children.get (1)) + ") then " + "(" +
                  treeOutput ((SBMLNode) rootNode.children.get (0)) + ")\n"; //\nendif";
            } else
               return "elseif (" + treeOutput ((SBMLNode) rootNode.children.get (1)) + ") then " + "(" +
                  treeOutput ((SBMLNode) rootNode.children.get (0)) + ")\n"; //\nendif";
         } else if (temp.equals ("otherwise"))
            return "else (" + treeOutput ((SBMLNode) rootNode.children.get (0)) + ")\n"; //endif";
         if (temp.equals ("eq")) {
            StringBuffer b = new StringBuffer ();
            for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i + 1 != ((SBMLNode) rootNode).children.size ())
                  b.append ("==");
            }

            return b.toString ();
         } else if (temp.equals ("neq")) { // binary
            if (((SBMLNode) rootNode).children.size () == 2)
               return treeOutput ((SBMLNode) rootNode.children.get (0)) + "!=" + treeOutput ((SBMLNode) rootNode.children.get (1));
         } else if (temp.equals ("gt")) { // n-ary
            StringBuffer b = new StringBuffer ();
            for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i + 1 != ((SBMLNode) rootNode).children.size ())
                  b.append (">");
            }

            return b.toString ();
         } else if (temp.equals ("lt")) { // n-ary
            StringBuffer b = new StringBuffer ();
            for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i + 1 != ((SBMLNode) rootNode).children.size ())
                  b.append ("<");
            }

            return b.toString ();
         } else if (temp.equals ("geq")) { // n-ary
            StringBuffer b = new StringBuffer ();
            for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i + 1 != ((SBMLNode) rootNode).children.size ())
                  b.append (">=");
            }

            return b.toString ();
         } else if (temp.equals ("leq")) { // n-ary
            StringBuffer b = new StringBuffer ();
            for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i + 1 != ((SBMLNode) rootNode).children.size ())
                  b.append ("<=");
            }

            return b.toString ();
         } else if (temp.equals ("plus")) { // n-ary
            StringBuffer b = new StringBuffer ("(");
            for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i + 1 != ((SBMLNode) rootNode).children.size ())
                  b.append ("+");
            }
            b.append (")");
            return b.toString ();
         } else if (temp.equals ("minus")) {
            if (((SBMLNode) rootNode).children.size () == 2)
               return "(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "-" + treeOutput ((SBMLNode) rootNode.children.get (1)) + ")";
            else // unary

               return "(-" + treeOutput ((SBMLNode) rootNode.children.get (0)) + ")";
         } else if (temp.equals ("times")) { // n-ary
            StringBuffer b = new StringBuffer ("(");
            for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i + 1 != ((SBMLNode) rootNode).children.size ())
                  b.append ("*");
            }
            b.append (")");
            return b.toString ();
         } else if (temp.equals ("divide")) {
            if (((SBMLNode) rootNode).children.size () == 2) {
               String child1,child2;
               child1 = treeOutput ((SBMLNode) rootNode.children.get (0));
               child2 = treeOutput ((SBMLNode) rootNode.children.get (1));
               return "(" + child1 + "/" + child2 + ")";
            } else if (((SBMLNode) rootNode).children.size () == 3) {
               String child1,child2,child3;
               child1 = treeOutput ((SBMLNode) rootNode.children.get (0));
               child2 = treeOutput ((SBMLNode) rootNode.children.get (1));
               child3 = treeOutput ((SBMLNode) rootNode.children.get (2));
               return "(" + child1 + "/" + child2 + ")";
            }
         } else if (temp.equals ("power")) {
            if (((SBMLNode) rootNode).children.size () == 2)
               return "(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "^" + treeOutput ((SBMLNode) rootNode.children.get (1)) + ")";
         } else if (temp.equals ("root")) {
            if (((SBMLNode) rootNode).children.size () == 2) {
               if ((((SBMLNode) rootNode.children.get (0))).element.toString ().equals ("2"))
                  return "(" + "sqrt" + "(" + treeOutput ((SBMLNode) rootNode.children.get (1)) + "))";
               else
                  return "(" + "root" + "(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "," +
                     treeOutput ((SBMLNode) rootNode.children.get (1)) + "))";
            }
         } else if (temp.equals ("abs"))
            return "(abs(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("exp"))
            return "(exp(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("ln"))
            return "(ln(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("log")) {
            if (((SBMLNode) rootNode).children.size () == 2) {
               if ((((SBMLNode) rootNode.children.get (0))).element.toString ().equals ("10"))
                  return "(" + "log10" + "(" + treeOutput ((SBMLNode) rootNode.children.get (1)) + "))";
               else
                  return "(" + "log" + "(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "," +
                     treeOutput ((SBMLNode) rootNode.children.get (1)) + "))";
            }
         } else if (temp.equals ("floor"))
            return "(floor(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("ceiling"))
            return "(ceiling(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("factorial"))
            return "(factorial" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("and")) { // n-ary
            StringBuffer b = new StringBuffer ();
            for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i + 1 != ((SBMLNode) rootNode).children.size ())
                  b.append (" AND ");
            }

            return b.toString ();
         } else if (temp.equals ("or")) { // n-ary
            StringBuffer b = new StringBuffer ();
            for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i + 1 != ((SBMLNode) rootNode).children.size ())
                  b.append (" OR ");
            }

            return b.toString ();
         } else if (temp.equals ("xor")) { // n-ary
            StringBuffer b = new StringBuffer ();
            for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i + 1 != ((SBMLNode) rootNode).children.size ())
                  b.append (" XOR ");
            }

            return b.toString ();
         } else if (temp.equals ("not"))
            return "(!" + treeOutput ((SBMLNode) rootNode.children.get (0)) + ")";
         else if (temp.equals ("logbase")) {}
         else if (temp.equals ("sin"))
            return "(sin(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("cos"))
            return "(cos(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("tan"))
            return "(tan(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("sec"))
            return "(sec(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("csc"))
            return "(csc(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("cot"))
            return "(cot(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("sinh"))
            return "(sinh(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("cosh"))
            return "(cosh(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("tanh"))
            return "(tanh(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("sech"))
            return "(sech(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("csch"))
            return "(csch(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("coth"))
            return "(coth(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arcsin"))
            return "(arcsin(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arccos"))
            return "(arccos(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arctan"))
            return "(arctan(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arcsinh"))
            return "(arcsinh(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arccosh"))
            return "(arccosh(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arctanh"))
            return "(arctanh(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arccot"))
            return "(arccot(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arccoth"))
            return "(arccoth(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arccsc"))
            return "(arccsc(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arccsch"))
            return "(arccsch(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arcsec"))
            return "(arcsec(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("arcsech"))
            return "(arcsech(" + treeOutput ((SBMLNode) rootNode.children.get (0)) + "))";
         else if (temp.equals ("true"))
            return "(true)";
         else if (temp.equals ("false"))
            return "(false)";
         else if (temp.equals ("notanumber"))
            return "(nan)";
         else if (temp.equals ("pi"))
            return "(pi)";
         else if (temp.equals ("infinity"))
            return "(infinity)";
         else if (temp.equals ("exponentiale"))
            return "(exponentiale)";
         else if (temp.equals ("semantics")) {}
         else if (temp.equals ("annotation")) {}
         else if (temp.equals ("annotation-xml")) {}
         else if (temp.equals ("time"))
            return "@time";
         else if (temp.equals ("delay")) { //delay
            StringBuffer b = new StringBuffer ();
            for (int i = 0; i < 2; i++) {
               b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
               if (i == 0)
                  b.append (",");
            }
            return "delay(" + b.toString () + ")";
         }
      } else if (rootNode.element.toString () != null && saxConvertor.getFunctions ().containsKey (rootNode.element.toString ())) {
         StringBuffer b = new StringBuffer (rootNode.element.toString () + "(");
         for (int i = 0; i < ((SBMLNode) rootNode).children.size (); i++) {
            b.append (treeOutput ((SBMLNode) rootNode.children.get (i)));
            if (i + 1 != ((SBMLNode) rootNode).children.size ())
               b.append (",");
         }
         b.append (")");
         return b.toString ();
      } else
         return rootNode.element.toString ();
      return "";
   }

   /**
    * Creates a new instance of MathMLConvertor
    */

   public MathMLConvertor () {
      saxConvertor = new MathMLConvertorSAX ();
      firstPiece = true;
      order = new Hashtable (30);
      order.put ("power", new Integer (1));
      order.put ("times", new Integer (2));
      order.put ("divide", new Integer (2));
      order.put ("plus", new Integer (3));
      order.put ("minus", new Integer (3));
      order.put ("eq", new Integer (4));
      order.put ("neq", new Integer (4));
      order.put ("gt", new Integer (4));
      order.put ("lt", new Integer (4));
      order.put ("geq", new Integer (4));
      order.put ("leq", new Integer (4));
      order.put ("not", new Integer (5));
      order.put ("and", new Integer (6));
      order.put ("or", new Integer (7));

      associative = new Hashtable (12);
      associative.put (new Integer (1), new Boolean (true));
      associative.put (new Integer (2), new Boolean (true));
      associative.put (new Integer (3), new Boolean (true));
      associative.put (new Integer (4), new Boolean (false));
      associative.put (new Integer (5), new Boolean (true));
      associative.put (new Integer (6), new Boolean (true));
      associative.put (new Integer (7), new Boolean (true));
   }
}
