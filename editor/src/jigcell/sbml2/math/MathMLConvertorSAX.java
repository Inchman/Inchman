package jigcell.sbml2.math;

import java.io.PrintStream;
import java.io.StringReader;
import java.util.Hashtable;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import org.xml.sax.Attributes;
import org.xml.sax.ErrorHandler;
import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;
import org.xml.sax.XMLReader;

/**
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * @author Marc Vass
 */

public class MathMLConvertorSAX extends DefaultHandler {
   private static Hashtable functions = new Hashtable ();
   private static String mathML;
   private boolean isLambda;
   private boolean isRatelaw;
   private boolean lambdaGiven;
   private Hashtable hashTable = new Hashtable ();
   private int BVARMODE = 16;
   private int CIMODE = 2;
   private int CNMODE = 1;
   private int CSYMBOLMODE = 4;
   private int DELAYMODE = 8;
   private int ENOTATION = 32;
   private int ENOTATIONSEP = 256;
   private int INTEGER = 64;
   private int INTEGERSEP = 512;
   private int NORMAL = 0;
   private int RATIONAL = 128;
   private int RATIONALSEP = 1024;
   private int argumentNum = 1;
   private int mode = NORMAL;
   private int numArgumentsExpected = 0;
   private SBMLNode currentNode;
   private SBMLNode root;
   private String type = "";

   // Error handler to report errors and warnings
   private static class MyErrorHandler implements ErrorHandler {

      /**
       * Error handler output goes here
       */

      private PrintStream out;

      public void error (SAXParseException spe) throws SAXException {
         String message = "Error: " + getParseExceptionInfo (spe);
         throw new SAXException (message);
      }

      public void fatalError (SAXParseException spe) throws SAXException {
         String message = "Fatal Error: " + getParseExceptionInfo (spe);
         throw new SAXException (message);
      }

      // The following methods are standard SAX ErrorHandler methods.
      // See SAX documentation for more info.
      public void warning (SAXParseException spe) throws SAXException {
         out.println ("Warning: " + getParseExceptionInfo (spe));
      }

      MyErrorHandler (PrintStream out) {
         this.out = out;
      }

      /**
       * Returns a string describing parse exception details
       */

      private String getParseExceptionInfo (SAXParseException spe) {
         String systemId = spe.getSystemId ();
         if (systemId == null)
            systemId = "null";
         String info = "URI=" + systemId + " Line=" + spe.getLineNumber () + ": " + spe.getMessage ();
         return info;
      }
   }

   public MathMLConvertorSAX () {}

   /**
    * Creates a new instance of MathMLConvertorSAX
    */

   public MathMLConvertorSAX (String mathML, boolean isLambda, boolean isRatelaw) throws Exception {
      MathMLConvertorSAX.mathML = mathML;

      this.isRatelaw = isRatelaw;
      root = null;

      currentNode = null;

      this.isLambda = isLambda;

      if (isLambda) {
         if (mathML.startsWith ("<lambda") && !mathML.endsWith ("</lambda>"))
            mathML = mathML + "</lambda>";
         else if (mathML.startsWith ("<math:lambda") && !mathML.endsWith ("</math:lambda>"))
            mathML = mathML + "</math:lambda>";
      }
      if (mathML.startsWith ("null"))
         mathML = mathML.substring (4);

      // Create a JAXP SAXParserFactory and configure it
      SAXParserFactory spf = SAXParserFactory.newInstance ();

      // Set namespaceAware to true to get a parser that corresponds to
      // the default SAX2 namespace feature setting.  This is necessary
      // because the default value from JAXP 1.0 was defined to be false.
      spf.setNamespaceAware (false);

      // Validation part 1: set whether validation is on
      spf.setValidating (false);

      // Create a JAXP SAXParser
      SAXParser saxParser = spf.newSAXParser ();

      // Get the encapsulated SAX XMLReader
      XMLReader xmlReader = saxParser.getXMLReader ();

      // Set the ContentHandler of the XMLReader
      xmlReader.setContentHandler (this);

      // Set an ErrorHandler before parsing
      xmlReader.setErrorHandler (new MyErrorHandler (System.err));

      // Tell the XMLReader to parse the XML document
      InputSource source = new InputSource (new StringReader (mathML));
      xmlReader.parse (source);
   }

   public void characters (char ch [], int start, int length) throws SAXException {
      try {
         String s = new String (ch, start, length);

         if (s.trim ().length () != 0 || s == null) {
            String trimmed = s.trim ();
            if (mode == CNMODE || (mode & CNMODE) == CNMODE) {
               if (mode == CNMODE + RATIONALSEP)
                  trimmed = type + "/" + trimmed;
               else if (mode == CNMODE + ENOTATIONSEP)
                  trimmed = type + "e" + trimmed;
               if (mode == CNMODE + RATIONAL) {
                  type = trimmed;
                  mode = CNMODE + RATIONALSEP;
               } else if (mode == CNMODE + ENOTATION) {
                  type = trimmed;
                  mode = CNMODE + ENOTATIONSEP;
               } else if (isLambda && !lambdaGiven) {
                  lambdaGiven = true;

                  SBMLNode n = new SBMLNode (currentNode);
                  n.setElement (new MathElement (trimmed));
                  currentNode.addChild (n);
               } else if (currentNode != null) {
                  SBMLNode n = new SBMLNode (currentNode);
                  n.setElement (new MathElement (trimmed));
                  currentNode.addChild (n);
               } else {
                  root = new SBMLNode (new MathElement ()); // the first operator
                  currentNode = root;
                  MathElement m = new MathElement (trimmed);
                  currentNode.setElement (m);
               }
            } else if (mode == CIMODE || mode == BVARMODE + CIMODE) {
               if (mode == BVARMODE + CIMODE) {
                  if (isLambda) {
                     String argument = trimmed;

                     // hashTable.put(argument,"A" + ++argumentNum);
                     String argument2 = "A" + argumentNum;
                     argumentNum++;
                     if (!isRatelaw)
                        hashTable.put (argument, argument2);
                     else
                        hashTable.put (argument, argument);
                  }
               } else if (isLambda && !lambdaGiven) {
                  String argument = (String) hashTable.get (trimmed);
                  if (currentNode != null) {
                     SBMLNode n = new SBMLNode (currentNode);
                     n.setElement (new MathElement (argument));
                     currentNode.addChild (n);
                  } else {
                     root = new SBMLNode (new MathElement (argument));
                     currentNode = root;
                  }
               } else {
                  if (functions.containsKey (trimmed))
                     currentNode.setElement (new MathElement (trimmed));
                  else if (currentNode != null) {
                     if (currentNode.getElement () == null) // || currentNode.getElement().toString()==null){

                        currentNode.setElement (new MathElement (trimmed));
                     else {
                        SBMLNode n = new SBMLNode (currentNode);
                        n.setElement (new MathElement (trimmed));
                        currentNode.addChild (n);
                     }
                  } else {
                     currentNode = new SBMLNode (new MathElement ());
                     root = currentNode;
                     currentNode.setElement (new MathElement (trimmed));
                  }
               }
            }
             // end cimode
         }
      } catch (Exception e) {
         e.printStackTrace ();
      }
   }

   // Parser calls this once after parsing a document
   public void endDocument () throws SAXException {}

   public void endElement (String uri, String localName, String qName) throws SAXException {
      if (qName.equals ("math:apply") || qName.equals ("apply"))
         currentNode = currentNode.getParent ();
      else if (qName.equals ("math:piecewise") || qName.equals ("piecewise"))
         currentNode = currentNode.parent;
      else if (qName.equals ("math:piece") || qName.equals ("piece"))
         currentNode = currentNode.parent;
      else if (qName.equals ("math:otherwise") || qName.equals ("otherwise"))
         currentNode = currentNode.parent;
      else if (qName.equals ("math:bvar") || qName.equals ("bvar"))
         mode = NORMAL;
   }

   public Hashtable getFunctions () {
      return functions;
   }

   /**
    * Getter for property root.
    *
    * @return Value of property root.
    */

   public jigcell.sbml2.math.SBMLNode getRoot () {
      return root;
   }

   public void insertFunction (String function) {
      functions.put (function, function);
   }

   public void reset () {
      functions.clear ();
   }

   // Parser calls this once at the beginning of a document
   public void startDocument () throws SAXException {}

   // Parser calls this for each element in a document
   public void startElement (String namespaceURI, String localName, String qName, Attributes atts) throws SAXException {
      try {
         if (qName.equals ("math:cn") || qName.equals ("cn")) {
            mode = CNMODE;
            if (atts.getIndex ("type") != -1) {
               if (atts.getValue ("type").equals ("e-notation"))
                  mode = CNMODE + ENOTATION;
               else if (atts.getValue ("type").equals ("rational"))
                  mode = CNMODE + RATIONAL;
            }
         } else if (qName.equals ("math:ci") || qName.equals ("ci")) {
            if (mode != CIMODE + BVARMODE)
               mode = CIMODE;
            if (atts.getValue ("type") != null) {}
         } else if (qName.equals ("math:csymbol") || qName.equals ("csymbol") || qName.equals ("/csymbol") || qName.equals ("/math:csymbol")) {
            mode = CSYMBOLMODE;
            if (atts.getValue ("definitionURL").equals ("http://www.sbml.org/sbml/symbols/time"))
               currentNode.setElement (new MathElement ("@time"));
            else if (atts.getValue ("definitionURL").equals ("http://www.sbml.org/sbml/symbols/delay"))
               currentNode.setElement (new MathElement ("@delay"));
         } else if (qName.equals ("math:apply") || qName.equals ("apply")) {
            if (isLambda) {
               if (root == null) {
                  root = new SBMLNode (new MathElement ()); // the first operator
                  currentNode = root;
               } else { // not the root
                  SBMLNode n = new SBMLNode (currentNode);
                  currentNode.addChild (n);
                  currentNode = n; // make this node the current one
               }
            } else {
               if (root == null) {
                  root = new SBMLNode (new MathElement ()); // the first operator will go here
                  currentNode = root;
               } else { // not the root
                  SBMLNode n = new SBMLNode (currentNode);
                  currentNode.addChild (n);
                  currentNode = n; // make this node the current one
               }
            }
         } else if (qName.equals ("math:piecewise") || qName.equals ("piecewise")) {
            if (root == null) {
               root = new SBMLNode (new MathElement ()); // the first operator will go here
               currentNode = root;
            } else { // not the root
               SBMLNode n = new SBMLNode (currentNode);
               currentNode.addChild (n);
               currentNode = n; // make this node the current one
            }
            currentNode.setElement (new MathElement ("@piecewise"));
         } else if (qName.equals ("math:piece") || qName.equals ("piece")) {
            SBMLNode n = new SBMLNode (currentNode);
            currentNode.addChild (n);
            currentNode = n; // make this node the current one
            currentNode.setElement (new MathElement ("@piece"));
         } else if (qName.equals ("math:otherwise") || qName.equals ("otherwise")) {
            SBMLNode n = new SBMLNode (currentNode);
            currentNode.addChild (n);
            currentNode = n; // make this node the current one
            currentNode.setElement (new MathElement ("@otherwise"));
         }
         else if (qName.equals ("math:eq") || qName.equals ("eq"))
            currentNode.setElement (new MathElement ("@eq"));
         else if (qName.equals ("math:neq") || qName.equals ("neq"))
            currentNode.setElement (new MathElement ("@neq"));
         else if (qName.equals ("math:gt") || qName.equals ("gt"))
            currentNode.setElement (new MathElement ("@gt"));
         else if (qName.equals ("math:lt") || qName.equals ("lt"))
            currentNode.setElement (new MathElement ("@lt"));
         else if (qName.equals ("math:geq") || qName.equals ("geq"))
            currentNode.setElement (new MathElement ("@geq"));
         else if (qName.equals ("math:leq") || qName.equals ("leq"))
            currentNode.setElement (new MathElement ("@leq"));
         else if (qName.equals ("math:plus") || qName.equals ("plus"))
            currentNode.setElement (new MathElement ("@plus"));
         else if (qName.equals ("math:minus") || qName.equals ("minus"))
            currentNode.setElement (new MathElement ("@minus"));
         else if (qName.equals ("math:times") || qName.equals ("times"))
            currentNode.setElement (new MathElement ("@times"));
         else if (qName.equals ("math:divide") || qName.equals ("divide"))
            currentNode.setElement (new MathElement ("@divide"));
         else if (qName.equals ("math:power") || qName.equals ("power"))
            currentNode.setElement (new MathElement ("@power"));
         else if (qName.equals ("math:root") || qName.equals ("root"))
            currentNode.setElement (new MathElement ("@root"));
         else if (qName.equals ("math:abs") || qName.equals ("abs"))
            currentNode.setElement (new MathElement ("@abs"));
         else if (qName.equals ("math:exp") || qName.equals ("exp"))
            currentNode.setElement (new MathElement ("@exp"));
         else if (qName.equals ("math:ln") || qName.equals ("ln"))
            currentNode.setElement (new MathElement ("@ln"));
         else if (qName.equals ("math:log") || qName.equals ("log"))
            currentNode.setElement (new MathElement ("@log"));
         else if (qName.equals ("math:floor") || qName.equals ("floor"))
            currentNode.setElement (new MathElement ("@floor"));
         else if (qName.equals ("math:ceiling") || qName.equals ("ceiling"))
            currentNode.setElement (new MathElement ("@ceiling"));
         else if (qName.equals ("math:factorial") || qName.equals ("factorial"))
            currentNode.setElement (new MathElement ("@factorial"));
         else if (qName.equals ("math:and") || qName.equals ("and"))
            currentNode.setElement (new MathElement ("@and"));
         else if (qName.equals ("math:or") || qName.equals ("or"))
            currentNode.setElement (new MathElement ("@or"));
         else if (qName.equals ("math:xor") || qName.equals ("xor"))
            currentNode.setElement (new MathElement ("@xor"));
         else if (qName.equals ("math:not") || qName.equals ("not"))
            currentNode.setElement (new MathElement ("@not"));
         else if (qName.equals ("math:degree") || qName.equals ("degree"))
            ;
         else if (qName.equals ("math:bvar") || qName.equals ("bvar"))
            mode = BVARMODE + CIMODE;
         else if (qName.equals ("math:logbase") || qName.equals ("logbase"))
            ;
         else if (qName.equals ("math:sin") || qName.equals ("sin"))
            currentNode.setElement (new MathElement ("@sin"));
         else if (qName.equals ("math:cos") || qName.equals ("cos"))
            currentNode.setElement (new MathElement ("@cos"));
         else if (qName.equals ("math:tan") || qName.equals ("tan"))
            currentNode.setElement (new MathElement ("@tan"));
         else if (qName.equals ("math:sec") || qName.equals ("sec"))
            currentNode.setElement (new MathElement ("@sec"));
         else if (qName.equals ("math:csc") || qName.equals ("csc"))
            currentNode.setElement (new MathElement ("@csc"));
         else if (qName.equals ("math:cot") || qName.equals ("cot"))
            currentNode.setElement (new MathElement ("@cot"));
         else if (qName.equals ("math:sinh") || qName.equals ("sinh"))
            currentNode.setElement (new MathElement ("@sinh"));
         else if (qName.equals ("math:cosh") || qName.equals ("cosh"))
            currentNode.setElement (new MathElement ("@cosh"));
         else if (qName.equals ("math:tanh") || qName.equals ("tanh"))
            currentNode.setElement (new MathElement ("@tanh"));
         else if (qName.equals ("math:sech") || qName.equals ("sech"))
            currentNode.setElement (new MathElement ("@sech"));
         else if (qName.equals ("math:csch") || qName.equals ("csch"))
            currentNode.setElement (new MathElement ("@csch"));
         else if (qName.equals ("math:coth") || qName.equals ("coth"))
            currentNode.setElement (new MathElement ("@coth"));
         else if (qName.equals ("math:arcsin") || qName.equals ("arcsin"))
            currentNode.setElement (new MathElement ("@arcsin"));
         else if (qName.equals ("math:arccos") || qName.equals ("arccos"))
            currentNode.setElement (new MathElement ("@arccos"));
         else if (qName.equals ("math:arctan") || qName.equals ("arctan"))
            currentNode.setElement (new MathElement ("@arctan"));
         else if (qName.equals ("math:arccosh") || qName.equals ("arccosh"))
            currentNode.setElement (new MathElement ("@arccosh"));
         else if (qName.equals ("math:arccot") || qName.equals ("arccot"))
            currentNode.setElement (new MathElement ("@arccot"));
         else if (qName.equals ("math:arccoth") || qName.equals ("arccoth"))
            currentNode.setElement (new MathElement ("@arccoth"));
         else if (qName.equals ("math:arccsc") || qName.equals ("arccsc"))
            currentNode.setElement (new MathElement ("@arccsc"));
         else if (qName.equals ("math:arccsch") || qName.equals ("arccsch"))
            currentNode.setElement (new MathElement ("@arccsch"));
         else if (qName.equals ("math:arcsec") || qName.equals ("arcsec"))
            currentNode.setElement (new MathElement ("@arcsec"));
         else if (qName.equals ("math:arcsech") || qName.equals ("arcsech"))
            currentNode.setElement (new MathElement ("@arcsech"));
         else if (qName.equals ("math:arcsinh") || qName.equals ("arcsinh"))
            currentNode.setElement (new MathElement ("@arcsinh"));
         else if (qName.equals ("math:arctanh") || qName.equals ("arctanh"))
            currentNode.setElement (new MathElement ("@arctanh"));
         else if (qName.equals ("math:true") || qName.equals ("true"))
            currentNode.setElement (new MathElement ("@true"));
         else if (qName.equals ("math:false") || qName.equals ("false"))
            currentNode.setElement (new MathElement ("@false"));
         else if (qName.equals ("math:notanumber") || qName.equals ("notanumber"))
            currentNode.setElement (new MathElement ("@notanumber"));
         else if (qName.equals ("math:pi") || qName.equals ("pi"))
            currentNode.setElement (new MathElement ("@pi"));
         else if (qName.equals ("math:infinity") || qName.equals ("infinity"))
            currentNode.setElement (new MathElement ("@infinity"));
         else if (qName.equals ("math:exponentiale") || qName.equals ("exponentiale"))
            currentNode.setElement (new MathElement ("@exponentiale"));
         else if (qName.equals ("math:semantics") || qName.equals ("semantics")) {

            // throw new Exception("JigCell does not support 'semantics'");
         } else if (qName.equals ("math:annotation") || qName.equals ("annotation")) {

            //throw new Exception("JigCell does not support 'annotation'");
         } else if (qName.equals ("math:annotation-xml") || qName.equals ("annotation-xml")) {

            // throw new Exception("JigCell does not support 'annotation-xml'");
         }
      } catch (Exception e) {
         e.printStackTrace ();
      }
   }
}
