package jigcell.sbml2;

import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import javax.xml.parsers.SAXParserFactory;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.SAXException;

public abstract class XMLReader extends DefaultHandler {
   protected final static int NORMAL = 0;

   protected final Map namespacePrefixToURI;
   protected int mode;
   protected final LinkedList stack;
   protected final PrintStream errorStream;
   protected StringBuffer localText;

   protected static void read (XMLReader handler, Reader reader) throws IOException {
      SAXParserFactory factory = SAXParserFactory.newInstance ();
      factory.setNamespaceAware (true);
      factory.setValidating (false);
      org.xml.sax.XMLReader xmlReader;
      try {
         xmlReader = factory.newSAXParser ().getXMLReader ();
      } catch (Exception e) {
         throw (IOException) new IOException ("Unable to create XML reader.").initCause (e);
      }
      xmlReader.setContentHandler (handler);
      try {
         xmlReader.parse (new InputSource (reader));
      } catch (SAXException e) {
         throw (IOException) new IOException ("Unable to parse input.").initCause (e);
      }
   }

   protected static String createAttributesText (Attributes attributes, Map namespacePrefixToURI) {
      int attributesCount = attributes.getLength ();
      StringBuffer text = new StringBuffer ("");
      if (attributesCount != 0)
         for (int attributesIndex = 0; attributesIndex < attributesCount; attributesIndex++)
            text.append (" " + attributes.getQName (attributesIndex) + "=\"" + attributes.getValue (attributesIndex) + "\"");
      if (namespacePrefixToURI != null)
         for (Iterator namespaceIterator = namespacePrefixToURI.entrySet ().iterator (); namespaceIterator.hasNext (); ) {
            Map.Entry entry = (Map.Entry) namespaceIterator.next ();
            String key = entry.getKey ().toString ().trim ();
            text.append ((key.length () == 0 ? " xmlns" : (" xmlns:" + key)) + "=\"" + entry.getValue () + "\"");
         }
      return text.toString ();
   }

   protected XMLReader () {
      errorStream = System.err;
      namespacePrefixToURI = new HashMap ();
      stack = new LinkedList ();
   }

   public void characters (char characters [], int start, int length) {
      if (localText == null)
         return;
      String text = new String (characters, start, length).trim ();
      if (text.length () != 0)
         localText.append (text);
   }

   public void startPrefixMapping (String prefix, String uri) {
      namespacePrefixToURI.put (prefix, uri);
   }

   protected String endLocalTextElement () {
      String text = localText.toString ();
      localText = null;
      mode = NORMAL;
      return text;
   }

   protected XMLElement peek () {
      return (XMLElement) stack.getLast ();
   }

   protected XMLElement peekPastContainer () {
      return (XMLElement) stack.get (stack.size () - 2);
   }

   protected XMLElement pop () {
      return (XMLElement) stack.removeLast ();
   }

   protected void push (XMLElement element) {
      stack.addLast (element);
   }
}
