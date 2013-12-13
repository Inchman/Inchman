package jigcell.sbml2;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.io.IOException;
import java.io.Writer;

/**
 * Used to format and print XML.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public final class XMLPrinter {
   private final List elements;
   private final Map attributes;
   private final String indent;
   private final String name;
   private String text;

   public static String quote (String cdata) {
      if (cdata == null)
         return null;
      cdata = cdata.replaceAll ("&", "&amp;");
      cdata = cdata.replaceAll ("<", "&lt;");
      cdata = cdata.replaceAll (">", "&gt;");
      cdata = cdata.replaceAll ("\"", "&quot;");
      return cdata.replaceAll ("'", "&apos;");
   }

   public static String unquote (String cdata) {
      if (cdata == null)
         return null;
      cdata = cdata.replaceAll ("&amp;", "&");
      cdata = cdata.replaceAll ("&lt;", "<");
      cdata = cdata.replaceAll ("&gt;", ">");
      cdata = cdata.replaceAll ("&quot;", "\"");
      return cdata.replaceAll ("&apos;", "'");
   }

   public XMLPrinter (String name) {
      this (null, name);
   }

   public XMLPrinter (XMLPrinter parent, String name) {
      this.name = name;
      indent = parent == null ? "" : parent.indent + "  ";
      attributes = new LinkedHashMap ();
      elements = new ArrayList ();
   }

   public void addAttribute (String attribute, String value) {
      if (value != null && value.length () > 0)
         attributes.put (attribute, value);
   }

   public void addCustomElement (XMLElement element, String name, String cdata) {
      XMLPrinter elementPrinter = element.print (this, name);
      elementPrinter.addText (cdata);
      elements.add (elementPrinter);
   }

   public void addElement (XMLElement element) {
      if (element != null)
         elements.add (element.print (this));
   }

   public void addElementList (XMLElement element, String name, Collection collection) {
      if (collection.isEmpty ())
         return;
      XMLPrinter elementPrinter = element.print (this, name);
      for (Iterator iterator = collection.iterator (); iterator.hasNext (); )
         elementPrinter.addElement ((XMLElement) iterator.next ());
      elements.add (elementPrinter);
   }

   public void addText (String cdata) {
      if (cdata != null)
         text = text == null ? cdata : text + cdata;
   }

   public String toString () {
      return toString (new StringBuffer ());
   }

   public void write (Writer writer) throws IOException {
      writer.write (indent + "<" + name);
      for (Iterator iterator = attributes.entrySet ().iterator (); iterator.hasNext (); ) {
         Map.Entry entry = (Map.Entry) iterator.next ();
         writer.write (" " + (String) entry.getKey () + "=\"" + quote ((String) entry.getValue ()) + "\"");
      }
      if (text == null && elements.isEmpty ()) {
         writer.write (" />\n");
         return;
      }
      writer.write (">\n");
      if (text != null) {
         writer.write (text);
         if (!text.endsWith ("\n"))
            writer.write ("\n");
      }
      for (Iterator iterator = elements.iterator (); iterator.hasNext (); )
         ((XMLPrinter) iterator.next ()).write (writer);
      writer.write (indent + "</" + name + ">\n");
   }

   private String toString (StringBuffer xml) {
      xml.append (indent + "<" + name);
      for (Iterator iterator = attributes.entrySet ().iterator (); iterator.hasNext (); ) {
         Map.Entry entry = (Map.Entry) iterator.next ();
         xml.append (" " + (String) entry.getKey () + "=\"" + quote ((String) entry.getValue ()) + "\"");
      }
      if (text == null && elements.isEmpty ())
         return xml.append (" />\n").toString ();
      xml.append (">\n");
      if (text != null) {
         xml.append (text);
         if (!text.endsWith ("\n"))
            xml.append ("\n");
      }
      for (Iterator iterator = elements.iterator (); iterator.hasNext (); )
         ((XMLPrinter) iterator.next ()).toString (xml);
      return xml.append (indent + "</" + name + ">\n").toString ();
   }
}
