package jigcell.sbml2;

import org.xml.sax.Attributes;
import java.io.StringWriter;

/**
 * Base class for elements used with the XML reader and printer.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public abstract class XMLElement {

   /**
    * Prints this element as XML.
    */

   public final String toString () {
      return print (null).toString ();
   }

   protected XMLElement () {}

   protected void parse (Attributes attributes) {}

   protected XMLPrinter print (XMLPrinter parent) {
      throw new UnsupportedOperationException ();
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      return new XMLPrinter (parent, name);
   }

   public String toXml() {
      XMLPrinter printer = print(null,"cool");
      StringWriter sw = new StringWriter();
      try {
         printer.write(sw);
      } catch (Exception e) {
         e.printStackTrace();
      }
      return sw.toString();
   }
}
