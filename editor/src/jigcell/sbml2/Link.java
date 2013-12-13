package jigcell.sbml2;

import org.xml.sax.Attributes;

/**
 * Links two components together for use in a composition hierarchy.  Link structures link together components within different models.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * <p>
 *
 * @author Ranjit Randhawa
 */

public final class Link extends SBaseId {

   private FromObjectReference fromObjectReference;
   private ToObjectReference toObjectReference;
   private final SBase fromObjectElement;
   private final SBase toObjectElement;

   public Link () {
      this (null, null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param link Link, must not be null
    */

   public Link (Link link) {
      this (link.getId (), link.isNameSet () ? link.getName () : null);
   }

   public Link (String id, String name) {
      super (id, name);
      fromObjectReference = new FromObjectReference ();
      toObjectReference = new ToObjectReference ();
      fromObjectElement = new SBase ();
      toObjectElement = new SBase ();
   }

   /**
    * From object methods
    */
   public void addFromObject (FromObjectReference reference) {
      if (reference == null)
         throw new IllegalArgumentException ();
      fromObjectReference = reference;         
   }
   
   public FromObjectReference getFromObjectReference () {
      return fromObjectReference;
   }

   public SBase getFromObjectElement () {
      return fromObjectElement;
   }

   /**
    * To object methods
    */
   public void addToObject (ToObjectReference reference) {
      if (reference == null)
         throw new IllegalArgumentException ();
      toObjectReference = reference;         
   }
   
   public ToObjectReference getToObjectReference () {
      return toObjectReference;
   }

   public SBase getToObjectElement () {
      return toObjectElement;
   }
   
   protected void parse (Attributes attributes) {
      super.parse (attributes);
    }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "link");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addText (getFromObjectReference ().toString ());
      printer.addText (getToObjectReference ().toString ());
      return printer;
   }
}
