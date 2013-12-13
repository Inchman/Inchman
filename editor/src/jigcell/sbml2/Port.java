package jigcell.sbml2;

import org.xml.sax.Attributes;

/**
 * An interface/terminal to an species or parameter.  Ports can be linked together to connect components within a model
 * to components in other models.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Ranjit Randhawa
 */

public final class Port extends SBaseId {

   private TargetObjectReference targetObjectReference;
   private final SBase targetObjectElement;
   private boolean input;
   
   public Port () {
      this (null, null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param Port Port, must not be null
    */

   public Port (Port port) {
      this (port.getId (), port.isNameSet () ? port.getName () : null);
      setInput (port.isInput ());      
   }

   public Port (String id, String name) {
      super (id, name);
      setInput (false);      
      targetObjectReference = new TargetObjectReference ();      
      targetObjectElement = new SBase ();
   }

   public boolean isInput () {
      return input;
   }
   
   public void setInput (boolean input) {
      this.input = input;
   }
   
   /**
    * Target object methods
    */
   public void addTargetObject (TargetObjectReference reference) {
      if (reference == null)
         throw new IllegalArgumentException ();
      targetObjectReference = reference;         
   }

   public TargetObjectReference getTargetObjectReference () {
      return targetObjectReference;
   }

   public SBase getTargetObjectElement () {
      return targetObjectElement;
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      if (attributes.getIndex ("input") != -1)
         setInput (Boolean.valueOf (attributes.getValue ("input")) == Boolean.TRUE);
   }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "port");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      if (isInput ())
         printer.addAttribute ("input", "true");
      printer.addText (getTargetObjectReference ().toString ());      
      return printer;
   }
}
