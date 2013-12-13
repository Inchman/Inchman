package jigcell.sbml2;

import org.xml.sax.Attributes;

/**
 * Instances indicate that a copy of a model is being instantiated within the current model.
 * Instance structures use XML Linking Language (XLink) to refer to local and remote submodels.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * <p>
 *
 * @author Ranjit Randhawa
 */

public final class Instance extends SBaseId {

   private String type;
   private String href;

   public Instance () {
      this (null, null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param instance Instance, must not be null
    */

   public Instance (Instance instance) {
      this (instance.getId (), instance.isNameSet () ? instance.getName () : null);
      setType (instance.getType ());
      setHref (instance.getHref ());
   }

   public Instance (Model model) {
      this (model.getId (), model.getName ());
   }

   public Instance (String id, String name) {
      super (id, name);
      setType (null);
   }

   public String getType () {
      return type;
   }

   public String getHref () {
      return href;
   }

   /**
    * Sets the instance's type, the default value is "simple" to indicate that only two resources
    * are associated to each other.
    *
    * @param type Type of xlink.
    */

   public void setType (String type) {
      if (type == null)
         type = "simple";
      this.type = type;
   }

   /**
    * Sets the instance's href. The href attribute contains a string that points to either an SBML
    * model document or a model element within the current SBML document.
    *
    * @param href contains an xpointer string.
    */
   public void setHref (String href) {
      this.href = href;
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      setType (attributes.getValue ("xlink:type"));
      setHref (attributes.getValue ("xlink:href"));
   }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "instance");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addAttribute ("xlink:type", getType ());
      printer.addAttribute ("xlink:href", getHref ());
      return printer;
   }
}
