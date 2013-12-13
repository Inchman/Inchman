package jigcell.sbml2;

import org.xml.sax.Attributes;

/**
 * A simple reference to an object contained in the model.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Ranjit Randhawa
 */

public abstract class SimpleObjectReference extends SBase {
   private String object;

   public String getObject () {
      return object;
   }

   public boolean isValid (Model model) {
      return super.isValid (model);
   }

   public final void setObject (String object) {
      if (object != null && !SBaseId.isValidId (object))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.object = object;
   }
   
   public final void setObject (SimpleObjectReference object) {
      if (object != null && !SBaseId.isValidId (object.getObject ()))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      if (object != null)
         this.object = object.getObject ();
   }
   
   protected void parse (Attributes attributes) {
      super.parse (attributes);
      setObject (attributes.getValue ("object"));
    }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addAttribute ("object", getObject ());
      return printer;
   }
}
