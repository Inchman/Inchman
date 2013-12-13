package jigcell.sbml2;

import org.xml.sax.Attributes;

/**
 * A reference to a species contained in the model.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public abstract class SimpleSpeciesReference extends SBase {
   private String species;

   public String getSpecies () {
      return species;
   }

   public final Species getSpecies (Model model) {
      return (Species) model.findElementWithId (getSpecies (), Model.SPECIES);
   }

   public boolean isValid (Model model) {
      return super.isValid (model) && getSpecies (model) != null;
   }

   public final void setSpecies (String species) {
      if (species != null && !SBaseId.isValidId (species))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.species = species;
   }

   public final void setSpecies (Species species) {
      setSpecies (species == null ? null : species.getId ());
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      setSpecies (attributes.getValue ("species"));
    }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addAttribute ("species", getSpecies ());
      return printer;
   }
}
