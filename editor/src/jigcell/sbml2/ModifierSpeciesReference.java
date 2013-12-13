package jigcell.sbml2;

/**
 * A reference to a species contained in the model.  This type of reference is used for the modifiers of a reaction.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public final class ModifierSpeciesReference extends SimpleSpeciesReference {
   public ModifierSpeciesReference () {
      this ((Species) null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param reference Reference, must not be null
    */

   public ModifierSpeciesReference (ModifierSpeciesReference reference) {
      this ();
      setSpecies (reference.getSpecies ());
   }

   public ModifierSpeciesReference (Species species) {
      setSpecies (species);
   }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "modifierSpeciesReference");
   }
}
