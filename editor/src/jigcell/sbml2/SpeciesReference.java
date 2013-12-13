package jigcell.sbml2;

import org.xml.sax.Attributes;

/**
 * A reference to a species contained in the model.  This type of reference is used for the products and reactants of a reaction and contains a
 * stoichiometry that describes how many copies of the species are involved in the reaction.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public final class SpeciesReference extends SimpleSpeciesReference {
   private double stoichiometry;
   private StoichiometryMath stoichiometryMath;

   public SpeciesReference () {
      this ((Species) null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param reference Reference, must not be null
    */

   public SpeciesReference (SpeciesReference reference) {
      this ();
      setSpecies (reference.getSpecies ());
      setStoichiometry (reference.getStoichiometry ());
      if (reference.getStoichiometryMath () != null)
         setStoichiometryMath (new StoichiometryMath (reference.getStoichiometryMath ()));
   }

   public SpeciesReference (Species species) {
      setSpecies (species);
      setStoichiometry (1.0);
   }

   public double getStoichiometry () {
      return stoichiometry;
   }

   public StoichiometryMath getStoichiometryMath () {
      return stoichiometryMath;
   }

   /**
    * Sets the rational stoichiometry of a reaction.  NaN indicates that the stoichiometry is not set.  If the stoichiometry is not NaN, the
    * stoichiometry math will be cleared.  Stoichiometries must be non-negative.  For maximum compatibility with other software, the
    * stoichiometry should be an integer.
    *
    * @param stoichiometry Stoichiometry
    */

   public void setStoichiometry (double stoichiometry) {
      if (stoichiometry < 0.0)
         throw new IllegalArgumentException ("Stoichiometry must be non-negative.");
      this.stoichiometry = stoichiometry;
      if (!Double.isNaN (stoichiometry))
         setStoichiometryMath (null);
   }

   /**
    * Sets the stoichiometry of a reaction.  If the stoichiometry math is not null, the rational stoichiometry will be cleared.  For maximum
    * compatibility with other software, use a rational stoichiometry instead of stoichiometry math.
    *
    * @param stoichiometryMath Stoichiometry math
    */

   public void setStoichiometryMath (StoichiometryMath stoichiometryMath) {
      this.stoichiometryMath = stoichiometryMath;
      if (stoichiometryMath != null)
         stoichiometry = Double.NaN;
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      if (attributes.getIndex ("stoichiometry") != -1)
         setStoichiometry (Double.parseDouble (attributes.getValue ("stoichiometry")));
     }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "speciesReference");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addElement (getStoichiometryMath ());
      if (!Double.isNaN (getStoichiometry ()) && getStoichiometry () != 1.0)
         printer.addAttribute ("stoichiometry", String.valueOf (getStoichiometry ()));
      return printer;
   }
}
