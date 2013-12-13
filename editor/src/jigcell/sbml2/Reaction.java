package jigcell.sbml2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.xml.sax.Attributes;

/**
 * A transformation, transport, or binding process that changes the amount of one or more species.  A reaction is defined in terms of the
 * participating reactants, products, and modifiers along with a kinetic law that describes the rate at which the reaction takes place.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class Reaction extends SBaseId {
   private Boolean fast;
   private boolean reversible;
   private KineticLaw kineticLaw;
   protected List modifier;
   protected List product;
   protected List reactant;
   private final SBase modifiersElement;
   private final SBase productsElement;
   private final SBase reactantsElement;

   public Reaction () {
      this (null, null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param reaction Reaction, must not be null
    */

   public Reaction (Reaction reaction) {
      this (reaction.getId (), reaction.isNameSet () ? reaction.getName () : null);
      if (reaction.isFastSet ())
         setFast (reaction.isFast ());
      setReversible (reaction.isReversible ());
      if (reaction.getKineticLaw () != null)
         setKineticLaw (new KineticLaw (reaction.getKineticLaw ()));
      for (Iterator iterator = reaction.getModifier ().iterator (); iterator.hasNext (); ) {
         ModifierSpeciesReference reference = (ModifierSpeciesReference) iterator.next ();
         if (reference != null)
            addModifier (new ModifierSpeciesReference (reference));
      }
      for (Iterator iterator = reaction.getProduct ().iterator (); iterator.hasNext (); ) {
         SpeciesReference reference = (SpeciesReference) iterator.next ();
         if (reference != null)
            addProduct (new SpeciesReference (reference));
      }
      for (Iterator iterator = reaction.getReactant ().iterator (); iterator.hasNext (); ) {
         SpeciesReference reference = (SpeciesReference) iterator.next ();
         if (reference != null)
            addReactant (new SpeciesReference (reference));
      }
   }

   public Reaction (String id, String name) {
      super (id, name);
      setReversible (true);
      unsetFast ();
      modifier = new ArrayList ();
      modifiersElement = new SBase ();
      product = new ArrayList ();
      productsElement = new SBase ();
      reactant = new ArrayList ();
      reactantsElement = new SBase ();
   }

   public void addModifier (Species reference) {
      addModifier (new ModifierSpeciesReference (reference));
   }

   public void addModifier (ModifierSpeciesReference reference) {
      if (reference == null)
         throw new IllegalArgumentException ();
      modifier.add (reference);
   }

   public void addProduct (Species reference) {
      addProduct (new SpeciesReference (reference));
   }

   public void addProduct (SpeciesReference reference) {
      if (reference == null)
         throw new IllegalArgumentException ();
      product.add (reference);
   }

   public void addReactant (Species reference) {
      addReactant (new SpeciesReference (reference));
   }

   public void addReactant (SpeciesReference reference) {
      if (reference == null)
         throw new IllegalArgumentException ();
      reactant.add (reference);
   }

   public KineticLaw getKineticLaw () {
      return kineticLaw;
   }

   public List getModifier () {
      return modifier;
   }

   public SBase getModifiersElement () {
      return modifiersElement;
   }

   public List getProduct () {
      return product;
   }

   public SBase getProductsElement () {
      return productsElement;
   }

   public List getReactant () {
      return reactant;
   }

   public SBase getReactantsElement () {
      return reactantsElement;
   }

   public boolean isFast () {
      assert isFastSet ();
      return isFastSet () && fast.booleanValue ();
   }

   public boolean isFastSet () {
      return fast != null;
   }

   public boolean isReversible () {
      return reversible;
   }

   public void setFast (boolean fast) {
      this.fast = fast ? Boolean.TRUE : Boolean.FALSE;
   }

   public void setKineticLaw (KineticLaw kineticLaw) {
      this.kineticLaw = kineticLaw;
   }

   /**
    * Sets whether the reaction is reversible.
    */

   public void setReversible (boolean reversible) {
      this.reversible = reversible;
   }

   public void unsetFast () {
      fast = null;
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      if (attributes.getIndex ("fast") != -1)
         setFast (Boolean.valueOf (attributes.getValue ("fast")) == Boolean.TRUE);
      if (attributes.getIndex ("reversible") != -1)
         setReversible (Boolean.valueOf (attributes.getValue ("reversible")) == Boolean.TRUE);
    }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "reaction");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      if (!isReversible ())
         printer.addAttribute ("reversible", "false");
      if (isFastSet ())
         printer.addAttribute ("fast", isFast () ? "true" : "false");
      printer.addElementList (getReactantsElement (), "listOfReactants", getReactant ());
      printer.addElementList (getProductsElement (), "listOfProducts", getProduct ());
      printer.addElementList (getModifiersElement (), "listOfModifiers", getModifier ());
      printer.addElement (getKineticLaw ());
      return printer;
   }
}
