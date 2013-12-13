package jigcell.sbml2.tests;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Iterator;
import java.util.List;
import jigcell.sbml2.Compartment;
import jigcell.sbml2.FunctionDefinition;
import jigcell.sbml2.Model;
import jigcell.sbml2.ModifierSpeciesReference;
import jigcell.sbml2.Parameter;
import jigcell.sbml2.Reaction;
import jigcell.sbml2.SBMLLevel2Document;
import jigcell.sbml2.Species;
import jigcell.sbml2.SpeciesReference;
import jigcell.sbml2.Unit;
import jigcell.sbml2.UnitDefinition;
import junit.framework.TestCase;

/**
 * A collection of test methods for pieces of an SBML model.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

abstract class SBMLModelTests extends TestCase {

   protected boolean cycle;
   protected Model model;
   protected String modelName;

   /**
    * Creates a new test for an SBML model.
    *
    * @param name Test name
    */

   public SBMLModelTests (String name) {
      super (name);
      cycle = false;
      modelName = null;
   }

   /*
   protected void setUp () throws Exception {
      SBMLLevel2Document document = SBMLLevel2Document.readDocument (SBMLModelTests.class.getResourceAsStream (modelName));
      model = document.getModel ();
      if (!cycle)
         return;
      modelName = "x-" + modelName;
      document.writeDocument (new BufferedWriter (new FileWriter (modelName)));
      model = SBMLLevel2Document.readDocument (modelName).getModel ();
   }
   */

   protected Compartment testCompartment (String name, boolean constant, String outside, double size, int dimensions, String units) {
      for (Iterator iterator = model.getCompartments ().iterator (); iterator.hasNext (); ) {
         Compartment compartment = (Compartment) iterator.next ();
         if (!compartment.getName ().equals (name))
            continue;
         assertTrue (compartment.isConstant () == constant);
         assertTrue (stringEquals (compartment.getOutside (), outside));
         assertTrue (size == compartment.getSize () || Double.isNaN (size) && Double.isNaN (compartment.getSize ()));
         assertTrue (dimensions == compartment.getSpatialDimensions ());
         assertTrue (stringEquals (compartment.getUnits (), units));
         return compartment;
      }
      fail ();
      return null;
   }

   protected FunctionDefinition testFunctionDefinition (String name) {
      for (Iterator iterator = model.getFunctionDefinitions ().iterator (); iterator.hasNext (); ) {
         FunctionDefinition definition = (FunctionDefinition) iterator.next ();
         if (!definition.getName ().equals (name))
            continue;
         return definition;
      }
      fail ();
      return null;
   }

   protected Parameter testParameter (String name, double value) {
      return testParameter (model.getParameters ().iterator (), name, value, null, true);
   }

   protected Parameter testParameter (String name, double value, String unit, boolean constant) {
      return testParameter (model.getParameters ().iterator (), name, value, unit, constant);
   }

   protected Parameter testParameter (Reaction reaction, String name, double value) {
      return testParameter (reaction.getKineticLaw ().getParameter ().iterator (), name, value, null, true);
   }

   protected Reaction testReaction (String name, Boolean fast, boolean reversible) {
      for (Iterator iterator = model.getReactions ().iterator (); iterator.hasNext (); ) {
         Reaction reaction = (Reaction) iterator.next ();
         if (!reaction.getName ().equals (name))
            continue;
         if (fast == null)
            assertTrue (!reaction.isFastSet ());
         else {
            assertTrue (reaction.isFastSet ());
            assertTrue (reaction.isFast () == fast.booleanValue ());
         }
         assertTrue (reaction.isReversible () == reversible);
         return reaction;
      }
      fail ();
      return null;
   }

   protected ModifierSpeciesReference testReactionModifier (List modifiers, String name) {
      for (Iterator iterator = modifiers.iterator (); iterator.hasNext (); ) {
         ModifierSpeciesReference reference = (ModifierSpeciesReference) iterator.next ();
         if (!reference.getSpecies ().equals (name))
            continue;
         return reference;
      }
      fail ();
      return null;
   }

   protected SpeciesReference testReactionModifier (List modifiers, String name, double stoichiometry) {
      for (Iterator iterator = modifiers.iterator (); iterator.hasNext (); ) {
         SpeciesReference reference = (SpeciesReference) iterator.next ();
         if (!reference.getSpecies ().equals (name))
            continue;
         assertTrue (stoichiometry == reference.getStoichiometry ());
         return reference;
      }
      fail ();
      return null;
   }

   protected Species testSpecies (String name, String compartment, double initialConcentration, boolean constant, boolean boundaryCondition) {
      for (Iterator iterator = model.getSpecies ().iterator (); iterator.hasNext (); ) {
         Species species = (Species) iterator.next ();
         if (!species.getName ().equals (name))
            continue;
         assertTrue (stringEquals (species.getCompartment (), compartment));
         assertTrue (doubleEquals (species.getInitialConcentration (), initialConcentration));
         assertTrue (species.isConstant () == constant);
         assertTrue (species.isBoundaryCondition () == boundaryCondition);
         return species;
      }
      fail ();
      return null;
   }

   protected void testUnit (Unit unit, int exponent, String kindName, double multiplier, double offset, int scale) {
      assertTrue (unit.getExponent () == exponent);
      assertTrue (unit.getKindName ().equals (kindName));
      assertTrue (unit.getMultiplier () == multiplier);
      assertTrue (unit.getOffset () == offset);
      assertTrue (unit.getScale () == scale);
   }

   protected UnitDefinition testUnitDefinition (String name, int size) {
      for (Iterator iterator = model.getUnitDefinitions ().iterator (); iterator.hasNext (); ) {
         UnitDefinition definition = (UnitDefinition) iterator.next ();
         if (definition.getName ().equals (name)) {
            assertTrue (definition.getUnits ().size () == size);
            return definition;
         }
      }
      fail ();
      return null;
   }

   private boolean doubleEquals (double d1, double d2) {
      return d1 == d2 || Double.isNaN (d1) && Double.isNaN (d2);
   }

   private boolean stringEquals (String s1, String s2) {
      return s1 == s2 || s1.equals (s2);
   }

   private Parameter testParameter (Iterator iterator, String name, double value, String unit, boolean constant) {
      while (iterator.hasNext ()) {
         Parameter parameter = (Parameter) iterator.next ();
         if (!parameter.getName ().equals (name))
            continue;
         assertTrue (doubleEquals (parameter.getValue (), value));
         assertTrue (stringEquals (parameter.getUnits (), unit));
         assertTrue (parameter.isConstant () == constant);
         return parameter;
      }
      fail ();
      return null;
   }
}
