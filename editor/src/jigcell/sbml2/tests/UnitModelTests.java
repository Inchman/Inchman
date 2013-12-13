package jigcell.sbml2.tests;

import java.util.List;
import jigcell.sbml2.Reaction;
import jigcell.sbml2.Unit;
import jigcell.sbml2.UnitDefinition;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.textui.TestRunner;

/**
 * A collection of tests to exercise the reference model for units.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class UnitModelTests extends SBMLModelTests {

   /**
    * Starts a new test suite run.
    *
    * @param args Program arguments
    */

   public static void main (String args []) {
      TestRunner.run (suite ());
   }

   /**
    * All of the tests in this test.
    */

   public static Test suite () {
      return new TestSuite (UnitModelTests.class);
   }

   /**
    * Creates a new test for the reference model for units.
    *
    * @param name Test name
    */

   public UnitModelTests (String name) {
      super (name);
      modelName = "l2v1-units.sbml";
   }

   /**
    * A test for compartments.
    */

   public void testCompartment0 () throws Exception {
      assertTrue (model.getCompartments ().size () == 1);
   }

   /**
    * A test for compartments.
    */

   public void testCompartment1 () throws Exception {
      testCompartment ("cell", true, null, Double.NaN, 3, "volume");
   }

   /**
    * A test for units.
    */

   public void testEvent0 () throws Exception {
      assertTrue (model.getEvents ().size () == 0);
   }

   /**
    * A test for model.
    */

   public void testModel1 () throws Exception {
      assertTrue (model.getId () == null);
   }

   /**
    * A test for parameters.
    */

   public void testParameter0 () throws Exception {
      assertTrue (model.getParameters ().size () == 2);
   }

   /**
    * A test for parameters.
    */

   public void testParameter1 () throws Exception {
      testParameter ("vm", 2.0, "mmls", true);
   }

   /**
    * A test for parameters.
    */

   public void testParameter2 () throws Exception {
      testParameter ("km", 2.0, "mm", true);
   }

   /**
    * A test for reactions.
    */

   public void testReaction0 () throws Exception {
      assertTrue (model.getReactions ().size () == 3);
   }

   /**
    * A test for reactions.
    */

   public void testReaction1 () throws Exception {
      Reaction reaction = testReaction ("v1", null, true);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "x0", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "s1", 1);
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 0);
   }

   /**
    * A test for reactions.
    */

   public void testReaction2 () throws Exception {
      Reaction reaction = testReaction ("v2", null, true);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "s1", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "s2", 1);
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 0);
   }

   /**
    * A test for reactions.
    */

   public void testReaction3 () throws Exception {
      Reaction reaction = testReaction ("v3", null, true);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "s2", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "x1", 1);
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 0);
   }

   /**
    * A test for rules.
    */

   public void testRule0 () throws Exception {
      assertTrue (model.getRules ().size () == 0);
   }

   /**
    * A test for species.
    */

   public void testSpecies0 () throws Exception {
      assertTrue (model.getSpecies ().size () == 4);
   }

   /**
    * A test for species.
    */

   public void testSpecies1 () throws Exception {
      testSpecies ("x0", "cell", 1.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies2 () throws Exception {
      testSpecies ("x1", "cell", 1.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies3 () throws Exception {
      testSpecies ("s1", "cell", 1.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies4 () throws Exception {
      testSpecies ("s2", "cell", 1.0, false, false);
   }

   /**
    * A test for units.
    */

   public void testUnit0 () throws Exception {
      assertTrue (model.getUnitDefinitions ().size () == 3);
   }

   /**
    * A test for units.
    */

   public void testUnit1 () throws Exception {
      UnitDefinition definition = testUnitDefinition ("substance", 1);
      List units = definition.getUnits ();
      testUnit ((Unit) units.get (0), 1, "mole", 1.0, 0.0, -3);
   }

   /**
    * A test for units.
    */

   public void testUnit2 () throws Exception {
      UnitDefinition definition = testUnitDefinition ("mmls", 3);
      List units = definition.getUnits ();
      testUnit ((Unit) units.get (0), 1, "mole", 1.0, 0.0, -3);
      testUnit ((Unit) units.get (1), -1, "litre", 1.0, 0.0, 0);
      testUnit ((Unit) units.get (2), -1, "second", 1.0, 0.0, 0);
   }

   /**
    * A test for units.
    */

   public void testUnit3 () throws Exception {
      UnitDefinition definition = testUnitDefinition ("mm", 1);
      List units = definition.getUnits ();
      testUnit ((Unit) units.get (0), 1, "mole", 1.0, 0.0, -3);
   }
}
