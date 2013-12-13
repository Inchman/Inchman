package jigcell.sbml2.tests;

import jigcell.sbml2.Reaction;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.textui.TestRunner;

/**
 * A collection of tests to exercise the reference model for odes.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class ODEModelTests extends SBMLModelTests {

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
      return new TestSuite (ODEModelTests.class);
   }

   /**
    * Creates a new test for the reference model for odes.
    *
    * @param name Test name
    */

   public ODEModelTests (String name) {
      super (name);
      modelName = "l2v1-mc-ode.sbml";
   }

   /**
    * A test for compartments.
    */

   public void testCompartment0 () throws Exception {
      assertTrue (model.getCompartments ().size () == 2);
   }

   /**
    * A test for compartments.
    */

   public void testCompartment1 () throws Exception {
      testCompartment ("V0", true, null, 10.0, 3, "volume");
   }

   /**
    * A test for compartments.
    */

   public void testCompartment2 () throws Exception {
      testCompartment ("V1", false, null, 1.0, 3, "volume");
   }

   /**
    * A test for odes.
    */

   public void testEvent0 () throws Exception {
      assertTrue (model.getEvents ().size () == 0);
   }

   /**
    * A test for model.
    */

   public void testModel1 () throws Exception {
      assertTrue (model.getId ().equals ("ODEExampleModel"));
   }

   /**
    * A test for parameters.
    */

   public void testParameter0 () throws Exception {
      assertTrue (model.getParameters ().size () == 6);
   }

   /**
    * A test for parameters.
    */

   public void testParameter1 () throws Exception {
      testParameter ("K0", 0.1);
   }

   /**
    * A test for parameters.
    */

   public void testParameter2 () throws Exception {
      testParameter ("K1", 0.5);
   }

   /**
    * A test for parameters.
    */

   public void testParameter3 () throws Exception {
      testParameter ("K2", 0.1);
   }

   /**
    * A test for parameters.
    */

   public void testParameter4 () throws Exception {
      testParameter ("K3", 0.5);
   }

   /**
    * A test for parameters.
    */

   public void testParameter5 () throws Exception {
      testParameter ("Kv", 0.5);
   }

   /**
    * A test for parameters.
    */

   public void testParameter6 () throws Exception {
      testParameter ("Kin", 0.1);
   }

   /**
    * A test for reactions.
    */

   public void testReaction0 () throws Exception {
      assertTrue (model.getReactions ().size () == 4);
   }

   /**
    * A test for reactions.
    */

   public void testReaction1 () throws Exception {
      Reaction reaction = testReaction ("reaction_1", null, false);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "X0", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "S1", 1);
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 0);
   }

   /**
    * A test for reactions.
    */

   public void testReaction2 () throws Exception {
      Reaction reaction = testReaction ("reaction_2", null, false);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "S1", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "S2", 1);
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 0);
   }

   /**
    * A test for reactions.
    */

   public void testReaction3 () throws Exception {
      Reaction reaction = testReaction ("reaction_3", null, false);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "S2", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "S3", 1);
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 0);
   }

   /**
    * A test for reactions.
    */

   public void testReaction4 () throws Exception {
      Reaction reaction = testReaction ("reaction_4", null, false);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "S3", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "X4", 1);
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 0);
   }

   /**
    * A test for rules.
    */

   public void testRule0 () throws Exception {
      assertTrue (model.getRules ().size () == 2);
   }

   /**
    * A test for species.
    */

   public void testSpecies0 () throws Exception {
      assertTrue (model.getSpecies ().size () == 5);
   }

   /**
    * A test for species.
    */

   public void testSpecies1 () throws Exception {
      testSpecies ("X0", "V0", 0.0, false, true);
   }

   /**
    * A test for species.
    */

   public void testSpecies2 () throws Exception {
      testSpecies ("S1", "V0", 0.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies3 () throws Exception {
      testSpecies ("S2", "V0", 0.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies4 () throws Exception {
      testSpecies ("S3", "V1", 0.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies5 () throws Exception {
      testSpecies ("X4", "V1", 0.0, true, true);
   }

   /**
    * A test for units.
    */

   public void testUnit0 () throws Exception {
      assertTrue (model.getUnitDefinitions ().size () == 0);
   }
}
