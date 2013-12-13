package jigcell.sbml2.tests;

import jigcell.sbml2.Reaction;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.textui.TestRunner;

/**
 * A collection of tests to exercise the reference model for function definitions.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class FunctionModelTests extends SBMLModelTests {

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
      return new TestSuite (FunctionModelTests.class);
   }

   /**
    * Creates a new test for the reference model for function definitions.
    *
    * @param name Test name
    */

   public FunctionModelTests (String name) {
      super (name);
      modelName = "l2v1-functions.sbml";
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
      testCompartment ("compartmentOne", true, null, 1.0, 3, "volume");
   }

   /**
    * A test for events.
    */

   public void testEvent0 () throws Exception {
      assertTrue (model.getEvents ().size () == 0);
   }

   /**
    * A test for function definitions.
    */

   public void testFunction0 () throws Exception {
      assertTrue (model.getFunctionDefinitions ().size () == 1);
   }

   /**
    * A test for function definitions.
    */

   public void testFunction1 () throws Exception {
      testFunctionDefinition ("f");
   }

   /**
    * A test for model.
    */

   public void testModel1 () throws Exception {
      assertTrue (model.getId ().equals ("Example"));
   }

   /**
    * A test for parameters.
    */

   public void testParameter0 () throws Exception {
      assertTrue (model.getParameters ().size () == 0);
   }

   /**
    * A test for reactions.
    */

   public void testReaction0 () throws Exception {
      assertTrue (model.getReactions ().size () == 1);
   }

   /**
    * A test for reactions.
    */

   public void testReaction1 () throws Exception {
      Reaction reaction = testReaction ("reaction_1", null, false);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "S1", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "S2", 1);
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
      assertTrue (model.getSpecies ().size () == 2);
   }

   /**
    * A test for species.
    */

   public void testSpecies1 () throws Exception {
      testSpecies ("S1", "compartmentOne", 0.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies2 () throws Exception {
      testSpecies ("S2", "compartmentOne", 0.0, false, false);
   }

   /**
    * A test for units.
    */

   public void testUnit0 () throws Exception {
      assertTrue (model.getUnitDefinitions ().size () == 0);
   }
}
