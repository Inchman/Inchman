package jigcell.sbml2.tests;

import jigcell.sbml2.Reaction;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.textui.TestRunner;

/**
 * A collection of tests to exercise the reference model for boundary conditions.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class BoundaryModelTests extends SBMLModelTests {

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
      return new TestSuite (BoundaryModelTests.class);
   }

   /**
    * Creates a new test for the reference model for boundary conditions.
    *
    * @param name Test name
    */

   public BoundaryModelTests (String name) {
      super (name);
      modelName = "l2v1-boundary.sbml";
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
      assertTrue (model.getFunctionDefinitions ().size () == 0);
   }

   /**
    * A test for model.
    */

   public void testModel1 () throws Exception {
      assertTrue (model.getId ().equals ("BoundaryCondExampleModel"));
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
      testParameter ("k1", 0.5);
   }

   /**
    * A test for parameters.
    */

   public void testParameter2 () throws Exception {
      testParameter ("k2", 0.1);
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
      assertTrue (reaction.getReactant ().size () == 2);
      testReactionModifier (reaction.getReactant (), "S1", 1);
      testReactionModifier (reaction.getReactant (), "S2", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "S4", 1);
      assertTrue (reaction.getModifier ().size () == 1);
      testReactionModifier (reaction.getModifier (), "S3");
   }

   /**
    * A test for rules.
    */

   public void testRule0 () throws Exception {
      assertTrue (model.getRules ().size () == 1);
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
      testSpecies ("S1", "compartmentOne", 0.0, false, true);
   }

   /**
    * A test for species.
    */

   public void testSpecies2 () throws Exception {
      testSpecies ("S2", "compartmentOne", 1.0, true, true);
   }

   /**
    * A test for species.
    */

   public void testSpecies3 () throws Exception {
      testSpecies ("S3", "compartmentOne", 3.0, true, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies4 () throws Exception {
      testSpecies ("S4", "compartmentOne", 0.0, false, false);
   }

   /**
    * A test for units.
    */

   public void testUnit0 () throws Exception {
      assertTrue (model.getUnitDefinitions ().size () == 0);
   }
}
