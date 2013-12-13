package jigcell.sbml2.tests;

import jigcell.sbml2.Reaction;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.textui.TestRunner;

/**
 * A collection of tests to exercise the reference model for algebraic rules.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class AlgebraicModelTests extends SBMLModelTests {

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
      return new TestSuite (AlgebraicModelTests.class);
   }

   /**
    * Creates a new test for the reference model for algebraic rules.
    *
    * @param name Test name
    */

   public AlgebraicModelTests (String name) {
      super (name);
      modelName = "l2v1-algebraic.sbml";
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
      assertTrue (model.getId () == null);
   }

   /**
    * A test for parameters.
    */

   public void testParameter0 () throws Exception {
      assertTrue (model.getParameters ().size () == 1);
   }

   /**
    * A test for parameters.
    */

   public void testParameter1 () throws Exception {
      testParameter ("Keq", 2.5);
   }

   /**
    * A test for reactions.
    */

   public void testReaction0 () throws Exception {
      assertTrue (model.getReactions ().size () == 2);
   }

   /**
    * A test for reactions.
    */

   public void testReaction1 () throws Exception {
      Reaction reaction = testReaction ("in", null, true);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "X0", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "T", 1);
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 1);
      testParameter (reaction, "k1", 0.1);
   }

   /**
    * A test for reactions.
    */

   public void testReaction2 () throws Exception {
      Reaction reaction = testReaction ("out", null, true);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "T", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "X1", 1);
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 1);
      testParameter (reaction, "k2", 0.15);
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
      testSpecies ("X0", "cell", 1.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies2 () throws Exception {
      testSpecies ("X1", "cell", 0.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies3 () throws Exception {
      testSpecies ("T", "cell", 0.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies4 () throws Exception {
      testSpecies ("S1", "cell", 0.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies5 () throws Exception {
      testSpecies ("S2", "cell", 0.0, false, false);
   }

   /**
    * A test for units.
    */

   public void testUnit0 () throws Exception {
      assertTrue (model.getUnitDefinitions ().size () == 0);
   }
}
