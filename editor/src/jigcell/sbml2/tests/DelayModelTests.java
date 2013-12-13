package jigcell.sbml2.tests;

import junit.framework.Test;
import junit.framework.TestSuite;
import junit.textui.TestRunner;

/**
 * A collection of tests to exercise the reference model for delay equations.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class DelayModelTests extends SBMLModelTests {

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
      return new TestSuite (DelayModelTests.class);
   }

   /**
    * Creates a new test for the reference model for delay equations.
    *
    * @param name Test name
    */

   public DelayModelTests (String name) {
      super (name);
      modelName = "l2v1-delay.sbml";
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
      testCompartment ("cell", true, null, 1.0, 3, "volume");
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
      assertTrue (model.getParameters ().size () == 4);
   }

   /**
    * A test for parameters.
    */

   public void testParameter1 () throws Exception {
      testParameter ("tau", 1.0);
   }

   /**
    * A test for parameters.
    */

   public void testParameter2 () throws Exception {
      testParameter ("m", 0.5);
   }

   /**
    * A test for parameters.
    */

   public void testParameter3 () throws Exception {
      testParameter ("q", 1.0);
   }

   /**
    * A test for parameters.
    */

   public void testParameter4 () throws Exception {
      testParameter ("delta_t", 1.0);
   }

   /**
    * A test for reactions.
    */

   public void testReaction0 () throws Exception {
      assertTrue (model.getReactions ().size () == 0);
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
      assertTrue (model.getSpecies ().size () == 1);
   }

   /**
    * A test for species.
    */

   public void testSpecies1 () throws Exception {
      testSpecies ("P", "cell", 0.0, false, false);
   }

   /**
    * A test for units.
    */

   public void testUnit0 () throws Exception {
      assertTrue (model.getUnitDefinitions ().size () == 0);
   }
}
