package jigcell.sbml2.tests;

import junit.framework.Test;
import junit.framework.TestSuite;
import junit.textui.TestRunner;

/**
 * A collection of tests to exercise the reference model for compartments.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class CompartmentModelTests2 extends CompartmentModelTests {

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
      return new TestSuite (CompartmentModelTests2.class);
   }

   /**
    * Creates a new test for the reference model for compartments.
    *
    * @param name Test name
    */

   public CompartmentModelTests2 (String name) {
      super (name);
      cycle = true;
   }
}
