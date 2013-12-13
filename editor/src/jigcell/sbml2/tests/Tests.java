package jigcell.sbml2.tests;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import junit.textui.TestRunner;

/**
 * Test harness for all SBML parser tests.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public abstract class Tests extends TestCase {

   /**
    * Starts a new test suite run.
    *
    * @param args Program arguments
    */

   public static void main (String args []) {
      TestRunner.run (suite ());
   }

   /**
    * A suite containing all test suites for the SBML parser.
    */

   public static Test suite () {
      TestSuite suite = new TestSuite ();
      suite.addTest (CompartmentModelTests.suite ());
      suite.addTest (CompartmentModelTests2.suite ());
      suite.addTest (AlgebraicModelTests.suite ());
      suite.addTest (AlgebraicModelTests2.suite ());
      suite.addTest (AssignmentModelTests.suite ());
      suite.addTest (AssignmentModelTests2.suite ());
      suite.addTest (BoundaryModelTests.suite ());
      suite.addTest (BoundaryModelTests2.suite ());
      suite.addTest (BranchModelTests.suite ());
      suite.addTest (BranchModelTests2.suite ());
      suite.addTest (DelayModelTests.suite ());
      suite.addTest (DelayModelTests2.suite ());
      suite.addTest (EventModelTests.suite ());
      suite.addTest (EventModelTests2.suite ());
      suite.addTest (FunctionModelTests.suite ());
      suite.addTest (FunctionModelTests2.suite ());
      suite.addTest (ODEModelTests.suite ());
      suite.addTest (ODEModelTests2.suite ());
      suite.addTest (UnitModelTests.suite ());
      suite.addTest (UnitModelTests2.suite ());
      return suite;
   }

   /**
    * Tests should never be instantiated.
    */

   private Tests () {}
}
