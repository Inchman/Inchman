package jigcell.sbml2.tests;

import java.util.List;
import jigcell.sbml2.Reaction;
import jigcell.sbml2.Unit;
import jigcell.sbml2.UnitDefinition;
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

public class CompartmentModelTests extends SBMLModelTests {

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
      return new TestSuite (CompartmentModelTests.class);
   }

   /**
    * Creates a new test for the reference model for compartments.
    *
    * @param name Test name
    */

   public CompartmentModelTests (String name) {
      super (name);
      modelName = "l2v1-2D-compartments.sbml";
   }

   /**
    * A test for compartments.
    */

   public void testCompartment0 () throws Exception {
      assertTrue (model.getCompartments ().size () == 3);
   }

   /**
    * A test for compartments.
    */

   public void testCompartment1 () throws Exception {
      testCompartment ("Extracellular", true, null, Double.NaN, 3, "volume");
   }

   /**
    * A test for compartments.
    */

   public void testCompartment2 () throws Exception {
      testCompartment ("PlasmaMembrane", true, "Extracellular", Double.NaN, 2, "area");
   }

   /**
    * A test for compartments.
    */

   public void testCompartment3 () throws Exception {
      testCompartment ("Cytosol", true, "PlasmaMembrane", Double.NaN, 3, "volume");
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
      assertTrue (model.getId ().equals ("facilitated_ca_diffusion"));
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
      assertTrue (model.getReactions ().size () == 4);
   }

   /**
    * A test for reactions.
    */

   public void testReaction1 () throws Exception {
      Reaction reaction = testReaction ("CalciumCalbindin_gt_BoundCytosol", Boolean.TRUE, true);
      assertTrue (reaction.getReactant ().size () == 2);
      testReactionModifier (reaction.getReactant (), "CaBP_C", 1);
      testReactionModifier (reaction.getReactant (), "Ca_C", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "CaBPB_C", 1);
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 2);
      testParameter (reaction, "Kf_CalciumCalbindin_BoundCytosol", 20.0);
      testParameter (reaction, "Kr_CalciumCalbindin_BoundCytosol", 8.6);
   }

   /**
    * A test for reactions.
    */

   public void testReaction2 () throws Exception {
      Reaction reaction = testReaction ("CalciumBuffer_gt_BoundCytosol", Boolean.TRUE, true);
      assertTrue (reaction.getReactant ().size () == 2);
      testReactionModifier (reaction.getReactant (), "Ca_C", 1);
      testReactionModifier (reaction.getReactant (), "B_C", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "CaB_C", 1);
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 2);
      testParameter (reaction, "Kf_CalciumBuffer_BoundCytosol", 0.1);
      testParameter (reaction, "Kr_CalciumBuffer_BoundCytosol", 1.0);
   }

   /**
    * A test for reactions.
    */

   public void testReaction3 () throws Exception {
      Reaction reaction = testReaction ("Ca_Pump", null, true);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "Ca_C", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "Ca_EC", 1);
      assertTrue (reaction.getModifier ().size () == 1);
      testReactionModifier (reaction.getModifier (), "CaPump_PM");
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 3);
      testParameter (reaction, "Vmax", -4000.0);
      testParameter (reaction, "kP", 0.25);
      testParameter (reaction, "Ca_Rest", 0.1);
   }

   /**
    * A test for reactions.
    */

   public void testReaction4 () throws Exception {
      Reaction reaction = testReaction ("Ca_channel", null, true);
      assertTrue (reaction.getReactant ().size () == 1);
      testReactionModifier (reaction.getReactant (), "Ca_EC", 1);
      assertTrue (reaction.getProduct ().size () == 1);
      testReactionModifier (reaction.getProduct (), "Ca_C", 1);
      assertTrue (reaction.getModifier ().size () == 1);
      testReactionModifier (reaction.getModifier (), "CaCh_PM");
      assertTrue (reaction.getKineticLaw ().getParameter ().size () == 2);
      testParameter (reaction, "J0", 0.014);
      testParameter (reaction, "Kc", 0.5);
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
      assertTrue (model.getSpecies ().size () == 8);
   }

   /**
    * A test for species.
    */

   public void testSpecies1 () throws Exception {
      testSpecies ("CaBPB_C", "Cytosol", 47.17, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies2 () throws Exception {
      testSpecies ("B_C", "Cytosol", 396.04, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies3 () throws Exception {
      testSpecies ("CaB_C", "Cytosol", 3.96, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies4 () throws Exception {
      testSpecies ("Ca_EC", "Extracellular", 1000.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies5 () throws Exception {
      testSpecies ("Ca_C", "Cytosol", 0.1, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies6 () throws Exception {
      testSpecies ("CaCh_PM", "PlasmaMembrane", 1.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies7 () throws Exception {
      testSpecies ("CaPump_PM", "PlasmaMembrane", 1.0, false, false);
   }

   /**
    * A test for species.
    */

   public void testSpecies8 () throws Exception {
      testSpecies ("CaBP_C", "Cytosol", 202.83, false, false);
   }

   /**
    * A test for units.
    */

   public void testUnit0 () throws Exception {
      assertTrue (model.getUnitDefinitions ().size () == 2);
   }

   /**
    * A test for units.
    */

   public void testUnit1 () throws Exception {
      UnitDefinition definition = testUnitDefinition ("substance", 1);
      List units = definition.getUnits ();
      testUnit ((Unit) units.get (0), 1, "mole", 1.0, 0.0, -6);
   }

   /**
    * A test for units.
    */

   public void testUnit2 () throws Exception {
      UnitDefinition definition = testUnitDefinition ("area", 1);
      List units = definition.getUnits ();
      testUnit ((Unit) units.get (0), 2, "metre", 1.0, 0.0, -6);
   }
}
