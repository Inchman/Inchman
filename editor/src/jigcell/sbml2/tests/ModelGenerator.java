package jigcell.sbml2.tests;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import jigcell.sbml2.Compartment;
import jigcell.sbml2.FunctionDefinition;
import jigcell.sbml2.KineticLaw;
import jigcell.sbml2.Model;
import jigcell.sbml2.Parameter;
import jigcell.sbml2.Reaction;
import jigcell.sbml2.SBMLLevel2Document;
import jigcell.sbml2.Species;

/**
 * Automated generators for creating random SBML models.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public final class ModelGenerator {

   /**
    * Generates a model with equal numbers of species and reactions.
    *
    * @param size Number of reactions to create
    */

   public static SBMLLevel2Document generate1 (int size) {
      Model model = new Model ("Test Model " + size);
      model.getNotes ().add ("Model with " + size + " reactions");
      FunctionDefinition massAction = new FunctionDefinition ("Function", "Function");
      massAction.setMath ("<math:math><math:lambda><math:bvar><math:ci>K</math:ci></math:bvar><math:bvar><math:ci>S</math:ci>" +
         "</math:bvar><math:apply><math:times/><math:ci>K</math:ci><math:ci>S</math:ci></math:apply></math:lambda></math:math>");
      model.addFunctionDefinition (massAction);
      Compartment compartment = new Compartment ("Compartment", "Compartment");
      compartment.getNotes ().add ("Everything");
      compartment.setSize (1);
      model.addCompartment (compartment);
      Species species = new Species ("Source", "Source");
      species.getNotes ().add (species.getName ());
      species.setCompartment (compartment);
      species.setInitialAmount (1.0);
      model.addSpecies (species);
      for (int i = 0; i < size; i++) {
         String reactionName = "Reaction" + i;
         Reaction reaction = new Reaction (reactionName, reactionName);
         reaction.setReversible (false);
         reaction.getNotes ().add (reactionName);
         reaction.addReactant (species);
         String parameterName = "Parameter" + i;
         Parameter parameter = new Parameter (parameterName, parameterName, 1.0 / (i + 2.0));
         parameter.getNotes ().add (parameterName);
         KineticLaw kineticLaw = new KineticLaw ("<math:math><math:apply><math:ci>" + massAction.getId () + "</math:ci><math:ci>" +
            parameter.getId () + "</math:ci><math:ci>" + species.getId () + "</math:ci></math:apply></math:math>");
         model.addParameter (parameter);
         reaction.setKineticLaw (kineticLaw);
         String speciesName = "Species" + i;
         species = new Species (speciesName, speciesName);
         species.getNotes ().add (speciesName);
         species.setCompartment (compartment);
         species.setInitialAmount (0.1);
         reaction.addProduct (species);
         model.addSpecies (species);
         model.addReaction (reaction);
      }
      return new SBMLLevel2Document (model);
   }

   /**
    * Generates a model with many more species than reactions.
    *
    * @param size Number of reactions to create
    */

   public static SBMLLevel2Document generate2 (int size) {
      Model model = new Model ("Test Model " + size);
      model.getNotes ().add ("Model with " + size + " reactions");
      FunctionDefinition function = new FunctionDefinition ("FourFunction", "FourFunction");
      function.setMath ("<math:math><math:lambda><math:bvar><math:ci>K</math:ci></math:bvar><math:bvar><math:ci>S1</math:ci></math:bvar>" +
         "<math:bvar><math:ci>S2</math:ci></math:bvar><math:bvar><math:ci>S3</math:ci></math:bvar><math:bvar><math:ci>S4</math:ci>" +
         "</math:bvar><math:apply><math:times/><math:ci>K</math:ci><math:apply><math:plus/><math:apply><math:plus/><math:apply><math:plus/>" +
         "<math:ci>S1</math:ci><math:ci>S2</math:ci></math:apply><math:ci>S3</math:ci></math:apply><math:ci>S4</math:ci></math:apply>" +
         "</math:apply></math:lambda></math:math>");
      model.addFunctionDefinition (function);
      Compartment compartment = new Compartment ("Compartment", "Compartment");
      compartment.getNotes ().add ("Everything");
      compartment.setSize (1);
      model.addCompartment (compartment);
      Species species = new Species ("Source", "Source");
      species.getNotes ().add (species.getName ());
      species.setCompartment (compartment);
      species.setInitialAmount (1.0);
      model.addSpecies (species);
      for (int i = 0; i < size; i++) {
         String reactionName = "Reaction" + i;
         Reaction reaction = new Reaction (reactionName, reactionName);
         reaction.setReversible (false);
         reaction.getNotes ().add (reactionName);
         reaction.addReactant (species);
         String parameterName = "Parameter" + i;
         Parameter parameter = new Parameter (parameterName, parameterName, 1.0 / (i + 2.0));
         parameter.getNotes ().add (parameterName);
         String speciesName2 = "Species" + (i * 4);
         Species species2 = new Species (speciesName2, speciesName2);
         species2.getNotes ().add (speciesName2);
         species2.setCompartment (compartment);
         species2.setInitialAmount (0.9);
         model.addSpecies (species2);
         String speciesName3 = "Species" + (i * 4 + 1);
         Species species3 = new Species (speciesName3, speciesName3);
         species3.getNotes ().add (speciesName3);
         species3.setCompartment (compartment);
         species3.setInitialAmount (0.8);
         model.addSpecies (species3);
         String speciesName4 = "Species" + (i * 4 + 2);
         Species species4 = new Species (speciesName4, speciesName4);
         species4.getNotes ().add (speciesName4);
         species4.setCompartment (compartment);
         species4.setInitialAmount (0.7);
         model.addSpecies (species4);
         KineticLaw kineticLaw = new KineticLaw ("<math:math><math:apply><math:ci>" + function.getId () + "</math:ci><math:ci>" +
            parameter.getId () + "</math:ci><math:ci>" + species.getId () + "</math:ci><math:ci>" + species2.getId () + "</math:ci><math:ci>" +
            species3.getId () + "</math:ci><math:ci>" + species4.getId () + "</math:ci></math:apply></math:math>");
         model.addParameter (parameter);
         reaction.setKineticLaw (kineticLaw);
         String speciesName = "Species" + (i * 4 + 3);
         species = new Species (speciesName, speciesName);
         species.getNotes ().add (speciesName);
         species.setCompartment (compartment);
         species.setInitialAmount (0.1);
         reaction.addProduct (species);
         model.addSpecies (species);
         model.addReaction (reaction);
      }
      return new SBMLLevel2Document (model);
   }

   /**
    * Generates a model with many more reactions than species.
    *
    * @param size Number of reactions to create
    */

   public static SBMLLevel2Document generate3 (int size) {
      Model model = new Model ("Test Model " + size);
      model.getNotes ().add ("Model with " + size + " reactions");
      FunctionDefinition massAction = new FunctionDefinition ("Function", "Function");
      massAction.setMath ("<math:math><math:lambda><math:bvar><math:ci>K</math:ci></math:bvar><math:bvar><math:ci>S</math:ci>" +
         "</math:bvar><math:apply><math:times/><math:ci>K</math:ci><math:ci>S</math:ci></math:apply></math:lambda></math:math>");
      model.addFunctionDefinition (massAction);
      Compartment compartment = new Compartment ("Compartment", "Compartment");
      compartment.getNotes ().add ("Everything");
      compartment.setSize (1);
      model.addCompartment (compartment);
      Species species = new Species ("Source", "Source");
      species.getNotes ().add (species.getName ());
      species.setCompartment (compartment);
      species.setInitialAmount (1.0);
      model.addSpecies (species);
      List oldSpecies = new ArrayList ();
      oldSpecies.add (species);
      while (size > 0) {
         String speciesName = "Species" + oldSpecies.size ();
         Species newSpecies = new Species (speciesName, speciesName);
         newSpecies.getNotes ().add (speciesName);
         newSpecies.setCompartment (compartment);
         newSpecies.setInitialAmount (0.1);
         for (Iterator iterator = oldSpecies.iterator (); size > 0 && iterator.hasNext (); size--) {
            species = (Species) iterator.next ();
            String reactionName = "Reaction" + size;
            Reaction reaction = new Reaction (reactionName, reactionName);
            reaction.setReversible (false);
            reaction.getNotes ().add (reactionName);
            reaction.addReactant (species);
            String parameterName = "Parameter" + size;
            Parameter parameter = new Parameter (parameterName, parameterName, 1.0 / (size + 2.0));
            parameter.getNotes ().add (parameterName);
            KineticLaw kineticLaw = new KineticLaw ("<math:math><math:apply><math:ci>" + massAction.getId () + "</math:ci><math:ci>" +
               parameter.getId () + "</math:ci><math:ci>" + species.getId () + "</math:ci></math:apply></math:math>");
            model.addParameter (parameter);
            reaction.setKineticLaw (kineticLaw);
            reaction.addProduct (newSpecies);
            model.addReaction (reaction);
         }
         model.addSpecies (newSpecies);
         oldSpecies.add (newSpecies);
      }
      return new SBMLLevel2Document (model);
   }
}
