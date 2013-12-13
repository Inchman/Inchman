package jigcell.sbml2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * The highest level construct in an SBML document.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 * @author Ranjit Randhawa
 */

public class Model extends SBaseId {
   // If you add to this list, be sure to update GLOBAL accordingly.  The SBML
   // spec says that FunctionDefinition, CompartmentType, SpeciesType,
   // Compartment, Species, Parameter, Reaction, SpeciesReference,
   // ModifierSpeciesReference, Event, and Model are in the same scope.
   // UnitDefinition is separate.
   // Each reaction also has it's own scope for local parameters.
   // The SUBMODEL, INSTANCE, and PORT are for aggregate models, not supported
   // in the SBML spec.  Therefore, we don't include them in GLOBAL.
   // ALL represents all flags.
   public final static int COMPARTMENT = 1;
   public final static int EVENT = 2;
   public final static int FUNCTIONDEFINITION = 4;
   public final static int PARAMETER = 8;
   public final static int REACTION = 16;
   public final static int SPECIES = 32;
   public final static int UNITDEFINITION = 64;
   public final static int SUBMODEL = 128;
   public final static int INSTANCE = 256;
   public final static int PORT = 516;
   public final static int MODEL = 1024;
   public final static int GLOBAL = COMPARTMENT | EVENT | FUNCTIONDEFINITION | PARAMETER | REACTION | SPECIES | MODEL;
   public final static int ALL = GLOBAL | UNITDEFINITION | SUBMODEL | INSTANCE | PORT;


   protected List compartments;
   protected List events;
   protected List functionDefinitions;
   protected List instances;
   protected List links;
   protected List models;
   protected List parameters;
   protected List ports;
   protected List reactions;
   protected List rules;
   protected List species;
   protected List unitDefinitions;
   protected SBase compartmentsElement;
   protected SBase eventsElement;
   protected SBase functionDefinitionsElement;
   protected SBase instancesElement;
   protected SBase linksElement;
   protected SBase modelsElement;
   protected SBase parametersElement;
   protected SBase portsElement;
   protected SBase reactionsElement;
   protected SBase rulesElement;
   protected SBase speciesElement;
   protected SBase unitDefinitionsElement;

   public Model () {
      this (null);
   }

   public Model (String name) {
      super (null, name);
      compartments = new ArrayList ();
      events = new ArrayList ();
      functionDefinitions = new ArrayList ();
      instances = new ArrayList ();
      links = new ArrayList ();
      models = new ArrayList ();
      parameters = new ArrayList ();
      ports = new ArrayList ();
      reactions = new ArrayList ();
      rules = new ArrayList ();
      species = new ArrayList ();
      unitDefinitions = new ArrayList ();
      compartmentsElement = new SBase ();
      eventsElement = new SBase ();
      functionDefinitionsElement = new SBase ();
      instancesElement = new SBase ();
      linksElement = new SBase ();
      modelsElement = new SBase ();
      parametersElement = new SBase ();
      portsElement = new SBase ();
      reactionsElement = new SBase ();
      rulesElement = new SBase ();
      speciesElement = new SBase ();
      unitDefinitionsElement = new SBase ();
   }

   public void addCompartment (Compartment element) {
      if (element == null)
         throw new IllegalArgumentException ();
      compartments.add (element);
   }

   public void addEvent (Event element) {
      if (element == null)
         throw new IllegalArgumentException ();
      events.add (element);
   }

   public void addFunctionDefinition (FunctionDefinition element) {
      if (element == null)
         throw new IllegalArgumentException ();
      functionDefinitions.add (element);
   }

   public void addInstance (Instance element) {
      if (element == null)
         throw new IllegalArgumentException ();
      instances.add (element);
   }

   public void addLink (Link element) {
      if (element == null)
         throw new IllegalArgumentException ();
      links.add (element);
   }
   
   public void addModel (Model element) {
      if (element == null)
         throw new IllegalArgumentException ();
      models.add (element);
   }
   
   public void addParameter (Parameter element) {
      if (element == null)
         throw new IllegalArgumentException ();
      parameters.add (element);
   }
   
   public void addPort (Port element) {
      if (element == null)
         throw new IllegalArgumentException ();
      ports.add (element);
   }
   
   public void addReaction (Reaction element) {
      if (element == null)
         throw new IllegalArgumentException ();
      reactions.add (element);
   }

   public void addRule (Rule element) {
      if (element == null)
         throw new IllegalArgumentException ();
      rules.add (element);
   }

   public void addSpecies (Species element) {
      if (element == null)
         throw new IllegalArgumentException ();
      species.add (element);
   }

   public void addUnitDefinition (UnitDefinition element) {
      if (element == null)
         throw new IllegalArgumentException ();
      String id = element.getId ();
      if (Unit.findBaseUnit (id) != null)
         throw new IllegalArgumentException ("Base units cannot be redefined.");
      unitDefinitions.add (element);
   }

   /**
    * Searches the model for an sbml element with the given id.  The different
    * types of elements to search for will be searched in the order of
    * compartments, events, function definitions, parameters, reactions, 
    * species, and then unit definitions.  Note that this method only
    * searches for global parameters.  To search for local parameters, use
    * {@link jigcell.sbml2.KineticLaw#findParameterWithId(java.lang.String)}.
    *
    * @param id SBML identifier
    * @param flags Bit flags specifying what types of elements to search for
    */

   public SBaseId findElementWithId (String id, int flags) {
      SBaseId element = null;
      if ((flags & COMPARTMENT) == COMPARTMENT)
         element = searchListForId (compartments, id);
      if (element != null)
         return element;
      if ((flags & EVENT) == EVENT)
         element = searchListForId (events, id);
      if (element != null)
         return element;
      if ((flags & FUNCTIONDEFINITION) == FUNCTIONDEFINITION)
         element = searchListForId (functionDefinitions, id);
      if (element != null)
         return element;
      if ((flags & PARAMETER) == PARAMETER)
         element = searchListForId (parameters, id);
      if (element != null)
         return element;
      if ((flags & PORT) == PORT)
         element = searchListForId (ports, id);
      if (element != null)
         return element;
      if ((flags & REACTION) == REACTION)
         element = searchListForId (reactions, id);
      if (element != null)
         return element;
      if ((flags & SPECIES) == SPECIES)
         element = searchListForId (species, id);
      if (element != null)
         return element;
      if ((flags & SUBMODEL) == SUBMODEL)
         element = searchListForId (models, id);
      if (element != null)
         return element;
      if ((flags & INSTANCE) == INSTANCE)
         element = searchListForId (instances, id);
      if (element != null)
         return element;
      if ((flags & MODEL) == MODEL)
         element = searchListForId (models, id);
      if (element != null)
         return element;
      if ((flags & UNITDEFINITION) == UNITDEFINITION)
         element = searchListForId (unitDefinitions, id);
      return element;
   }

   /**
    * Searches the SBML model for an element with the given name.
    * @see {@link jigcell.sbml2.Model#findParameterWithId(java.lang.String,int)}.
    */
   public SBaseId findElementWithName (String name, int flags) {
      SBaseId element = null;
      if ((flags & COMPARTMENT) == COMPARTMENT)
         element = searchListForName (compartments, name);
      if (element != null)
         return element;
      if ((flags & EVENT) == EVENT)
         element = searchListForName (events, name);
      if (element != null)
         return element;
      if ((flags & FUNCTIONDEFINITION) == FUNCTIONDEFINITION)
         element = searchListForName (functionDefinitions, name);
      if (element != null)
         return element;
      if ((flags & PARAMETER) == PARAMETER)
         element = searchListForName (parameters, name);
      if (element != null)
         return element;
      if ((flags & REACTION) == REACTION)
         element = searchListForName (reactions, name);
      if (element != null)
         return element;
      if ((flags & SPECIES) == SPECIES)
         element = searchListForName (species, name);
      if (element != null)
         return element;
      if ((flags & SUBMODEL) == SUBMODEL)
         element = searchListForName (models, name);
      if (element != null)
         return element;
      if ((flags & INSTANCE) == INSTANCE)
         element = searchListForName (instances, name);
      if (element != null)
         return element;
      if ((flags & PORT) == PORT)
         element = searchListForName (ports, name);
      if (element != null)
         return element;
      if ((flags & MODEL) == MODEL)
         element = searchListForName (models, name);
      if (element != null)
         return element;
      if ((flags & UNITDEFINITION) == UNITDEFINITION)
         element = searchListForName (unitDefinitions, name);
      return element;
   }

   public List findReactionsWithModifierId( String id ) {
      List reactions = getReactions();
      List reactionsWithId = new ArrayList();
      for (Iterator i=reactions.iterator(); i.hasNext(); ) {
         Reaction reaction = (Reaction)i.next();
         List modifiers = reaction.getModifier();
         for ( Iterator j=modifiers.iterator(); j.hasNext(); ) {
            ModifierSpeciesReference modifier =
               (ModifierSpeciesReference)j.next();
            if ( modifier.getSpecies().equals(id) ) {
               reactionsWithId.add(reaction);
               break;
            }
         }
      }
      return reactionsWithId;
   }

   public List findReactionsWithProductId( String id ) {
      List reactions = getReactions();
      List reactionsWithId = new ArrayList();
      for (Iterator i=reactions.iterator(); i.hasNext(); ) {
         Reaction reaction = (Reaction)i.next();
         List products = reaction.getProduct();
         for ( Iterator j=products.iterator(); j.hasNext(); ) {
            SpeciesReference product = (SpeciesReference)j.next();
            if ( product.getSpecies().equals(id) ) {
               reactionsWithId.add(reaction);
               break;
            }
         }
      }
      return reactionsWithId;
   }

   public List findReactionsWithReactantId( String id ) {
      List reactions = getReactions();
      List reactionsWithId = new ArrayList();
      for (Iterator i=reactions.iterator(); i.hasNext(); ) {
         Reaction reaction = (Reaction)i.next();
         List reactants = reaction.getReactant();
         for ( Iterator j=reactants.iterator(); j.hasNext(); ) {
            SpeciesReference reactant = (SpeciesReference)j.next();
            if ( reactant.getSpecies().equals(id) ) {
               reactionsWithId.add(reaction);
               break;
            }
         }
      }
      return reactionsWithId;
   }

   /**
    * Searches all the rules for an AssignmentRule or RateRule that defines
    * species with <i>id</i>.
    */
   public VariableRule findRuleWithVariableId( String id ) {
      for ( Iterator i = rules.iterator(); i.hasNext(); ) {
         Rule rule = (Rule)i.next();
         if ( !(rule instanceof VariableRule) ) continue;
         VariableRule vr = (VariableRule)rule;
         if ( id.equals(vr.getVariable()) ) return vr;
      }
      return null;
   }

   public List getCompartments () {
      return compartments;
   }

   public List getCompartmentsWithIds () {
      return getElementsWithIds(compartments);
   }

   public SBase getCompartmentsElement () {
      return compartmentsElement;
   }

   /**
    * Returns a list of elements that have non-null ids.  JigCell creates
    * elements with null ids as place holders for empty rows in the
    * ModelBuilder.  This allows the code to easily get a list of just those
    * elements with ids.
    */
   private List getElementsWithIds ( List elements ) {
      List results = new ArrayList();
      for ( Iterator i = elements.iterator(); i.hasNext(); ) {
         SBaseId element = (SBaseId)i.next();
         if ( element.getId() != null ) results.add(element);
      }
      return results;
   }

   public List getEvents () {
      return events;
   }

   public List getEventsWithIds () {
      return getElementsWithIds(events);
   }

   public SBase getEventsElement () {
      return eventsElement;
   }

   public List getFunctionDefinitions () {
      return functionDefinitions;
   }

   public List getFunctionDefinitionsWithIds () {
      return getElementsWithIds(functionDefinitions);
   }

   public SBase getFunctionDefinitionsElement () {
      return functionDefinitionsElement;
   }

   public List getInstances () {
      return instances;
   }

   public SBase getInstancesElement () {
      return instancesElement;
   }
   
   public List getLinks () {
      return links;
   }

   public SBase getLinksElement () {
      return linksElement;
   }
   
   public List getModels () {
      return models;
   }

   public SBase getModelsElement () {
      return modelsElement;
   }

   public List getParameters () {
      return parameters;
   }

   public List getParametersWithIds () {
      return getElementsWithIds(parameters);
   }

   public SBase getParametersElement () {
      return parametersElement;
   }

   public List getPorts () {
      return ports;
   }

   public SBase getPortsElement () {
      return portsElement;
   }
   
   public List getReactions () {
      return reactions;
   }

   public List getReactionsWithIds () {
      return getElementsWithIds(reactions);
   }

   public SBase getReactionsElement () {
      return reactionsElement;
   }

   public List getRules () {
      return rules;
   }

   public SBase getRulesElement () {
      return rulesElement;
   }

   public List getSpecies () {
      return species;
   }

   public List getSpeciesWithIds () {
      return getElementsWithIds(species);
   }

   public SBase getSpeciesElement () {
      return speciesElement;
   }

   public List getUnitDefinitions () {
      return unitDefinitions;
   }

   public List getUnitDefinitionsWithIds () {
      return getElementsWithIds(unitDefinitions);
   }

   public SBase getUnitDefinitionsElement () {
      return unitDefinitionsElement;
   }

   public boolean hasLocalParameters () {
      for (Iterator iterator = reactions.iterator (); iterator.hasNext (); ) {
         KineticLaw kineticLaw = ((Reaction) iterator.next ()).getKineticLaw ();
         if (kineticLaw != null && kineticLaw.getParameter ().size () > 0)
            return true;
      }
      return false;
   }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "model");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addElementList (getFunctionDefinitionsElement (), "listOfFunctionDefinitions", getFunctionDefinitions ());
      printer.addElementList (getUnitDefinitionsElement (), "listOfUnitDefinitions", getUnitDefinitions ());
      printer.addElementList (getCompartmentsElement (), "listOfCompartments", getCompartments ());
      printer.addElementList (getSpeciesElement (), "listOfSpecies", getSpecies ());
      printer.addElementList (getParametersElement (), "listOfParameters", getParameters ());
      printer.addElementList (getRulesElement (), "listOfRules", getRules ());
      printer.addElementList (getReactionsElement (), "listOfReactions", getReactions ());
      printer.addElementList (getEventsElement (), "listOfEvents", getEvents ());
      printer.addElementList (getModelsElement (), "listOfSubmodels", getModels ());
      printer.addElementList (getInstancesElement (), "listOfInstances", getInstances ());
      printer.addElementList (getPortsElement (), "listOfPorts", getPorts ());
      printer.addElementList (getLinksElement (), "listOfLinks", getLinks ());
      return printer;
   }
}
