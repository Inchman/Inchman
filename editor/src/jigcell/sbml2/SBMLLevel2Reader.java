package jigcell.sbml2;

import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.Map;
import jigcell.sbml2.jep.JEP;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * Reads an SBML Level 2 file.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

final class SBMLLevel2Reader extends XMLReader {
   private final static int SBML = NORMAL;
   private final static int ANNOTATIONS = 1;
   private final static int NOTES = 2;
   private final static int RDF = 3;
   private final static int MATH = 4;
   private final static Map elementNameToClass;
   private final static String NAMESPACE_MATHML = "http://www.w3.org/1998/Math/MathML";
   private final static String NAMESPACE_RDF = "http://www.w3.org/1999/02/22-rdf-syntax-ns#";
   private final static String NAMESPACE_SBML = "http://www.sbml.org/sbml/level2";

   private final SBMLLevel2Document document;

   static {
      elementNameToClass = new HashMap ();
      elementNameToClass.put ("algebraicRule", AlgebraicRule.class);
      elementNameToClass.put ("assignmentRule", AssignmentRule.class);
      elementNameToClass.put ("compartment", Compartment.class);
      elementNameToClass.put ("event", Event.class);
      elementNameToClass.put ("eventAssignment", EventAssignment.class);
      elementNameToClass.put ("functionDefinition", FunctionDefinition.class);
      elementNameToClass.put ("kineticLaw", KineticLaw.class);
      elementNameToClass.put ("modifierSpeciesReference", ModifierSpeciesReference.class);
      elementNameToClass.put ("parameter", Parameter.class);
      elementNameToClass.put ("port", Port.class);
      elementNameToClass.put ("rateRule", RateRule.class);
      elementNameToClass.put ("reaction", Reaction.class);
      elementNameToClass.put ("species", Species.class);
      elementNameToClass.put ("speciesReference", SpeciesReference.class);
      elementNameToClass.put ("stoichiometryMath", StoichiometryMath.class);
      elementNameToClass.put ("unit", Unit.class);
      elementNameToClass.put ("unitDefinition", UnitDefinition.class);
      elementNameToClass.put ("instance", Instance.class);
      elementNameToClass.put ("link", Link.class);
      elementNameToClass.put ("target", TargetObjectReference.class);
      elementNameToClass.put ("from", FromObjectReference.class);
      elementNameToClass.put ("to", ToObjectReference.class);
      elementNameToClass.put ("subobject", SubObjectReference.class);
      elementNameToClass.put ("model", Model.class);
   }

   public static SBMLLevel2Document read (Reader reader) throws IOException {
      SBMLLevel2Document document = new SBMLLevel2Document ();
      read (new SBMLLevel2Reader (document), reader);
      return document;
   }

   public static SBMLLevel2Document read (Reader reader,Model model) throws IOException {
      SBMLLevel2Document document = new SBMLLevel2Document (model);
      read (new SBMLLevel2Reader (document), reader);
      return document;
   }

   /**
    * Provides special handling of characters for NOTES mode.  Overrides the
    * XMLReader.characters method which cannot correctly handle NOTES mode.
    */
   public void characters (char characters [], int start, int length) {
      if (localText == null)
         return;
      String text = new String (characters, start, length).trim ();
      if (text.length () != 0)
         if ( mode == NOTES )
            localText.append (XMLPrinter.quote(text));
         else
            localText.append (text);
   }

   public void endElement (String namespace, String localName, String qName) throws SAXException {
      if (mode == SBML && namespace.startsWith(NAMESPACE_SBML)) // changed from .equals to .startsWith - to allow for different versions of level 2 (Aidan, Sep 2010)
         endSBMLElement (localName);
      else if (mode == ANNOTATIONS)
         if (localName.equals ("annotation") && namespace.startsWith (NAMESPACE_SBML)) // changed from .equals to .startsWith - to allow for different versions of level 2 (Aidan, Sep 2010)
            ((SBase) peek ()).getAnnotations ().add (endLocalTextElement () + "</annotation>");
         else
            localText.append ("</" + qName + ">");
      else if (mode == NOTES)
         if (localName.equals ("notes") && namespace.startsWith (NAMESPACE_SBML)) // changed from .equals to .startsWith - to allow for different versions of level 2 (Aidan, Sep 2010)
            ((SBase) peek ()).getNotes ().add (endLocalTextElement () + "</notes>");
         else
            localText.append ("</" + qName + ">");
      else if (mode == RDF)
         if (localName.equals ("RDF") && namespace.startsWith (NAMESPACE_RDF))
            ((SBase) peek ()).setRDF (endLocalTextElement () + "</rdf:RDF>");
         else
            localText.append ("</" + qName + ">");
      else if (mode == MATH && namespace.startsWith (NAMESPACE_MATHML))
         if (localName.equals ("math"))
            setTextOfMathElement (endLocalTextElement () + "</math:math>\n");
         else
            localText.append (JEP.isMathMLFunction (localName) ? ("<math:" + localName + "/>\n") : ("</math:" + localName + ">\n"));
   }

   public void startElement (String namespace, String localName, String qName, Attributes attributes) throws SAXException {
      if (mode == SBML && namespace.startsWith (NAMESPACE_SBML)) // changed from .equals to .startsWith - to allow for different versions of level 2 (Aidan, Sep 2010)
         startSBMLElement (localName, attributes);
      else if (mode == NOTES || mode == ANNOTATIONS || mode == RDF)
         localText.append ("<" + qName + createAttributesText (attributes, namespacePrefixToURI) + ">");
      else if (namespace.equals (NAMESPACE_RDF)) {
         if (mode != SBML || !localName.equals ("RDF"))
            throw new SAXException ("Unexpected RDF tag " + localName + " encountered.");
         mode = RDF;
         ((SBase) peek ()).setRDF ("<rdf:RDF" + createAttributesText (attributes, namespacePrefixToURI) + ">");
      } else if (namespace.equals (NAMESPACE_MATHML)) {
         if (mode != MATH)
            if (localName.equals ("math")) {
               mode = MATH;
               localText = new StringBuffer ("<math:math>\n");
            } else
               throw new SAXException ("Unexpected MathML tag " + localName + " encountered.");
         else if (mode == MATH && !JEP.isMathMLFunction (localName)) {
            localText.append ("<math:" + localName + createAttributesText (attributes, null) + ">");
            if (!localName.equals ("ci") && !localName.equals ("cn"))
               localText.append ("\n");
         }
      } else
         throw new SAXException ("Found unexpected element [" + namespace + "]" + localName + ".");
      namespacePrefixToURI.clear ();
   }

   private SBMLLevel2Reader (SBMLLevel2Document document) {
      super ();
      this.document = document;
   }

   private void endSBMLElement (String localName) throws SAXException {
      XMLElement element = pop ();
      if (localName.equals ("algebraicRule"))
         ((Model) peekPastContainer ()).addRule ((Rule) element);
      else if (localName.equals ("assignmentRule"))
         ((Model) peekPastContainer ()).addRule ((Rule) element);
      else if (localName.equals ("compartment"))
         ((Model) peekPastContainer ()).addCompartment ((Compartment) element);
      else if (localName.equals ("event"))
         ((Model) peekPastContainer ()).addEvent ((Event) element);
      else if (localName.equals ("eventAssignment"))
         ((Event) peekPastContainer ()).getEventAssignment ().add ((EventAssignment) element);
      else if (localName.equals ("functionDefinition"))
         ((Model) peekPastContainer ()).addFunctionDefinition ((FunctionDefinition) element);
      else if (localName.equals ("kineticLaw"))
         ((Reaction) peek ()).setKineticLaw ((KineticLaw) element);
      else if (localName.equals ("modifierSpeciesReference")) {
         XMLElement localElement = peek ();
         Reaction reaction = (Reaction) peekPastContainer ();
         if (reaction.getModifiersElement () == localElement)
            reaction.getModifier ().add ((ModifierSpeciesReference) element);
         else
            throw new SAXException ("ModifierSpeciesReference can only be used in the list of modifiers.");
      } else if (localName.equals ("parameter")) {
         XMLElement localElement = peekPastContainer ();
         if (localElement instanceof KineticLaw)
            ((KineticLaw) localElement).getParameter ().add ((Parameter) element);
         else
            ((Model) localElement).addParameter ((Parameter) element);
      } else if (localName.equals ("rateRule"))
         ((Model) peekPastContainer ()).addRule ((Rule) element);
      else if (localName.equals ("reaction"))
         ((Model) peekPastContainer ()).addReaction ((Reaction) element);
      else if (localName.equals ("species"))
         ((Model) peekPastContainer ()).addSpecies ((Species) element);
      else if (localName.equals ("speciesReference")) {
         XMLElement localElement = peek ();
         Reaction reaction = (Reaction) peekPastContainer ();
         if (reaction.getReactantsElement () == localElement)
            reaction.getReactant ().add ((SpeciesReference) element);
         else if (reaction.getProductsElement () == localElement)
            reaction.getProduct ().add ((SpeciesReference) element);
         else
            throw new SAXException ("SpeciesReference can only be used in the list of reactants or list of products.");
      } else if (localName.equals ("stoichiometryMath"))
         ((SpeciesReference) peek ()).setStoichiometryMath ((StoichiometryMath) element);
      else if (localName.equals ("unit"))
         ((UnitDefinition) peekPastContainer ()).getUnits ().add ((Unit) element);
      else if (localName.equals ("unitDefinition"))
         ((Model) peekPastContainer ()).addUnitDefinition ((UnitDefinition) element);
      else if (localName.equals ("instance"))
         ((Model) peekPastContainer ()).addInstance ((Instance) element);
      else if (localName.equals ("port"))
         ((Model) peekPastContainer ()).addPort ((Port)element);
      else if (localName.equals ("target"))
         ((Port) peek ()).addTargetObject ((TargetObjectReference) element);
      else if (localName.equals ("link"))
         ((Model) peekPastContainer ()).addLink ((Link)element);
      else if (localName.equals ("from"))
         ((Link) peek ()).addFromObject ((FromObjectReference) element);
      else if (localName.equals ("to"))
         ((Link) peek ()).addToObject ((ToObjectReference) element);
      else if (localName.equals ("subobject")) {
         XMLElement localElement = peek ();
         if (localElement instanceof FromObjectReference)
            ((FromObjectReference)localElement).setSubObjectReference ((SubObjectReference) element);
         else if (localElement instanceof TargetObjectReference)
            ((TargetObjectReference)localElement).setSubObjectReference ((SubObjectReference) element);
         else if (localElement instanceof ToObjectReference)
            ((ToObjectReference)localElement).setSubObjectReference ((SubObjectReference) element);
         else
            throw new SAXException ("SubObjectReference can only be used in the from and to structures of a link or port.");
      }
      else if (localName.equals ("model")) {
         // Check if this model is a submodel and add it to the list of submodels
         if (!(peek () instanceof SBMLLevel2Document))
            ((Model) peekPastContainer ()).addModel ((Model)element);
      }
   }

   private void setTextOfMathElement (String text) throws SAXException {
      XMLElement element = peek ();
      if (element instanceof MathElement) {
         ((MathElement) element).setMath (text);
         return;
      }
      Event event = (Event) peekPastContainer ();
      if (event.getDelayElement () == element)
         event.setDelay (text);
      else if (event.getTriggerElement () == element)
         event.setTrigger (text);
      else
         throw new SAXException ("Can only set math for an event when in delay, trigger, or assignment.");
   }

   private void startSBMLElement (String localName, Attributes attributes) throws SAXException {
      if (localName.equals ("notes")) {
         mode = NOTES;
         localText = new StringBuffer ("<notes" + createAttributesText (attributes, namespacePrefixToURI) + ">");
         return;
      }
      if (localName.equals ("annotation")) {
         mode = ANNOTATIONS;
         localText = new StringBuffer ("<annotation" + createAttributesText (attributes, namespacePrefixToURI) + ">");
         return;
      }
      SBase element;
      Class elementClass = (Class) elementNameToClass.get (localName);
      // Check for the main root model by checking whether localName equals "model"
      if (elementClass != null && !localName.equals ("model"))
         try {
            element = (SBase) elementClass.newInstance ();
         } catch (Exception e) {
            throw (SAXException) new SAXException ("Unable to create " + localName + " element.").initCause (e);
         }
      else if (localName.equals ("delay"))
         element = ((Event) peek ()).getDelayElement ();
      else if (localName.equals ("listOfCompartments"))
         element = ((Model) peek ()).getCompartmentsElement ();
      else if (localName.equals ("listOfEventAssignments"))
         element = ((Event) peek ()).getAssignmentsElement ();
      else if (localName.equals ("listOfEvents"))
         element = ((Model) peek ()).getEventsElement ();
      else if (localName.equals ("listOfFunctionDefinitions"))
         element = ((Model) peek ()).getFunctionDefinitionsElement ();
      else if (localName.equals ("listOfModifiers"))
         element = ((Reaction) peek ()).getModifiersElement ();
      else if (localName.equals ("listOfParameters")) {
         XMLElement localElement = peek ();
         element = localElement instanceof KineticLaw ? ((KineticLaw) localElement).getParameterElement () :
            ((Model) peek ()).getParametersElement ();
      } else if (localName.equals ("listOfProducts"))
         element = ((Reaction) peek ()).getProductsElement ();
      else if (localName.equals ("listOfReactants"))
         element = ((Reaction) peek ()).getReactantsElement ();
      else if (localName.equals ("listOfReactions"))
         element = ((Model) peek ()).getReactionsElement ();
      else if (localName.equals ("listOfRules"))
         element = ((Model) peek ()).getRulesElement ();
      else if (localName.equals ("listOfSpecies"))
         element = ((Model) peek ()).getSpeciesElement ();
      else if (localName.equals ("listOfUnitDefinitions"))
         element = ((Model) peek ()).getUnitDefinitionsElement ();
      else if (localName.equals ("listOfUnits"))
         element = ((UnitDefinition) peek ()).getUnitsElement ();
      else if (localName.equals ("listOfInstances"))
         element = ((Model) peek ()).getInstancesElement ();
      else if (localName.equals ("listOfPorts"))
         element = ((Model) peek ()).getPortsElement ();
      else if (localName.equals ("listOfLinks"))
         element = ((Model) peek ()).getLinksElement ();
      else if (localName.equals ("listOfSubmodels"))
         element = ((Model) peek ()).getModelsElement ();
      else if (localName.equals ("model")) {
         Model model;
         // Only the main root model should be set as the document's model
         if (peek () instanceof SBMLLevel2Document) {
            model = document.getModel();
            if ( model == null ) {
               model = new Model();
               document.setModel (model);
            }
         } else {
            model = new Model ();
         }
         element = model;

      } else if (localName.equals ("sbml"))
         element = document;
      else if (localName.equals ("trigger"))
         element = ((Event) peek ()).getTriggerElement ();
      else
         throw new SAXException ("Found unknown SBML element " + localName + ".");
      element.parse (attributes);
      push (element);
   }
}
