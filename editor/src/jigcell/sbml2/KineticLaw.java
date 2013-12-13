package jigcell.sbml2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.xml.sax.Attributes;

/**
 * Describes the rate at which a reaction takes place.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class KineticLaw extends SBase implements MathElement {
   protected List parameter;
   private final SBase parameterElement;
   private String math;
   private String substanceUnits;
   private String timeUnits;

   public KineticLaw () {
      this ((String) null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param law Kinetic law, must not be null
    */

   public KineticLaw (KineticLaw law) {
      this (law.getMath ());
      setSubstanceUnits (law.getSubstanceUnits ());
      setTimeUnits (law.getTimeUnits ());
      for (Iterator iterator = law.getParameter ().iterator (); iterator.hasNext (); ) {
         Parameter parameter = (Parameter) iterator.next ();
         if (parameter != null)
            addParameter (new Parameter (parameter));
      }
   }

   public KineticLaw (String math) {
      setMath (math);
      parameterElement = new SBase ();
      parameter = new ArrayList ();
   }

   public void addParameter (Parameter element) {
      if (element == null)
         throw new IllegalArgumentException ();
      parameter.add (element);
   }

   public Parameter findParameterWithId (String id) {
      return (Parameter) searchListForId (parameter, id);
   }

   public String getMath () {
      return math;
   }

   public List getParameter () {
      return parameter;
   }

   public SBase getParameterElement () {
      return parameterElement;
   }

   public String getSubstanceUnits () {
      return substanceUnits == null ? "substance" : substanceUnits;
   }

   public UnitDefinition getSubstanceUnits (Model model) {
      return (UnitDefinition) model.findElementWithId (getSubstanceUnits (), Model.UNITDEFINITION);
   }

   public String getTimeUnits () {
      return timeUnits == null ? "time" : timeUnits;
   }

   public UnitDefinition getTimeUnits (Model model) {
      return (UnitDefinition) model.findElementWithId (getTimeUnits (), Model.UNITDEFINITION);
   }

   public boolean isValid (Model model) {
      return super.isValid (model) && UnitDefinition.isValidTimeUnit (getTimeUnits (model)) &&
         UnitDefinition.isValidSubstanceUnit (getSubstanceUnits (model));
   }

   /**
    * Sets the MathML expression for this {@link KineticLaw}.
    *
    * @param math New value of property math.
    */

   public void setMath (String math) {
      assert math == null || math.startsWith ("<math:math>");
      this.math = math;
   }

   public void setSubstanceUnits (String substanceUnits) {
      if (substanceUnits != null && !SBaseId.isValidId (substanceUnits))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.substanceUnits = substanceUnits;
   }

   public void setSubstanceUnits (UnitDefinition units) {
      setSubstanceUnits (units == null ? null : units.getId ());
   }

   public void setTimeUnits (String timeUnits) {
      if (timeUnits != null && !SBaseId.isValidId (timeUnits))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.timeUnits = timeUnits;
   }

   public void setTimeUnits (UnitDefinition units) {
      setTimeUnits (units == null ? null : units.getId ());
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      setSubstanceUnits (attributes.getValue ("substanceUnits"));
      setTimeUnits (attributes.getValue ("timeUnits"));
    }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "kineticLaw");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      if (!getTimeUnits ().equals ("time"))
         printer.addAttribute ("timeUnits", getTimeUnits ());
      if (!getSubstanceUnits ().equals ("substance"))
         printer.addAttribute ("substanceUnits", getSubstanceUnits ());
      printer.addText (getMath ());
      printer.addElementList (getParameterElement (), "listOfParameters", getParameter ());
      return printer;
   }
}
