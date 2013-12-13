package jigcell.sbml2;

import org.xml.sax.Attributes;

/**
 * Declares a variable for use in mathematical formulas in an SBML model definition.  Parameters have a constant value by default.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * <p>
 *
 * @author Nicholas Allen
 */

public final class Parameter extends SBaseId {
   private boolean constant;
   private double value;
   private String units;

   public Parameter () {
      this (null, null, Double.NaN);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param parameter Parameter, must not be null
    */

   public Parameter (Parameter parameter) {
      this (parameter.getId (), parameter.isNameSet () ? parameter.getName () : null);
      setConstant (parameter.isConstant ());
      setValue (parameter.getValue ());
      setUnits (parameter.getUnits ());
   }

   public Parameter (String id, String name) {
      this (id, name, Double.NaN);
   }

   public Parameter (String id, String name, double value) {
      super (id, name);
      setValue (value);
      setConstant (true);
   }

   public String getUnits () {
      return units;
   }

   public UnitDefinition getUnits (Model model) {
      return (UnitDefinition) model.findElementWithId (getUnits (), Model.UNITDEFINITION);
   }

   public double getValue () {
      return value;
   }

   public boolean isConstant () {
      return constant;
   }

   public boolean isValid (Model model) {
      return super.isValid (model) && getUnits (model) != null;
   }

   /**
    * Sets whether the parameter's value is constant throughout a simulation.
    *
    * @param constant New value of property constant.
    */

   public void setConstant (boolean constant) {
      this.constant = constant;
   }

   public void setUnits (String units) {
      if (units != null && !isValidId (units))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.units = units;
   }

   public void setUnits (UnitDefinition units) {
      setUnits (units == null ? null : units.getId ());
   }

   public void setValue (double value) {
      this.value = value;
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      if (attributes.getIndex ("constant") != -1)
         setConstant (Boolean.valueOf (attributes.getValue ("constant")) == Boolean.TRUE);
      setUnits (attributes.getValue ("units"));
      if (attributes.getIndex ("value") != -1)
         setValue (Double.parseDouble (attributes.getValue ("value")));
    }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "parameter");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      if (!Double.isNaN (getValue ()))
         printer.addAttribute ("value", String.valueOf (getValue ()));
      if (!isConstant ())
         printer.addAttribute ("constant", "false");
      printer.addAttribute ("units", getUnits ());
      return printer;
   }
}
