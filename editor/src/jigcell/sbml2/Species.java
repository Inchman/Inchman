package jigcell.sbml2;

import java.util.Iterator;
import org.xml.sax.Attributes;

/**
 * A chemical entity that takes part in a reaction.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public final class Species extends SBaseId {
   private boolean boundaryCondition;
   private boolean constant;
   private boolean hasOnlySubstanceUnits;
   private double initialAmount;
   private double initialConcentration;
   private Integer charge;
   private String compartment;
   private String spatialSizeUnits;
   private String substanceUnits;

   public Species () {
      this (null, null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param species Species, must not be null
    */

   public Species (Species species) {
      this (species.getId (), species.isNameSet () ? species.getName () : null);
      setBoundaryCondition (species.isBoundaryCondition ());
      setConstant (species.isConstant ());
      setHasOnlySubstanceUnits (species.getHasOnlySubstanceUnits ());
      setInitialAmount (species.getInitialAmount ());
      setInitialConcentration (species.getInitialConcentration ());
      if (species.isChargeSet ())
         setCharge (species.getCharge ());
      setCompartment (species.getCompartment ());
      setSpatialSizeUnits (species.getSpatialSizeUnits ());
      setSubstanceUnits (species.getSubstanceUnits ());
   }

   public Species (String id, String name) {
      super (id, name);
      setInitialAmount (Double.NaN);
      setInitialConcentration (Double.NaN);
      setBoundaryCondition (false);
      setConstant (false);
      setHasOnlySubstanceUnits (false);
      unsetCharge ();
   }

   public int getCharge () {
      assert isChargeSet ();
      return isChargeSet () ? charge.intValue () : 0;
   }

   public String getCompartment () {
      return compartment;
   }

   public Compartment getCompartment (Model model) {
      return (Compartment) model.findElementWithId (getCompartment (), Model.COMPARTMENT);
   }

   public boolean getHasOnlySubstanceUnits () {
      return hasOnlySubstanceUnits;
   }

   public double getInitialAmount () {
      return initialAmount;
   }

   public double getInitialConcentration () {
      return getHasOnlySubstanceUnits () ? Double.NaN : initialConcentration;
   }

   public String getSpatialSizeUnits () {
      return spatialSizeUnits;
   }

   public UnitDefinition getSpatialSizeUnits (Model model) {
      UnitDefinition units = (UnitDefinition) model.findElementWithId (getSpatialSizeUnits (), Model.UNITDEFINITION);
      if (units != null)
         return units;
      Compartment compartment = getCompartment (model);
      return compartment == null ? null : compartment.getUnits (model);
   }

   public String getSubstanceUnits () {
      return substanceUnits == null ? "substance" : substanceUnits;
   }

   public UnitDefinition getSubstanceUnits (Model model) {
      return (UnitDefinition) model.findElementWithId (getSubstanceUnits (), Model.UNITDEFINITION);
   }

   public boolean isBoundaryCondition () {
      return boundaryCondition;
   }

   public boolean isChargeSet () {
      return charge != null;
   }

   public boolean isConstant () {
      return constant;
   }

   public boolean isSetByEvent (Model model) {
      String id = getId ();
      if (id == null)
         return false;
      for (Iterator eventIterator = model.getEvents ().iterator (); eventIterator.hasNext (); )
         for (Iterator assignmentIterator = ((Event) eventIterator.next ()).getEventAssignment ().iterator (); assignmentIterator.hasNext ();
            ) {
         if (id.equals (((EventAssignment) assignmentIterator.next ()).getVariable ()))
            return true;
      }
      return false;
   }

   public boolean isSetByRule (Model model) {
      String id = getId ();
      if (id == null)
         return false;
      for (Iterator iterator = model.getRules ().iterator (); iterator.hasNext (); ) {
         Rule rule = (Rule) iterator.next ();
         if (rule instanceof VariableRule && id.equals (((VariableRule) rule).getVariable ()))
            return true;
      }
      return false;
   }

   public boolean isValid (Model model) {
      if (!super.isValid (model) || !UnitDefinition.isValidSubstanceUnit (getSubstanceUnits (model)) || getCompartment (model) == null)
         return false;
      Compartment compartment = getCompartment (model);
      if (compartment == null)
         return false;
      int dimensions = compartment.getSpatialDimensions ();
      if (!UnitDefinition.isValidSpatialSizeUnit (getSpatialSizeUnits (model), dimensions))
         return false;
      return dimensions != 0 || Double.isNaN (getInitialConcentration ());
   }

   /**
    * Sets whether the species is on the boundary of the reaction system.
    *
    * @param boundaryCondition New value of property boundaryCondition.
    */

   public void setBoundaryCondition (boolean boundaryCondition) {
      this.boundaryCondition = boundaryCondition;
   }

   public void setCharge (int charge) {
      this.charge = new Integer (charge);
   }

   public void setCompartment (String compartment) {
      if (compartment != null && !isValidId (compartment))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.compartment = compartment;
   }

   public void setCompartment (Compartment compartment) {
      setCompartment (compartment == null ? null : compartment.getId ());
   }

   public void setConstant (boolean constant) {
      this.constant = constant;
   }

   public void setHasOnlySubstanceUnits (boolean hasOnlySubstanceUnits) {
      this.hasOnlySubstanceUnits = hasOnlySubstanceUnits;
   }

   public void setInitialAmount (double initialAmount) {
      this.initialAmount = initialAmount;
      if (!Double.isNaN (initialAmount))
         setInitialConcentration (Double.NaN);
   }

   public void setInitialConcentration (double initialConcentration) {
      this.initialConcentration = initialConcentration;
      if (!Double.isNaN (initialConcentration))
         setInitialAmount (Double.NaN);
   }

   public void setSpatialSizeUnits (String spatialSizeUnits) {
      if (spatialSizeUnits != null && !isValidId (spatialSizeUnits))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.spatialSizeUnits = spatialSizeUnits;
   }

   public void setSpatialSizeUnits (UnitDefinition units) {
      setSpatialSizeUnits (units == null ? null : units.getId ());
   }

   public void setSubstanceUnits (String substanceUnits) {
      if (substanceUnits != null && !isValidId (substanceUnits))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.substanceUnits = substanceUnits;
   }

   public void setSubstanceUnits (UnitDefinition units) {
      setSubstanceUnits (units == null ? null : units.getId ());
   }

   public void unsetCharge () {
      charge = null;
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      if (attributes.getIndex ("boundaryCondition") != -1)
         setBoundaryCondition (Boolean.valueOf (attributes.getValue ("boundaryCondition")) == Boolean.TRUE);
      if (attributes.getIndex ("charge") != -1)
         setCharge (Integer.parseInt (attributes.getValue ("charge")));
      setCompartment (attributes.getValue ("compartment"));
      if (attributes.getIndex ("constant") != -1)
         setConstant (Boolean.valueOf (attributes.getValue ("constant")) == Boolean.TRUE);
      if (attributes.getIndex ("hasOnlySubstanceUnits") != -1)
         setHasOnlySubstanceUnits (Boolean.valueOf (attributes.getValue ("hasOnlySubstanceUnits")) == Boolean.TRUE);
      if (attributes.getIndex ("initialAmount") != -1)
         setInitialAmount (Double.parseDouble (attributes.getValue ("initialAmount")));
      if (attributes.getIndex ("initialConcentration") != -1)
         setInitialConcentration (Double.parseDouble (attributes.getValue ("initialConcentration")));
      setSpatialSizeUnits (attributes.getValue ("spatialSizeUnits"));
      setSubstanceUnits (attributes.getValue ("substanceUnits"));
    }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "species");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addAttribute ("compartment", getCompartment ());
      if (!Double.isNaN (getInitialAmount ()))
         printer.addAttribute ("initialAmount", String.valueOf (getInitialAmount ()));
      if (!Double.isNaN (getInitialConcentration ()))
         printer.addAttribute ("initialConcentration", String.valueOf (getInitialConcentration ()));
      if (isBoundaryCondition ())
         printer.addAttribute ("boundaryCondition", "true");
      if (isChargeSet ())
         printer.addAttribute ("charge", String.valueOf (getCharge ()));
      if (isConstant ())
         printer.addAttribute ("constant", "true");
      if (getHasOnlySubstanceUnits ())
         printer.addAttribute ("hasOnlySubstanceUnits", "true");
      printer.addAttribute ("spatialSizeUnits", getSpatialSizeUnits ());
      if (!getSubstanceUnits ().equals ("substance"))
         printer.addAttribute ("substanceUnits", getSubstanceUnits ());
      return printer;
   }
}
