package jigcell.sbml2;

import org.xml.sax.Attributes;

/**
 * A bounded space in which species are located.  Compartments do not have to correspond to actual cellular structures.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public final class Compartment extends SBaseId {
   private boolean constant;
   private double size;
   private int spatialDimensions;
   private String outside;
   private String units;

   public Compartment () {
      this (null, null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param compartment Compartment, must not be null
    */

   public Compartment (Compartment compartment) {
      this (compartment.getId (), compartment.isNameSet () ? compartment.getName () : null);
      setConstant (compartment.isConstant ());
      setSize (compartment.getSize ());
      setSpatialDimensions (compartment.getSpatialDimensions ());
      setOutside (compartment.getOutside ());
      setUnits (compartment.getUnits ());
   }

   public Compartment (String id, String name) {
      super (id, name);
      setConstant (true);
      setSize (Double.NaN);
      setSpatialDimensions (3);
   }

   public String getOutside () {
      return outside;
   }

   public Compartment getOutside (Model model) {
      return (Compartment) model.findElementWithId (getOutside (), Model.COMPARTMENT);
   }

   public double getSize () {
      return spatialDimensions == 0 ? Double.NaN : size;
   }

   public int getSpatialDimensions () {
      return spatialDimensions;
   }

   public String getUnits () {
      return units == null ? getDefaultUnits () : units;
   }

   public UnitDefinition getUnits (Model model) {
      return (UnitDefinition) model.findElementWithId (getUnits (), Model.UNITDEFINITION);
   }

   public boolean isConstant () {
      return spatialDimensions == 0 || constant;
   }

   public boolean isValid (Model model) {
      return super.isValid (model) && UnitDefinition.isValidSpatialSizeUnit (getUnits (model), getSpatialDimensions ()) &&
         getOutside (model) != null;
   }

   /**
    * Sets whether the compartment's volume may vary during simulation.
    */

   public void setConstant (boolean constant) {
      this.constant = constant;
   }

   public void setOutside (String outside) {
      this.outside = outside;
   }

   public void setOutside (Compartment compartment) {
      setOutside (compartment == null ? null : compartment.getId ());
   }

   public void setSize (double size) {
      this.size = size;
   }

   /**
    * Sets the number of spatial dimensions for the compartment.  Must be 0, 1, 2, or 3.
    *
    * @param dimensions Spatial dimensions
    */

   public void setSpatialDimensions (int dimensions) {
      if (spatialDimensions < 0 || spatialDimensions > 3)
         throw new IllegalArgumentException ("Spatial dimensions must be 0, 1, 2, or 3.");
      spatialDimensions = dimensions;
   }

   public void setUnits (String units) {
      if (units != null && !isValidId (units))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.units = units;
   }

   public void setUnits (UnitDefinition units) {
      setUnits (units == null ? null : units.getId ());
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      if (attributes.getIndex ("constant") != -1)
         setConstant (Boolean.valueOf (attributes.getValue ("constant")) == Boolean.TRUE);
      setOutside (attributes.getValue ("outside"));
      if (attributes.getIndex ("spatialDimensions") != -1)
         setSpatialDimensions (Integer.parseInt (attributes.getValue ("spatialDimensions")));
      if (attributes.getIndex ("size") != -1)
         setSize (Double.parseDouble (attributes.getValue ("size")));
      setUnits (attributes.getValue ("units"));
    }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "compartment");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      if (getSpatialDimensions () != 3)
         printer.addAttribute ("spatialDimensions", String.valueOf (getSpatialDimensions ()));
      if (!Double.isNaN (getSize ()))
         printer.addAttribute ("size", String.valueOf (getSize ()));
      if (!getUnits ().equals (getDefaultUnits ()))
         printer.addAttribute ("units", getUnits ());
      printer.addAttribute ("outside", getOutside ());
      if (!isConstant ())
         printer.addAttribute ("constant", "false");
      return printer;
   }

   private String getDefaultUnits () {
      switch (getSpatialDimensions ()) {
         case 0 : return "dimensionless";
         case 1 : return "length";
         case 2 : return "area";
         case 3 : return "volume";
      }
      throw new IllegalStateException ("Spatial dimensions must be 0, 1, 2, or 3.");
   }
}
