package jigcell.sbml2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Gives an abbreviated name to a combination of units.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public final class UnitDefinition extends SBaseId {
   private final List units;
   private final SBase unitsElement;

   public static boolean isValidAreaUnit (UnitDefinition areaUnit) {
      if (areaUnit == null)
         return false;
      String id = areaUnit.getId ();
      if (id.equals ("area"))
         return true;
      List units = areaUnit.getUnits ();
      if (units.size () != 1)
         return false;
      Unit unit = (Unit) units.get (0);
      return unit.getKindName ().equals ("metre") && unit.getExponent () == 2;
   }

   public static boolean isValidLengthUnit (UnitDefinition lengthUnit) {
      if (lengthUnit == null)
         return false;
      String id = lengthUnit.getId ();
      if (id.equals ("length"))
         return true;
      List units = lengthUnit.getUnits ();
      if (units.size () != 1)
         return false;
      Unit unit = (Unit) units.get (0);
      return unit.getKindName ().equals ("metre") && unit.getExponent () == 1;
   }

   public static boolean isValidSpatialSizeUnit (UnitDefinition spatialSizeUnit, int spatialDimensions) {
      if (spatialSizeUnit == null)
         return false;
      switch (spatialDimensions) {
         case 0 :
            return isValidSubstanceUnit (spatialSizeUnit);
         case 1 :
            return isValidLengthUnit (spatialSizeUnit);
         case 2 :
            return isValidAreaUnit (spatialSizeUnit);
         case 3 :
            return isValidVolumeUnit (spatialSizeUnit);
      }
      return false;
   }

   public static boolean isValidSubstanceUnit (UnitDefinition substanceUnit) {
      if (substanceUnit == null)
         return false;
      String id = substanceUnit.getId ();
      if (id.equals ("substance") || id.equals ("moles") || id.equals ("item"))
         return true;
      List units = substanceUnit.getUnits ();
      if (units.size () != 1)
         return false;
      Unit unit = (Unit) units.get (0);
      return (unit.getKindName ().equals ("moles") || unit.getKindName ().equals ("item")) && unit.getExponent () == 1;
   }

   public static boolean isValidTimeUnit (UnitDefinition timeUnit) {
      if (timeUnit == null)
         return false;
      String id = timeUnit.getId ();
      if (id.equals ("time") || id.equals ("second"))
         return true;
      List units = timeUnit.getUnits ();
      if (units.size () != 1)
         return false;
      Unit unit = (Unit) units.get (0);
      return unit.getKindName ().equals ("second") && unit.getExponent () == 1;
   }

   public static boolean isValidVolumeUnit (UnitDefinition volumeUnit) {
      if (volumeUnit == null)
         return false;
      String id = volumeUnit.getId ();
      if (id.equals ("volume"))
         return true;
      List units = volumeUnit.getUnits ();
      if (units.size () != 1)
         return false;
      Unit unit = (Unit) units.get (0);
      return unit.getKindName ().equals ("metre") && unit.getExponent () == 3 ||
         unit.getKindName ().equals ("litre") && unit.getExponent () == 1;
   }

   public UnitDefinition () {
      this (null, null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param definition Unit definition, must not be null
    */

   public UnitDefinition (UnitDefinition definition) {
      this (definition.getId (), definition.isNameSet () ? definition.getName () : null);
      for (Iterator iterator = definition.getUnits ().iterator (); iterator.hasNext (); ) {
         Unit unit = (Unit) iterator.next ();
         if (unit != null)
            addUnit (new Unit (unit));
      }
   }

   public UnitDefinition (String id, String name) {
      super (id, name);
      unitsElement = new SBase ();
      units = new ArrayList ();
   }

   public void addUnit (Unit unit) {
      if (unit == null)
         throw new IllegalArgumentException ();
      units.add (unit);
   }

   public SBase getUnitsElement () {
      return unitsElement;
   }

   /**
    * List[Unit] of units combined in this definition.
    */

   public List getUnits () {
      return units;
   }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "unitDefinition");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addElementList (getUnitsElement (), "listOfUnits", getUnits ());
      return printer;
   }
}
