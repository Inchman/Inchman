package jigcell.sbml2;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.xml.sax.Attributes;

/**
 * Units of measurement for quantities in an SBML model.  A unit is defined by unit = (multiplier * 10^scale * kind^exponent) + offset.
 *
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 *
 * @author Nicholas Allen
 */

public final class Unit extends SBase {
   private final static Map baseUnits = new HashMap ();

   private double multiplier;
   private double offset;
   private int exponent;
   private int scale;
   private String name;
   private Unit kind;

   static {
      addBaseUnit ("ampere");
      addBaseUnit ("becquerel");
      addBaseUnit ("candela");
      addBaseUnit ("Celsius");
      addBaseUnit ("coulomb");
      addBaseUnit ("dimensionless");
      addBaseUnit ("farad");
      addBaseUnit ("gram");
      addBaseUnit ("gray");
      addBaseUnit ("henry");
      addBaseUnit ("hertz");
      addBaseUnit ("item");
      addBaseUnit ("joule");
      addBaseUnit ("katal");
      addBaseUnit ("kelvin");
      addBaseUnit ("kilogram");
      addBaseUnit ("litre");
      addBaseUnit ("lumen");
      addBaseUnit ("lux");
      addBaseUnit ("metre");
      addBaseUnit ("mole");
      addBaseUnit ("newton");
      addBaseUnit ("ohm");
      addBaseUnit ("pascal");
      addBaseUnit ("radian");
      addBaseUnit ("second");
      addBaseUnit ("siemens");
      addBaseUnit ("sievert");
      addBaseUnit ("steradian");
      addBaseUnit ("tesla");
      addBaseUnit ("volt");
      addBaseUnit ("watt");
      addBaseUnit ("weber");
   }

   public static Unit findBaseUnit (String name) {
      return (Unit) baseUnits.get (name);
   }

   /**
    * Set[Unit] of the SBML base units.
    */

   public static Set getBaseUnits () {
      return new HashSet (baseUnits.values ());
   }

   private static void addBaseUnit (String name) {
      baseUnits.put (name, new Unit (name));
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param unit Unit, must not be null
    */

   public Unit (Unit unit) {
      this (unit.getKind (), unit.getMultiplier (), unit.getScale (), unit.getExponent (), unit.getOffset ());
   }

   /**
    * Creates a new unit.  unit = (multiplier * 10^scale * kind^exponent) + offset
    */

   public Unit (Unit kind, double multiplier, int scale, int exponent, double offset) {
      this ();
      setKind (kind);
      setMultiplier (multiplier);
      setScale (scale);
      setExponent (exponent);
      setOffset (offset);
   }

   public boolean equals (Object o) {
      if (!(o instanceof Unit))
         return false;
      if (o == this)
         return true;
      Unit unit = (Unit) o;
      return getKindName ().equals (unit.getKindName ()) && getMultiplier () == unit.getMultiplier () && getScale () == unit.getScale () &&
         getExponent () == unit.getExponent () && getOffset () == unit.getOffset ();
   }

   public int getExponent () {
      return exponent;
   }

   public Unit getKind () {
      return isBaseUnit () ? this : kind;
   }

   public String getKindName () {
      return isBaseUnit () ? name : kind.getKindName ();
   }

   public double getMultiplier () {
      return multiplier;
   }

   public double getOffset () {
      return offset;
   }

   public int getScale () {
      return scale;
   }

   public int hashCode () {
      int sum = getKindName ().hashCode ();
      long multiplierBits = Double.doubleToLongBits (getMultiplier ());
      sum += (int) (multiplierBits ^ (multiplierBits >>> 32));
      sum += getScale ();
      sum += getExponent () << 16;
      long offsetBits = Double.doubleToLongBits (getOffset ());
      return sum + (int) (offsetBits ^ (offsetBits >>> 32));
   }

   public boolean isBaseUnit () {
      return kind == null;
   }

   public void setExponent (int exponent) {
      if (isBaseUnit ())
         throw new UnsupportedOperationException ("Base units cannot be modified.");
      this.exponent = exponent;
   }

   public void setKind (Unit kind) {
      if (isBaseUnit ())
         throw new UnsupportedOperationException ("Base units cannot be modified.");
      if (kind == null)
         throw new IllegalArgumentException ("Must specify a unit kind");
      if (!kind.isBaseUnit ())
         throw new IllegalArgumentException ("Units cannot be derived from a user defined unit");
      this.kind = kind;
   }

   public void setMultiplier (double multiplier) {
      if (isBaseUnit ())
         throw new UnsupportedOperationException ("Base units cannot be modified.");
      this.multiplier = multiplier;
   }

   public void setOffset (double offset) {
      if (isBaseUnit ())
         throw new UnsupportedOperationException ("Base units cannot be modified.");
      this.offset = offset;
   }

   public void setScale (int scale) {
      if (isBaseUnit ())
         throw new UnsupportedOperationException ("Base units cannot be modified.");
      this.scale = scale;
   }

   Unit () {
      super ();
      this.kind = this;
      setMultiplier (1.0);
      setScale (0);
      setExponent (1);
      setOffset (0.0);
   }

   /**
    * Creates a new base unit.
    *
    * @param name Unit name
    */

   Unit (String name) {
      this ();
      this.name = name;
      kind = null;
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      if (attributes.getIndex ("exponent") != -1)
         setExponent (Integer.parseInt (attributes.getValue ("exponent")));
      setKind (findBaseUnit (attributes.getValue ("kind")));
      if (attributes.getIndex ("multiplier") != -1)
         setMultiplier (Double.parseDouble (attributes.getValue ("multiplier")));
      if (attributes.getIndex ("offset") != -1)
         setOffset (Double.parseDouble (attributes.getValue ("offset")));
      if (attributes.getIndex ("scale") != -1)
         setScale (Integer.parseInt (attributes.getValue ("scale")));
    }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "unit");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addAttribute ("kind", getKindName ());
      if (getExponent () != 1)
         printer.addAttribute ("exponent", String.valueOf (getExponent ()));
      if (getScale () != 0)
         printer.addAttribute ("scale", String.valueOf (getScale ()));
      if (getMultiplier () != 1.0)
         printer.addAttribute ("multiplier", String.valueOf (getMultiplier ()));
      if (getOffset () != 0.0)
         printer.addAttribute ("offset", String.valueOf (getOffset ()));
      return printer;
   }
}
