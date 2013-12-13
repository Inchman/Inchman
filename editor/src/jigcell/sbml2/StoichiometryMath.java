package jigcell.sbml2;

/**
 * Special element used when the stoichiometry of a species reference is given in MathML.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public final class StoichiometryMath extends SBase implements MathElement {
   private String math;

   public StoichiometryMath () {
      super ();
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param stoichiometry Stoichiometry math, must not be null
    */

   public StoichiometryMath (StoichiometryMath stoichiometry) {
      this ();
      setMath (stoichiometry.getMath ());
   }

   public String getMath () {
      return math;
   }

   public void setMath (String math) {
      assert math == null || math.startsWith ("<math:math>");
      this.math = math;
   }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "stoichiometryMath");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addText (getMath ());
      return printer;
   }
}
