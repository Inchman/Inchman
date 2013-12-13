package jigcell.sbml2;

/**
 * A rule constrains a model variable according to the given equation.  The three types of rules are algebraic, assignment, and rate.  Algebraic
 * rules define the constrained variable implicitly; assignment and rate rules have an explicit variable.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public abstract class Rule extends SBase implements MathElement {
   private String math;

   public String getMath () {
      return math;
   }

   public void setMath (String math) {
      assert math == null || math.startsWith ("<math:math>");
      this.math = math;
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addText (getMath ());
      return printer;
   }
}
