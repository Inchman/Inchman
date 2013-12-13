package jigcell.sbml2;

/**
 * A rule dx/dt=f(W) where x is a variable and f is an arbitrary function.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * <p>
 *
 * @author Nicholas Allen
 */

public class RateRule extends VariableRule {
   public RateRule () {
      super ();
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param rule Rule, must not be null
    */

   public RateRule (RateRule rule) {
      this ();
      setMath (rule.getMath ());
      setVariable (rule.getVariable ());
   }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "rateRule");
   }
}
