package jigcell.sbml2;

/**
 * A rule 0=f(W) where f is an arbitrary function and a suitable W must be found.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class AlgebraicRule extends Rule {
   public AlgebraicRule () {
      super ();
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param rule Rule, must not be null
    */

   public AlgebraicRule (AlgebraicRule rule) {
      this ();
      setMath (rule.getMath ());
   }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "algebraicRule");
   }
}
