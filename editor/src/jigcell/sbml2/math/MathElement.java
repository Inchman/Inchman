package jigcell.sbml2.math;

/**
 * A class that holds a user or program defined value. This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more
 * details.
 *
 * @author Marc Vass
 */

public class MathElement {
   private String value;

   /**
    * Creates a new instance of MathElement
    */

   public MathElement () {}

   public MathElement (String s) {
      value = s;
   }

   public void setValue (String s) {
      value = s;
   }

   public String toString () {
      return value;
   }
}
