package jigcell.sbml2;

/**
 * Associates an identifier with a function definition.  The identifier can then be used in a MathML apply element.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class FunctionDefinition extends SBaseId implements MathElement {
   private String math;

   public FunctionDefinition () {
      this (null, null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param definition Function definition, must not be null
    */

   public FunctionDefinition (FunctionDefinition definition) {
      this (definition.getId (), definition.isNameSet () ? definition.getName () : null);
      setMath (definition.getMath ());
   }

   public FunctionDefinition (String id, String name) {
      super (id, name);
   }

   public String getMath () {
      return math;
   }

   public void setMath (String math) {
      assert math == null || math.startsWith ("<math:math>");
      this.math = math;
   }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "functionDefinition");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addText (getMath ());
      return printer;
   }
}
