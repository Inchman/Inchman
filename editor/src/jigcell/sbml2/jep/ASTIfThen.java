/* Generated By:JJTree: Do not edit this line. ASTIfThen.java */
package jigcell.sbml2.jep;

public class ASTIfThen extends SimpleNode {
   public ASTIfThen (int id) {
      super (id);
   }

   public ASTIfThen (Parser p, int id) {
      super (p, id);
   }

   /**
    * Accept the visitor.
    */

   public Object jjtAccept (ParserVisitor visitor, Object data) {
      return visitor.visit (this, data);
   }
}
