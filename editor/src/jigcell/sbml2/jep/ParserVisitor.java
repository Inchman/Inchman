/* Generated By:JJTree: Do not edit this line. ParserVisitor.java */
package jigcell.sbml2.jep;

public interface ParserVisitor {
   Object visit (SimpleNode node, Object data);

   Object visit (ASTStart node, Object data);

   Object visit (ASTFunNode node, Object data);

   Object visit (ASTVarNode node, Object data);

   Object visit (ASTIfThen node, Object data);

   Object visit (ASTElseIfThen node, Object data);

   Object visit (ASTElse node, Object data);

   Object visit (ASTPiecewise node, Object data);

   Object visit (ASTConstant node, Object data);
}
