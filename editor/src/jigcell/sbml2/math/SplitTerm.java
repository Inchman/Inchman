package jigcell.sbml2.math;

public class SplitTerm {
   /* This will be null if the coefficient could not be calculated as a number */
   public Double coefficient    = null;
   /* This is the tree with id removed (if id was found, otherwise it's unchanged with id not present). */
   public Node   termWithoutVar = null;
   /* This is the id removed from the tree */
   public String id             = null;
}
