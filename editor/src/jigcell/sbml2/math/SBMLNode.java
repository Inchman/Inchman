package jigcell.sbml2.math;

import java.util.Vector;

/**
 * A class representing a node in a parse tree. This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more
 * details.
 *
 * @author Marc Vass
 * @author Nicholas Allen
 */

public class SBMLNode {
   public MathElement element;
   public SBMLNode parent;
   public Vector children = new Vector ();

   /**
    * Creates a new instance of SBMLNode
    */

   public SBMLNode (SBMLNode parent) {
      this.parent = parent;
   }

   /**
    * @param parent
    * @param m
    */

   public SBMLNode (SBMLNode parent, MathElement m) {
      this.parent = parent;
      this.element = m;
   }

   /**
    * @param m
    */

   public SBMLNode (MathElement m) {
      this.parent = null;
      this.element = m;
   }

   /**
    * @param n
    */

   public void addChild (SBMLNode n) {
      children.add (n);
   }

   /**
    * Getter for property element.
    *
    * @return Value of property element.
    */

   public MathElement getElement () {
      return element;
   }

   /**
    * @return
    */

   public SBMLNode getParent () {
      return parent;
   }

   /**
    * @param m
    */

   public void setElement (MathElement m) {
      element = m;
   }
}
