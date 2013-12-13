package jigcell.sbml2;

/**
 * An element that contains MathML.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public interface MathElement {
   String getMath ();

   void setMath (String math);
}

