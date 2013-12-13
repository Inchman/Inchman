/**
 * JEP - Java Expression Parser JEP is a Java package for parsing and evaluating mathematical expressions. It currently supports user defined
 * variables, constant, and functions. A number of common mathematical functions and constants are included. Author: Nathan Funk Copyright (C)
 * 2000 Nathan Funk JEP is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. JEP is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with
 * JEP; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

package jigcell.sbml2.jep.function;

import java.util.Stack;
import jigcell.sbml2.FunctionDefinition;

/**
 * Function classes extend this class. It is an implementation of the PostfixMathCommandI interface.
 *
 * <p>
 * It includes a numberOfParameters member, that is checked when parsing the expression. This member should be initialized to an appropriate
 * value for all classes extending this class. If an arbitrary number of parameters should be allowed, initialize this member to -1.
 * </p>
 */

public class PostfixMathCommand extends FunctionDefinition implements PostfixMathCommandI {

   /**
    * Number of parameters to be used for the next run() invocation. Applies only if the required umber of parameters is variable
    * (numberOfParameters = -1).
    */

   protected int curNumberOfParameters;

   /**
    * Number of parameters a the function requires. Initialize this value to -1 if any number of parameters should be allowed.
    */

   protected int numberOfParameters;

   /**
    * Creates a new PostfixMathCommand class.
    */

   public PostfixMathCommand (int numParameters) {
      super ();
      this.numberOfParameters = numParameters;
      curNumberOfParameters = numParameters;
   }

   /**
    * Return the required number of parameters.
    */

   public int getNumberOfParameters () {
      return numberOfParameters;
   }

   /**
    * Throws an exception because this method should never be called under normal circumstances. Each function should use it's own run() method
    * for evaluating the function. This includes popping off the parameters from the stack, and pushing the result back on the stack.
    */

   public void run (Stack s) throws Exception {
      throw new Exception ("run() method of PostfixMathCommand called");
   }

   /**
    * Sets the number of current number of parameters used in the next call of run(). This method is only called when the reqNumberOfParameters
    * is -1.
    */

   public void setCurNumberOfParameters (int n) {
      curNumberOfParameters = n;
   }

   /**
    * Convenience method for use with the parser to set the parameters if they have changed.
    */

   public void setNumberOfParameters (int n) {
      numberOfParameters = n;
      curNumberOfParameters = n;
   }

   /**
    * Check whether the stack is not null, throw an Exception if it is.
    */

   protected void checkStack (Stack inStack) throws Exception {

      /* Check if stack is null */
      if (null == inStack)
         throw new Exception ("Stack argument null");
   }
}
