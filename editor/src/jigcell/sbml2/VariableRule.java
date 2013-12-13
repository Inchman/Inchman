package jigcell.sbml2;

import org.xml.sax.Attributes;

/**
 * A rule whose left-hand side is a variable.  This is not an official type in SBML.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public abstract class VariableRule extends Rule {
   private String variable;

   public final String getVariable () {
      return variable;
   }

   public final SBaseId getVariable (Model model) {
      return (SBaseId) model.findElementWithId (variable, Model.COMPARTMENT | Model.PARAMETER | Model.SPECIES);
   }

   public boolean isValid (Model model) {
      return super.isValid (model) && getVariable (model) != null;
   }

   public void setVariable (String variable) {
      if (variable != null && !SBaseId.isValidId (variable))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.variable = variable;
   }

   public final void setVariable (Compartment variable) {
      setVariable (variable == null ? null : variable.getId ());
   }

   public final void setVariable (Parameter variable) {
      setVariable (variable == null ? null : variable.getId ());
   }

   public final void setVariable (Species variable) {
      setVariable (variable == null ? null : variable.getId ());
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      setVariable (attributes.getValue ("variable"));
    }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addAttribute ("variable", getVariable ());
      return printer;
   }
}
