package jigcell.sbml2;

import org.xml.sax.Attributes;

/**
 * A variable assignment that takes place when an event is executed.  The variable is a compartment, species, or parameter of the model.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class EventAssignment extends SBase implements MathElement {
   private String math;
   private String variable;

   public EventAssignment () {
      super ();
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param assignment Event assignment, must not be null
    */

   public EventAssignment (EventAssignment assignment) {
      this ();
      setMath (assignment.getMath ());
      setVariable (assignment.getVariable ());
   }

   public String getMath () {
      return math;
   }

   public String getVariable () {
      return variable;
   }

   public SBaseId getVariable (Model model) {
      return (SBaseId) model.findElementWithId (variable, Model.COMPARTMENT | Model.PARAMETER | Model.SPECIES);
   }

   public boolean isValid (Model model) {
      return super.isValid (model) && getVariable (model) != null;
   }

   public void setMath (String math) {
      assert math == null || math.startsWith ("<math:math>");
      this.math = math;
   }

   public void setVariable (String variable) {
      if (variable != null && !SBaseId.isValidId (variable))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.variable = variable;
   }

   public void setVariable (Compartment variable) {
      setVariable (variable == null ? null : variable.getId ());
   }

   public void setVariable (Parameter variable) {
      setVariable (variable == null ? null : variable.getId ());
   }

   public void setVariable (Species variable) {
      setVariable (variable == null ? null : variable.getId ());
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      setVariable (attributes.getValue ("variable"));
    }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "eventAssignment");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addAttribute ("variable", getVariable ());
      printer.addText (getMath ());
      return printer;
   }
}
