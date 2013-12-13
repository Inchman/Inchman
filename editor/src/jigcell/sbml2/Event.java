package jigcell.sbml2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.xml.sax.Attributes;

/**
 * Instantaneous discontinuous state changes to the model that are applied when a criterion is met.  The state changes are given by the list
 * of event assignments.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class Event extends SBaseId {
   protected List eventAssignment;
   private final SBase assignmentsElement;
   private final SBase delayElement;
   private final SBase triggerElement;
   private String delay;
   private String timeUnits;
   private String trigger;

   public Event () {
      this (null, null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param event Event, must not be null
    */

   public Event (Event event) {
      this (event.getId (), event.isNameSet () ? event.getName () : null);
      setDelay (event.getDelay ());
      setTimeUnits (event.getTimeUnits ());
      setTrigger (event.getTrigger ());
      for (Iterator iterator = event.getEventAssignment ().iterator (); iterator.hasNext (); ) {
         EventAssignment assignment = (EventAssignment) iterator.next ();
         if (assignment != null)
            addEventAssignment (new EventAssignment (assignment));
      }
   }

   public Event (String id, String name) {
      super (id, name);
      assignmentsElement = new SBase ();
      delayElement = new SBase ();
      triggerElement = new SBase ();
      eventAssignment = new ArrayList ();
      setDelay (null);
   }

   public void addEventAssignment (EventAssignment assignment) {
      if (assignment == null)
         throw new IllegalArgumentException ();
      eventAssignment.add (assignment);
   }

   public SBase getAssignmentsElement () {
      return assignmentsElement;
   }

   public String getDelay () {
      return delay;
   }

   public SBase getDelayElement () {
      return delayElement;
   }

   public List getEventAssignment () {
      return eventAssignment;
   }

   public String getTimeUnits () {
      return timeUnits == null ? "time" : timeUnits;
   }

   public UnitDefinition getTimeUnits (Model model) {
      return (UnitDefinition) model.findElementWithId (getTimeUnits (), Model.UNITDEFINITION);
   }

   public String getTrigger () {
      return trigger;
   }

   public SBase getTriggerElement () {
      return triggerElement;
   }

   public boolean isValid (Model model) {
      return super.isValid (model) && UnitDefinition.isValidTimeUnit (getTimeUnits (model));
   }

   /**
    * Sets the length of time after the event has fired that the event is executed.
    *
    * @param delay New value of property delay.
    */

   public void setDelay (String delay) {
      if (delay == null)
         delay = "<math:math><math:cn>0</math:cn></math:math>";
      assert delay.startsWith ("<math:math>");
      this.delay = delay;
   }

   /**
    * Sets the units of time that apply to the delay field.
    *
    * @param timeUnits New value of property timeUnits.
    */

   public void setTimeUnits (String timeUnits) {
      if (timeUnits != null && !isValidId (timeUnits))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.timeUnits = timeUnits;
   }

   public void setTimeUnits (UnitDefinition units) {
      setTimeUnits (units == null ? null : units.getId ());
   }

   /**
    * Sets the MathML boolean expression that defines when an event is fired on the transition from false to true.
    *
    * @param trigger New value of property trigger.
    */

   public void setTrigger (String trigger) {
      assert trigger == null || trigger.startsWith ("<math:math>");
      this.trigger = trigger;
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      setTimeUnits (attributes.getValue ("timeUnits"));
    }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "event");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      if (!getTimeUnits ().equals ("time"))
         printer.addAttribute ("timeUnits", getTimeUnits ());
      printer.addCustomElement (getTriggerElement (), "trigger", getTrigger ());
      printer.addCustomElement (getDelayElement (), "delay", getDelay ());
      printer.addElementList (getAssignmentsElement (), "listOfEventAssignments", getEventAssignment ());
      return printer;
   }
}
