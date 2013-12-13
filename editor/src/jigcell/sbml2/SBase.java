package jigcell.sbml2;

import java.util.Iterator;
import java.util.List;
import org.xml.sax.Attributes;

/**
 * Contains the notes, annotations, and metadata for an SBML element.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public class SBase extends XMLElement {
   protected Annotations annotations;
   protected Notes notes;
   private String metaid;
   private String rdf;

   final static SBaseId searchListForId (List elements, String id) {
      if (id == null)
         return null;
      for (Iterator iterator = elements.iterator (); iterator.hasNext (); ) {
         SBaseId element = (SBaseId) iterator.next ();
         if (id.equals (element.getId ()))
            return element;
      }
      return null;
   }

   final static SBaseId searchListForName (List elements, String name) {
      if (name == null)
         return null;
      for (Iterator iterator = elements.iterator (); iterator.hasNext (); ) {
         SBaseId element = (SBaseId) iterator.next ();
         if (name.equals (element.getName ()))
            return element;
      }
      return null;
   }

   public final Annotations getAnnotations () {
      return annotations;
   }

   public final String getMetaid () {
      return metaid;
   }

   public final Notes getNotes () {
      return notes;
   }

   public final String getRDF () {
      return rdf;
   }

   public boolean isValid (Model model) {
      return true;
   }

   public void setMetaid (String metaid) {
      this.metaid = metaid;
   }

   public void setRDF (String rdf) {
      this.rdf = rdf;
   }

   public void setRDF (String rdf, String metaid) {
      setRDF (rdf);
      setMetaid (metaid);
   }

   protected SBase () {
      annotations = new Annotations ();
      notes = new Notes ();
   }

   protected void parse (Attributes attributes) {
      setMetaid (attributes.getValue ("metaid"));
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      String notesText = getNotes ().toString ().trim ();
      if (notesText.length () > 0)
         printer.addText (notesText);
      String annotationsText = getAnnotations ().toString ().trim ();
      if (annotationsText.length () > 0)
         printer.addText (annotationsText);
      String rdfText = getRDF ();
      if (rdfText != null) {
         rdfText = rdfText.trim ();
         if (rdfText.length () > 0)
            printer.addText (rdfText);
      }
      printer.addAttribute ("metaid", getMetaid ());
      return printer;
   }
}
