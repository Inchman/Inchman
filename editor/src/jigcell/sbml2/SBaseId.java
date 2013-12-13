package jigcell.sbml2;

import java.util.regex.Pattern;
import org.xml.sax.Attributes;

/**
 * An SBML element that has a name and id.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public abstract class SBaseId extends SBase {
   private final static Pattern PATTERN_ID = Pattern.compile ("\\A[a-zA-Z_][a-zA-Z0-9_]*\\z");

   private String id;
   private String name;

   public final static boolean isValidId (String id) {
      return id != null && PATTERN_ID.matcher (id).matches ();
   }

   public SBaseId () {
      this (null, null);
   }

   public SBaseId (String id, String name) {
      if (id != null)
         setId (id);
      if (name != null)
         setName (name);
   }

   public final String getId () {
      return id;
   }

   public final String getName () {
      return name == null ? id : name;
   }

   public boolean isNameSet () {
      return name != null;
   }

   public boolean isValid (Model model) {
      if (!super.isValid (model))
         return false;
      String id = getId ();
      return this instanceof Event 
         || this instanceof Link 
         || this instanceof Model 
         || id != null;
   }

   public final void setId (String id) {
      if (id != null && !isValidId (id))
         throw new IllegalArgumentException ("Invalid SBML identifier.");
      this.id = id;
   }

   public final void setName (String name) {
      this.name = name;
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      setId (attributes.getValue ("id"));
      setName (attributes.getValue ("name"));
    }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      if (id != null)
         printer.addAttribute ("id", getId ());
      if (this.name != null)
         printer.addAttribute ("name", this.name);
      return printer;
   }
}
