package jigcell.sbml2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Stores the note metadata for an SBML node.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 * @author Ranjit Randhawa
 */

public final class Notes {
   private ArrayList notes;

   private final static Pattern checkNoteFormat = Pattern.compile("(?s).*<[sbml:]?notes>.?<p xmlns=\"http://www.w3.org/1999/xhtml\">(.*)</p>.?</[sbml:]?notes>.*");
   private final static Pattern getNotes = Pattern.compile("<notes>(.*)</notes>");
   
   private static String wrapItem (String item) {
      Matcher checkNotes = checkNoteFormat.matcher (item);
      if (checkNotes.matches ())
         return item;
      checkNotes = getNotes.matcher (item);
      if (checkNotes.matches ())
         item = checkNotes.group (1);
      return "<notes><p xmlns=\"http://www.w3.org/1999/xhtml\">" + XMLPrinter.quote(item) + "</p></notes>\n";
   }

   public void add (String note) {
      notes.add (wrapItem (note));
   }

   public void clear () {
      notes.clear();
   }

   public String get (int index) {
      return (String) notes.get (index);
   }

   public int size () {
      return notes.size ();
   }

   public Iterator iterator () {
      return notes.iterator ();
   }

   public void remove (int index) {
      notes.remove (index);
   }

   public void set (int index, String note) {
      notes.set (index, wrapItem (note));
   }

   /**
    * The SBML for this element.
    */

   public String toString () {
      StringBuffer buffer = new StringBuffer  ();
      for (Iterator iterator = notes.iterator (); iterator.hasNext (); ) {
         String note = (String) iterator.next ();
         if (note.split ("\\s*</??notes>\\s*").length != 0)
            buffer.append (note + "\n");
      }
      return buffer.toString ();
   }

   public Notes ( Notes src ) {
      notes = (ArrayList)src.notes.clone();
   }

   public Notes () {
      notes = new ArrayList();
   }
}
