package jigcell.sbml2;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.Writer;
import org.xml.sax.Attributes;

/**
 * Represents the SBML Level 2 document that contains the model.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public final class SBMLLevel2Document extends SBase {
   private Model model;

   public static SBMLLevel2Document readDocument (InputStream stream) throws IOException {
      return SBMLLevel2Reader.read (new InputStreamReader (stream));
   }

   public static SBMLLevel2Document readDocument (Reader reader) throws IOException {
      return SBMLLevel2Reader.read (reader);
   }

   public static SBMLLevel2Document readDocument (Reader reader, Model model) throws IOException {
      return SBMLLevel2Reader.read (reader,model);
   }

   public static SBMLLevel2Document readDocument (String fileName) throws IOException {
      return SBMLLevel2Reader.read (new BufferedReader (new FileReader (fileName)));
   }

   public SBMLLevel2Document () {
      this (null);
   }

   public SBMLLevel2Document (Model model) {
      this.model = model;
   }

   public int getLevel () {
      return 2;
   }

   public Model getModel () {
      return model;
   }

   public int getVersion () {
      return 1;
   }

   public boolean isValid (Model model) {
      return super.isValid (model) && getModel () != null;
   }

   public void setLevel (int level) {
      if (level != 2)
         throw new IllegalArgumentException ("Level is required to be 2.");
   }

   public void setModel (Model model) {
      this.model = model;
   }

   public void setVersion (int version) {
      //if (version != 1) -- allow newer versions (Aidan, Sep 2010)
      //   throw new IllegalArgumentException ("Version is required to be 1.");
   }

   public void writeDocument (OutputStream stream) throws IOException {
      writeDocument (new OutputStreamWriter (stream));
   }

   /* removed to gain compatibility with Google App Engine
   public void writeDocument (String fileName) throws IOException {
      writeDocument (new BufferedWriter (new FileWriter (fileName)));
   }
   */

   /**
    * Writes an SBML document and closes the writer.
    */

   public void writeDocument (Writer writer) throws IOException {
      writer.write ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
      print (null).write (writer);
      writer.close ();
   }

   protected void parse (Attributes attributes) {
      super.parse (attributes);
      setLevel (Integer.parseInt (attributes.getValue ("level")));
      setVersion (Integer.parseInt (attributes.getValue ("version")));
    }

   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "sbml");
   }

   protected XMLPrinter print (XMLPrinter parent, String name) {
      XMLPrinter printer = super.print (parent, name);
      printer.addAttribute ("xmlns", "http://www.sbml.org/sbml/level2");
      printer.addAttribute ("xmlns:html", "http://www.w3.org/1999/xhtml");
      printer.addAttribute ("xmlns:jigcell", "http://www.sbml.org/2001/ns/jigcell");
      printer.addAttribute ("xmlns:math", "http://www.w3.org/1998/Math/MathML");
      printer.addAttribute ("xmlns:rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#");
      printer.addAttribute ("xmlns:sbml", "http://www.sbml.org/sbml/level2");
      printer.addAttribute ("xmlns:xlink", "http://www.w3.org/1999/xlink");
      printer.addAttribute ("level", String.valueOf (getLevel ()));
      printer.addAttribute ("version", String.valueOf (getVersion ()));
      printer.addElement (getModel ());
      return printer;
   }
}
