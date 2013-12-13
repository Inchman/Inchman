package jigcell.sbml2.math;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.xml.sax.Attributes;

/**
 * Represents a node in a parse tree based on XML (in particular MathML). Node is publicly a non-mutable class.  Therefore, it is safe to give
 * many Objects a references to the same Node.
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * @author Marc Vass
 */

public class Node {
   protected List children;
   protected Attributes attributes;
   protected Node parent;
   protected String localName;
   protected String qName;
   protected String uri;
   protected StringBuffer value;
   protected boolean error; // if there is an error in the parse tree at this
                            // Node, then error is set to true.

   public Node (Node oldNode) {
      attributes = oldNode.attributes;
      children = new ArrayList (oldNode.children);
      localName = oldNode.localName;
      parent = oldNode.parent;
      qName = oldNode.qName;
      uri = oldNode.uri;
      value = new StringBuffer (oldNode.value.toString ());
   }

   public Node (Node oldNode, String value, Node children[]) {
      this(oldNode);
      this.children = Arrays.asList(children);
      this.value = new StringBuffer(value);
   }

   public Node (Node parent, String uri, String localName, String qName, Attributes attributes) {
      this.parent = parent;
      this.attributes = attributes;
      this.localName = localName;
      this.qName = qName;
      this.uri = uri;
      value = new StringBuffer ("");
      children = new ArrayList ();
   }

   // add by Pengyuan
   public Node(Node parent, String qName, String value) {
	   this(parent, null, null, qName, null);
	   this.value = new StringBuffer(value);
   }
   
   // add by Pengyuan
   public Node(Node parent, Node oldNode) {
	   this(oldNode);
	   this.parent = parent;
   }
   
   public Node (String uri, String localName, String qName, Attributes attributes) {
      this (null, uri, localName, qName, attributes);
   }

   public Node (String qName, String value, Node children []) {
      this (null, null, null, qName, null);
      this.value = new StringBuffer (value);
      for (int i = 0; i < children.length; i++)
         addChild (children [i]);
   }

   public Node (String qName, String value, Node child1, Node child2) {
      this (null, null, null, qName, null);
      this.value = new StringBuffer (value);
      addChild (child1);
      addChild (child2);
   }

   public Node (String qName, String value, Node child) {
      this (null, null, null, qName, null);
      this.value = new StringBuffer (value);
      addChild (child);
   }

   public Node (String qName, String value) {
      this (null, null, null, qName, null);
      this.value = new StringBuffer (value);
   }

   public Node () {
      this (null, "", "", "", null);
   }

   // Getters follow
   public Attributes getAttributes () {
      return attributes;
   }

   public Node getChild (int index) {
      return (Node) children.get (index);
   }

   public Node [] getChildren () {
      return (Node []) children.toArray (new Node [0]);
   }
   
   public List getListOfChildren () {
      return children;
   }

   public String getLocalName () {
      return localName;
   }

   public int getNumChildren () {
      return children.size ();
   }

   public Node getParent () {
      return parent;
   }

   public String getQName () {
      return qName;
   }

   // Returns the qName with out a prepended "math:"
   public String getSimpleName () {
      if (qName.startsWith ("math:"))
         return qName.substring (5);
      return qName;
   }

   public String getUri () {
      return uri;
   }

   public String getValue () {
      return value.toString ();
   }

   public void addChild (Node child) {
      child.setParent (this);
      children.add (child);
   }

   public void addChild (int index, Node child) {
      child.setParent (this);
      children.add (index,child);
   }

   protected void appendToValue (String s) {
      value.append (s);
   }

   // End Getters
   protected void removeAllChildren () {
      children = new ArrayList ();
   }

   protected Node removeChild (int index) {
      return (Node) children.remove (index);
   }

   protected int removeChild( Node child ) {
      for ( int i = 0; i < children.size(); i++ ) {
         if ( children.get(i) == child ) {
            children.remove(i);
            return i;
         }
      }
      return -1;
   }

   // Setters follow
   protected void setAttributes (Attributes attributes) {
      this.attributes = attributes;
   }

   public Node setChild( int index, Node child ) {
      child.setParent(this);
      return (Node)children.set(index,child);
   }

   public void setError( boolean error ) {
      this.error = error;
   }

   protected void setLocalName (String localName) {
      this.localName = localName;
   }

   protected void setParent (Node parent) {
      this.parent = parent;
   }

   protected void setQName (String qName) {
      this.qName = qName;
   }

   protected void setUri (String uri) {
      this.uri = uri;
   }

   protected void setValue (String value) {
      this.value = new StringBuffer (value);
   }

   // End Setters

   public boolean hasError() {
      return error;
   }

   public String toString() {
      return toString("");
   }

   private String toString( String indent ) {
      String thisNodeS = indent+qName+": "+value+"\n";
      indent += "  ";
      String childrenS = "";
      for ( int i=0; i < children.size(); i++ ) {
         childrenS += ((Node)children.get(i)).toString(indent);
      }
      return thisNodeS+childrenS;
   }

   public boolean equals( Node other ) {
      boolean equal = true;
      String thisSimpleQName = this.qName.toString().replaceFirst("^math:","");
      String otherSimpleQName =
          other.qName.toString().replaceFirst("^math:","");
      equal &= thisSimpleQName.equalsIgnoreCase(otherSimpleQName);
      equal &= this.value.toString().trim().equals(
         other.value.toString().trim());
      if ( this.children.size() != other.children.size() ) {
          return false;
      }
      for ( int i = 0; i < children.size(); i++ ) {
          equal &= ((Node)this.children.get(i)).equals(
              (Node)other.children.get(i));
      }
      return equal;
   }

}
