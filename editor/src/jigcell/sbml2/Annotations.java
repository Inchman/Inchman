package jigcell.sbml2;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * Stores the annotation metadata for an SBML node.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public final class Annotations {
   private ArrayList annotations;
   
   // Begin Ted
   private double convFactorValue = 0.0;
   private String convFactorUnit = "";
   private String strValue = "Conversion_Factor_Value=";
   private String strUnit = "Conversion_Factor_Unit=";
   
   private String defTimeUnit = "";
   private String defVolUnit  = "";
   private String defConUnit  = "";
   private String strDefTime = "Defalut_Time=";
   private String strDefVol  = "Defalut_Vol=";
   private String strDefCon  = "Defalut_Con=";
   // End Ted

   private static String wrapItem (String item) {
      return item.startsWith ("<annotation") || item.startsWith ("<sbml:annotation")  ? item : "<annotation>" + item + "</annotation>";
   }
   
   // Begin Ted
   public double getConvFactorValue() {
	   convFactorValue = getAnnoConValue();
	   return convFactorValue;
   }
   
   public void setConvFactorValue(double p) {
	   convFactorValue = p;
   }
   
   public String getConvFactorUnit() {
	   convFactorUnit = getAnnoConUnit();
	   return convFactorUnit;
   }
   
   public void setConvFactorUnit(String p) {
	   convFactorUnit = p;
   }
   
   public String getDefTimeUnit() {
	   defTimeUnit = getAnnoDefTimeUnit();
	   return defTimeUnit;
   }
   
   public void setDefTimeUnit(String p) {
	   defTimeUnit = p;
   }
   
   public String getDefVolUnit() {
	   defVolUnit = getAnnoDefVolUnit();
	   return defVolUnit;
   }
   
   public void setDefVolUnit(String p) {
	   defVolUnit = p;
   }
   
   public String getDefConUnit() {
	   defConUnit = getAnnoDefConUnit();
	   return defConUnit;
   }
   
   public void setDefConUnit(String p) {
	   defConUnit = p;
   }
   // End Ted

   public void add (String annotation) {
      annotations.add (wrapItem (annotation));
   }

   public void add (int index, String annotation) {
      annotations.add (index, wrapItem (annotation));
   }

   public void clear () {
      annotations.clear();
   }

   public String get (int index) {
      return (String) annotations.get (index);
   }

   public int size () {
      return annotations.size ();
   }

   public Iterator iterator () {
      return annotations.iterator ();
   }

   public void remove (int index) {
      annotations.remove (index);
   }

   public void set (int index, String annotation) {
      annotations.set (index, wrapItem (annotation));
   }

   /**
    * The SBML for this element.
    */

   public String toString () {
      StringBuffer buffer = new StringBuffer  ();
      for (Iterator iterator = iterator (); iterator.hasNext (); )
         buffer.append ((String) iterator.next () + "\n");
      return buffer.toString ();
   }

   /**
    * Creates a new Annotations object that copies all the attributes
    * (protecting the objects of this Annotations object).
    */
   public Annotations ( Annotations src ) {
      annotations = (ArrayList)src.annotations.clone();
   }

   public Annotations () {
      annotations = new ArrayList();
   }
   
   // Begin Ted
   public double getAnnoConValue()
   {
	   String annoStr = toString().trim();
	   if (annoStr == "") {
		   return 0.0;
	   }
	   else {
		   int startpoint = annoStr.indexOf(strValue);
		   int endpoint   = annoStr.indexOf(",", startpoint);
		   if(startpoint == -1)
			   return 0.0;
		   String conValue = annoStr.substring(startpoint + strValue.length(), endpoint);
		   return Double.parseDouble(conValue);
	   }
   }
   
   public String getAnnoConUnit()
   {
	   String annoStr = toString().trim();
	   if (annoStr == "") {
		   return "";
	   }
	   else {
		   int startpoint = annoStr.indexOf(strUnit);
		   int endpoint   = annoStr.indexOf("\"", startpoint);
		   if(startpoint == -1)
			   return "";
		   String conUnit = annoStr.substring(startpoint + strUnit.length(), endpoint);
		   return conUnit;
	   }
   }
   
   public String getAnnoDefTimeUnit()
   {
	   String annoStr = toString().trim();
	   if (annoStr == "") {
		   return "";
	   }
	   else {
		   int startpoint = annoStr.indexOf(strDefTime);
		   int endpoint   = annoStr.indexOf(",", startpoint);
		   if(startpoint == -1)
			   return "";
		   String timeUnit = annoStr.substring(startpoint + strDefTime.length(), endpoint);
		   return timeUnit;
	   }
   }
   
   public String getAnnoDefVolUnit()
   {
	   String annoStr = toString().trim();
	   if (annoStr == "") {
		   return "";
	   }
	   else {
		   int startpoint = annoStr.indexOf(strDefVol);
		   int endpoint   = annoStr.indexOf(",", startpoint);
		   if(startpoint == -1)
			   return "";
		   String volUnit = annoStr.substring(startpoint + strDefVol.length(), endpoint);
		   return volUnit;
	   }
   }
   
   public String getAnnoDefConUnit()
   {
	   String annoStr = toString().trim();
	   if (annoStr == "") {
		   return "";
	   }
	   else {
		   int startpoint = annoStr.indexOf(strDefCon);
		   int endpoint   = annoStr.indexOf(",", startpoint);
		   if(startpoint == -1)
			   return "";
		   String conUnit = annoStr.substring(startpoint + strDefCon.length(), endpoint);
		   return conUnit;
	   }
   }
   // End Ted
}
