package jigcell.sbml2;

/**
 * A port reference to a object contained in the model.  This type of reference is used to speficy the target object in a port.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Ranjit Randhawa
 */

public final class TargetObjectReference extends SimpleObjectReference {

   private SubObjectReference subObjectReference;
   private final SBase subObjectElement;

   public TargetObjectReference () {
      this ((SimpleObjectReference) null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param reference Reference, must not be null
    */

   public TargetObjectReference (TargetObjectReference reference) {
      this ();
      setObject (reference.getObject ());
   }

   public TargetObjectReference (String object) {
      super ();
      setObject (object);
      subObjectReference = new SubObjectReference ();
      subObjectElement = new SBase ();
   }
   
   public TargetObjectReference (SimpleObjectReference object) {
      super ();
      setObject (object);
      subObjectReference = new SubObjectReference ();
      subObjectElement = new SBase ();
   }
   
   public SubObjectReference getSubObjectReference () {
      return subObjectReference;
   }

   public SBase getSubObjectElement () {
      return subObjectElement;
   }
   
   public void setSubObjectReference (String reference) {
      if (reference != null)
         subObjectReference.setObject (reference);
   }

   public void setSubObjectReference (SubObjectReference reference) {
      if (reference != null)
         subObjectReference = reference;
   }

   protected XMLPrinter print (XMLPrinter parent) {
      XMLPrinter printer = super.print (parent, "target");      
      if (subObjectReference.getObject () != null)
         printer.addText (getSubObjectReference ().toString ());      
      return printer;
   }
}
