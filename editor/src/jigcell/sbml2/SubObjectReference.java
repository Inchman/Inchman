package jigcell.sbml2;

/**
 * A reference to a subObject contained within an object in the model.  This type of reference is used to specify the subObject in a link.
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Ranjit Randhawa
 */

public final class SubObjectReference extends SimpleObjectReference {

   public SubObjectReference () {
      this ((SimpleObjectReference) null);
   }

   /**
    * Creates a copy of this existing object.  This does not copy metadata, such as annotations or notes.
    *
    * @param reference Reference, must not be null
    */

   public SubObjectReference (SubObjectReference reference) {
      this ();
      setObject (reference.getObject ());
   }
   
   public SubObjectReference (SimpleObjectReference object) {
      super ();
      setObject (object);
   }
   
   public SubObjectReference (String object) {
      super ();
      setObject (object);
   }
   
   protected XMLPrinter print (XMLPrinter parent) {
      return print (parent, "subobject");
   }
}
