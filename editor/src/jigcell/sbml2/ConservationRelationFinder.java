package jigcell.sbml2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Detects, enumerates, and validates the conservation relations in a model.  The conservation relations detected by this class are not
 * updated in response to modifications to the model.  After changing the model, a new ConservationRelationFinder must be created.
 * <p>
 * augmentedStoichiometryMatrix has the following structure (where the first index is for
 * the vertical dimension and the second index is for the horizontal dimension):
 *
 * <pre>
 *     `,          reactions      |   species     |   species     |
 *       `,         (size=        |    (size=     |    index      |
 *         `,        num          |     num       |    (size=1)   |
 *           `,      reactions)   |     species)  |               |
 *             `r-----------------+---------------+---------------+-
 *              |                 |               |               |
 *     species  |                 |               |               |
 *              |                 |               |               |
 *              |                 |               |               |
 *              |                 |               |               |
 *    __________L_________________L_______________L_______________L_
 * </pre>
 * </p>
 *
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See LICENSE for more details.
 * </p>
 *
 * @author Nicholas Allen
 */

public final class ConservationRelationFinder {
   /**
    * Values of true at index i means species i cannot be used in a conservation
    * relation.  This is the case if the species is set by a rule or event.  Any
    * generating vector that has an illegal species cannot be used as
    * conservation relations.  Any user entered conservation relations also may
    * not contain illegal species.
    */
   private final boolean illegalSpecies [];
   private final double augmentedStoichiometryMatrix [][];
   private final double stoichiometryMatrix [][];
   private final double generatingVectors [][];
   private final int conservationRelationCount;
   private final List speciesIds;
   private boolean stoichiometryMatrixFilled = false;

   /**
    * Represents a single conservation relation.
    */

   public final static class ConservationRelation {
      private final String dependentSpeciesId;
      private final Map conservationRelation;

      /**
       * Creates a new conservation relation.
       *
       * @param conservationRelation Map[String,Double] from species ids to the quantity of the species involved in the relation
       * @param dependentSpeciesId Id of the species being eliminated by this
       * conservation relation, or null if this conservation relation is simply
       * recorded for possible future use.
       */

      public ConservationRelation (Map conservationRelation, String dependentSpeciesId) {
         if (conservationRelation == null)
            throw new IllegalArgumentException ();
         if (dependentSpeciesId != null && !conservationRelation.containsKey (dependentSpeciesId))
            throw new IllegalArgumentException ("Dependent species must be included in the conservation relation.");
         for (Iterator iterator = conservationRelation.keySet ().iterator (); iterator.hasNext (); ) {
            String speciesId = (String) iterator.next ();
            Double value = (Double) conservationRelation.get (speciesId);
            if (value.isNaN () || value.isInfinite () || value.doubleValue () == 0.0)
               throw new IllegalArgumentException ("Coefficient for species " + speciesId + " is illegal.");
         }
         this.conservationRelation = Collections.unmodifiableMap (new TreeMap (conservationRelation));
         this.dependentSpeciesId = dependentSpeciesId;
      }

      public Map getConservationRelation () {
         return conservationRelation;
      }

      public String getDependentSpeciesId () {
         return dependentSpeciesId;
      }

      public String toString () {
         return "(" + dependentSpeciesId + ")" + conservationRelation.toString ();
      }
   }

   public static double calcGCF( double a, double b ) {
      int smaller = Math.abs((int)Math.rint(a));
      int larger  = Math.abs((int)Math.rint(b));
      if ( smaller == 0 ) return larger;
      if ( larger == 0  ) return smaller;
      int facs = 1; // factor small
      int facl;     // factor large
      int gcf = 1;
      if ( smaller > larger ) { int temp = smaller; smaller = larger; larger = temp; }
      facl = smaller / facs;
      while ( !commonFactor(smaller,larger,facl) && facl > facs ) {
         facs++;
         facl = smaller / facs;
         if ( commonFactor(smaller,larger,facs) ) gcf = facs;
      }
      if ( facl > facs ) return (double)facl;
      return (double)gcf;
   }

   public static boolean commonFactor(int a,int b,int factor) {
      return a/factor*factor == a && b/factor*factor == b;
   }

   private static int countDependentRows (double matrix [][], int rowCount, int columnCount) {
      for (int rowIndex = rowCount - 1; rowIndex >= 0; rowIndex--)
         for (int reactionIndex = 0; reactionIndex < columnCount; reactionIndex++)
            if (matrix [rowIndex][reactionIndex] != 0.0)
               return rowCount - rowIndex - 1;
      return rowCount;
   }

   private static void multiplyRow (double matrix [][], int entryCount, int row1, double coefficient) {
      double row [] = matrix [row1];
      for (int entryIndex = 0; entryIndex < entryCount; entryIndex++)
         row [entryIndex] *= coefficient;
   }

   private static void multiplyAndAddToRow (double matrix [][], int entryCount, int row1, int row2, double coefficient) {
      double fromRow [] = matrix [row1];
      double toRow [] = matrix [row2];
      for (int entryIndex = 0; entryIndex < entryCount; entryIndex++)
         toRow [entryIndex] += coefficient * fromRow [entryIndex];
   }

   /**
    * Takes a stoichiometry matrix and calculates all the generating vectors for
    * conservation relations.  The result is returned as double[][] where each
    * double[i] is a generating vector and each double[i][j] is the jth element
    * in the ith generating vector.  Note, there may be more generating vectors
    * than there are dependent rows in the stoichiometry matrix.  This would
    * mean that not all generating vectors can be used as a conservation
    * relations.  For example, the system
    * <pre>
    *    a+b -&gt; c+d
    *    c+d+e -&gt; f
    * </pre>
    * Would have generating vectors of:
    * <pre>
    *    a b c d e f
    *    0 0 0 0 1 1
    *    1 0 1 0 0 1
    *    1 0 0 1 0 1
    *    0 1 1 0 0 1
    *    0 1 0 1 0 1
    * </pre>
    * But the stoichiometry matrix has 6 rows and rank 2.  This means only 4
    * conservation relations can be used, yet 5 can be generated by the
    * generating vectors.
    * <p>
    * The algorithm employed here is described in "Determining All Extreme
    * Semi-positive Conservation Relations in Chemical Reaction Systems: A Test
    * Criteria for Conservativity" by Stefan Schuster and Thomas Höfer, J. Chem.
    * Soc. Faraday Trans., 1991 87(16) pp 2561-2566.
    */
   private static double[][] findAllGeneratingVectors( double stoichMat[][] ) {
      int reactionCount = stoichMat[0].length;
      int speciesCount  = stoichMat.length;
      int columnCount   = reactionCount+speciesCount;
      double mat[][] = new double[stoichMat.length][];
      // Copy the matrixArg in to the work array mat and augment with identity
      // matrix.
      for ( int i = 0; i < stoichMat.length; i++ ) {
         mat[i] = new double[columnCount];
         System.arraycopy(stoichMat[i],0,mat[i],0,reactionCount);
         // makes identity matrix in augmented section, Java initializes double
         // arrays to 0.0
         mat[i][reactionCount+i] = 1.0;
      }
      // Reduce system and find generating vectors.
      List<Set> zeroSets = new ArrayList<Set>(speciesCount);
      for ( int j = 0; j < reactionCount; j++ ) {
         // calculate all zero sets:
         // S(i) = { h | mat[i][h] == 0.0, for h = reactionCount ... columnCount }
         zeroSets.clear();
         for ( int i = 0; i < mat.length; i++ ) {
            Set<Integer> zeroSet = new HashSet<Integer>(speciesCount*2);
            for ( int h = reactionCount; h < columnCount; h++ ) {
               if ( mat[i][h] == 0.0 ) zeroSet.add(new Integer(h));
            }
            zeroSets.add(zeroSet);
         }
         double[][] newRows = new double[0][]; // there may be no generating vectors
         for ( int i = 0; i < mat.length; i++ ) {
            for ( int k = i+1; k < mat.length; k++ ) {
               if ( mat[i][j] * mat[k][j] >= 0.0 ) continue;
               // If for each row l the row k and i have a zero in a column that
               // l doesn't (the column could be different for different l's),
               // then use this row as a new vector.
               boolean hasUniqueZeroColumn = true;
               for ( int l = 0; l < mat.length; l++ ) {
                  if ( l==i || l==k ) continue;
                  Set<Integer> working = new HashSet<Integer>(zeroSets.get(i));
                  working.retainAll(zeroSets.get(k)); // intersection
                  working.removeAll(zeroSets.get(l)); // difference
                  if ( working.size() == 0 ) hasUniqueZeroColumn = false;
               }
               if ( !hasUniqueZeroColumn ) continue;
               double[][] temp = newRows;
               newRows = new double[temp.length+1][];
               System.arraycopy(temp,0,newRows,0,temp.length);
               double [] newRow = new double[columnCount];
               double gcf = 0.0;
               for ( int ik = 0; ik < columnCount; ik++ ) {
                  newRow[ik] = Math.abs(mat[i][j])*mat[k][ik] +
                               Math.abs(mat[k][j])*mat[i][ik];
                  gcf = calcGCF( gcf, newRow[ik] );
               }
               for ( int l = 0; l < columnCount; l++ ) {
                  newRow[l] = Math.rint(newRow[l]/gcf);
               }
               newRows[newRows.length-1] = newRow;
            }
         }
         int zeroRowCount = 0;
         for ( int l = 0; l < mat.length; l++ ) {
            if ( mat[l][j] == 0.0 ) zeroRowCount++;
         }
         double[][] temp = mat;
         mat = new double[zeroRowCount+newRows.length][];
         int m = 0;
         for ( int l = 0; l < temp.length; l++ ) {
            if ( temp[l][j] == 0.0 ) {
               mat[m] = temp[l];
               m++;
            }
         }
         System.arraycopy(newRows,0,mat,zeroRowCount,newRows.length);
      }
      // Copy the generating vectors out of the working matrix mat.
      double[][] genVectors = new double[mat.length][];
      for ( int i = 0; i < mat.length; i++ ) {
         genVectors[i] = new double[speciesCount];
         System.arraycopy(mat[i],reactionCount,genVectors[i],0,speciesCount);
      }
      return genVectors;
   }

   /**
    * Reduces an augmented matrix to reduced row echelon form.  The matrix up to
    * columnCount is taken as the matrix to reduce to row echelon form.  The
    * augmented part of the matrix (columnCount to entryCount) will have all the
    * same row operations performed on it, but the augmented part will not
    * influence the algorithm.  This
    * method can be used for a few purposes, one is to reduce the stoichiometry
    * matrix augmented with the indentity matrix.  In this case dependent rows
    * will be all zero at the bottom of the reduced stoichiometry matrix and
    * conservation relations will appear at the same bottom part in the
    * augmented (originally the identity matrix) portion of the matrix.  Note,
    * using this method to generate conservation relations may not find all
    * conserved moieties and may not find physically meaningful conserved
    * quantities (i.e., negative values may appear as conserved quantities).
    * @see findAllGeneratingVectors
    */
   private static void reduceMatrix (double matrix [][], int rowCount, int columnCount, int entryCount) {
      int linkCount = Math.min (rowCount, columnCount);
      int activeRow = 0;
      for (int columnIndex = 0; columnIndex < columnCount; columnIndex++) {
         if (activeRow >= linkCount)
            break;
         double coefficient = 0.0;
         int targetRow = -1;
         for (int rowIndex = activeRow; rowIndex < rowCount; rowIndex++) {
            coefficient = matrix [rowIndex][columnIndex];
            if (coefficient != 0.0) {
               targetRow = rowIndex;
               break;
            }
         }
         if (targetRow == -1)
            continue;
         if (activeRow != targetRow)
            swapRows (matrix, activeRow, targetRow);
         if (coefficient != 1.0)
            multiplyRow (matrix, entryCount, activeRow, 1.0 / coefficient);
         for (int rowIndex = activeRow + 1; rowIndex < rowCount; rowIndex++) {
            coefficient = matrix [rowIndex][columnIndex];
            if (coefficient != 0.0)
               multiplyAndAddToRow (matrix, entryCount, activeRow, rowIndex, -coefficient);
         }
         activeRow++;
      }
   }

   private static void swapRows (double matrix [][], int row1, int row2) {
      double row [] = matrix [row1];
      matrix [row1] = matrix [row2];
      matrix [row2] = row;
   }

   /**
    * Computes the conservation relations for a model.
    *
    * @param model SBML model
    */

   public ConservationRelationFinder (Model model) {
      if (model == null)
         throw new IllegalArgumentException ("Model is null.");
      List speciesList = model.getSpecies ();
      if (speciesList == null)
         throw new IllegalArgumentException ("Model does not have a list of species.");
      List reactions = model.getReactions ();
      if (reactions == null)
         throw new IllegalArgumentException ("Model does not have a list of reactions.");
      int speciesCount = speciesList.size ();
      int reactionCount = reactions.size ();
      int entryCount = speciesCount + reactionCount;
      illegalSpecies               = new boolean [speciesCount];
      augmentedStoichiometryMatrix = new double [speciesCount][];
      stoichiometryMatrix          = new double [speciesCount][];
      speciesIds                   = new ArrayList (speciesCount);
      for ( int i = 0; i < speciesCount; i++ ) { stoichiometryMatrix[i] = new double[reactionCount]; }
      fillStoichiometryMatrix (model, speciesList, reactions, reactionCount);
      fillAugmentedStoichiometryMatrix (model, speciesList, reactions, speciesCount, reactionCount, entryCount);
      reduceMatrix (augmentedStoichiometryMatrix, speciesCount, reactionCount, entryCount);
      conservationRelationCount = countDependentRows (augmentedStoichiometryMatrix, speciesCount, reactionCount);
      generatingVectors = findAllGeneratingVectors(stoichiometryMatrix);
   }

   /**
    * Computes the List[ConservationRelation] of conservation relations in the
    * model from the generating vectors.
    */

   /*
    * This is the old method.  This one uses the reduceMatrix output to create
    * conservation relations.  This method has been replaced because it often
    * gives physically meaningless (albeit correct) conservation relations;
    * these physically meaningless relations have negative values.
   public List getConservationRelations () {
      List equations = new ArrayList ();
      int speciesCount = augmentedStoichiometryMatrix.length;
      if (speciesCount == 0)
         return equations;
      int speciesColIndex = augmentedStoichiometryMatrix [0].length - 1;
      int reactionCount = speciesColIndex - speciesCount;
outer:
      for (int rowIndex = speciesCount - 1, lastRow = rowIndex - conservationRelationCount; rowIndex > lastRow; rowIndex--) {
         double row [] = augmentedStoichiometryMatrix [rowIndex];
         Map conservationRelation = new TreeMap ();
         for (int inverseIndex = reactionCount; inverseIndex < reactionCount + speciesCount; inverseIndex++) {
            double coefficient = row [inverseIndex];
            if (coefficient == 0.0)
               continue;
            int speciesIndex = inverseIndex - reactionCount;
            if (illegalSpecies [speciesIndex])
               continue outer;
            conservationRelation.put (speciesIds.get (speciesIndex), new Double (coefficient));
         }
         if (conservationRelation.size () > 1)
            equations.add (new ConservationRelation (conservationRelation,
               (String) speciesIds.get ((int) augmentedStoichiometryMatrix [rowIndex][speciesColIndex])));
      }
      assert validateConservationRelations (equations);
      return equations;
   }
   */

   /**
    * Uses the generatingVectors to create conservation relations.  Not all the
    * generatingVectors will necessarilly be used.  Only up to
    * conservationRelationCount will be used.
    */
   public List getConservationRelations () {
      Set<String> usedAsDependentSpecies = new HashSet<String>();
      List equations = new ArrayList ();
      int speciesCount = stoichiometryMatrix.length;
      if (speciesCount == 0)
         return equations;
      int reactionCount = stoichiometryMatrix[0].length;
      int generatingVectorCount = generatingVectors.length;
      if ( generatingVectorCount < conservationRelationCount )
         throw new RuntimeException("There are less generating vectors than conservation relations.  There must be a bug in the code.");
outer:
      for ( int rowIndex = 0; rowIndex < generatingVectorCount; rowIndex++ ) {
         double gv [] = generatingVectors [rowIndex];
         Map conservationRelation = new TreeMap ();
         for (int colIndex = 0; colIndex < speciesCount; colIndex++) {
            double coefficient = gv [colIndex];
            if (coefficient == 0.0)
               continue;
            if (illegalSpecies [colIndex])
               continue outer;
            conservationRelation.put (speciesIds.get (colIndex), new Double (coefficient));
         }
         if (conservationRelation.size () > 1) {
            System.out.println("Got a conservation relation.");
            equations.add (new ConservationRelation (conservationRelation,
               getNextAvailableDependentSpecies(usedAsDependentSpecies,conservationRelation)));
         }
      }
      assert validateConservationRelations (equations);
      return equations;
   }

   private static String getNextAvailableDependentSpecies(
      Set usedAsDependentSpecies,
      Map conservationRelation
   ) {
      for ( Iterator i = conservationRelation.keySet().iterator(); i.hasNext(); ) {
         String id = (String)i.next();
         if ( !usedAsDependentSpecies.contains(id) ) {
            usedAsDependentSpecies.add(id);
            System.out.println("   Using dependent species: "+id);
            return id;
         }
      }
      return null;
//      throw new RuntimeException("Could not find an available species to use as the dependent species.");
   }

   /**
    * Checks that the given conservation relations are consistent with the
    * conservation relations found in the model.  All ids must be ids of
    * species.  These relations don't include the parameter representing the
    * conserved quantity.
    *
    * @param userRelations List[ConservationRelation] of proposed conservation relations
    */

   public boolean validateConservationRelations (List userRelations) {
      if (userRelations == null)
         throw new IllegalArgumentException ("List of conservation relations is null.");
      int userRelationCount = userRelations.size ();
      if (userRelationCount > conservationRelationCount)
         throw new IllegalArgumentException ("Too many conservation relations.");
      if (userRelationCount == 0)
         return true;
      int totalRelationCount = userRelationCount + conservationRelationCount;
      double userSystem [][] = new double [userRelationCount + conservationRelationCount][];
      int speciesCount = augmentedStoichiometryMatrix.length;
      Set usedSpeciesIds = new HashSet ();
      for (int userRelationIndex = 0; userRelationIndex < userRelationCount; userRelationIndex++) {
         ConservationRelation userRelation = (ConservationRelation) userRelations.get (userRelationIndex);
         String dependentSpeciesId = userRelation.dependentSpeciesId;
         if (!speciesIds.contains (dependentSpeciesId))
            throw new IllegalArgumentException ("Species with id '" + dependentSpeciesId + "' is not present in the model.");
         if (!usedSpeciesIds.add (dependentSpeciesId))
            throw new IllegalArgumentException ("Species with id '" + dependentSpeciesId + "' is used as a dependent variable twice.");
         double relationVector [] = new double [speciesCount];
         for (Iterator entryIterator = userRelation.getConservationRelation ().entrySet ().iterator (); entryIterator.hasNext (); ) {
            Map.Entry entry = (Map.Entry) entryIterator.next ();
            String speciesId = (String) entry.getKey ();
            int speciesIndex = speciesIds.indexOf (speciesId);
            if (speciesIndex == -1)
               throw new IllegalArgumentException ("Species with id '" + speciesId + "' is not present in the model.");
            if (illegalSpecies [speciesIndex])
               throw new IllegalArgumentException ("Species with id '" + speciesId + "' cannot be used in a conservation relation.");
            if (relationVector [speciesIndex] != 0.0)
               throw new IllegalArgumentException ("Species with id '" + speciesId + "' is used in the conservation relation twice.");
            relationVector [speciesIndex] = ((Double) entry.getValue ()).doubleValue ();
         }
         userSystem [userRelationIndex] = relationVector;
      }
      reduceMatrix (userSystem, userRelationCount, speciesCount, speciesCount);
      if (countDependentRows (userSystem, userRelationCount, speciesCount) > 0)
      {
         System.err.println("ERROR: User conservation relations are not linearly independent.");
         return false;
      }
      int reactionCount = augmentedStoichiometryMatrix [0].length - speciesCount - 1;
      for (int conservationRelationIndex = 0; conservationRelationIndex < conservationRelationCount; conservationRelationIndex++) {
         double relationVector [] = new double [speciesCount];
         System.arraycopy (augmentedStoichiometryMatrix [speciesCount - 1 - conservationRelationIndex], reactionCount, relationVector, 0, speciesCount);
         userSystem [userRelationCount + conservationRelationIndex] = relationVector;
      }
      reduceMatrix (userSystem, totalRelationCount, speciesCount, speciesCount);
      return countDependentRows (userSystem, totalRelationCount, speciesCount) == userRelationCount;
   }

   private void addSpeciesStoichiometry (Model model, List modelSpecies, int reactionIndex, SpeciesReference reference, boolean positive) {
      if (reference == null)
         throw new IllegalArgumentException (" contains an incomplete species reference.");
      Species species = reference.getSpecies (model);
      if (species == null)
         throw new IllegalArgumentException (" contains an incomplete species reference.");
      if (species.isConstant () || species.isBoundaryCondition ())
         return;
      int speciesIndex = modelSpecies.indexOf (species);
      if (speciesIndex == -1)
         throw new IllegalArgumentException (" contains species " + species.getName () + " not defined in the model.");
      double stoichiometry = reference.getStoichiometry ();
      if (Double.isNaN (stoichiometry))
         throw new IllegalArgumentException (" contains species " + species.getName () + " without fixed stoichiometry.");
      // multiply stoichiometry by compartment volume.  In case this reaction
      // involves multiple compartments, this will make the conservation
      // relations correct for multiple compartments.
      double compartmentSize = species.getCompartment(model).getSize();
      if ( !Double.isNaN(compartmentSize) && compartmentSize != 1.0 ) {
         stoichiometry /= compartmentSize;
      }
      stoichiometryMatrix [speciesIndex][reactionIndex] += positive ? stoichiometry : -stoichiometry;
   }

   private void fillStoichiometryMatrix( Model model, List speciesList, List reactions, int reactionCount ) {
      for (int reactionIndex = 0; reactionIndex < reactionCount; reactionIndex++) {
         Reaction reaction = (Reaction) reactions.get (reactionIndex);
         if (reaction == null)
            continue;
         try {
            for (Iterator productIterator = reaction.getProduct ().iterator (); productIterator.hasNext (); )
               addSpeciesStoichiometry (model, speciesList, reactionIndex, (SpeciesReference) productIterator.next (), true);
            for (Iterator reactantIterator = reaction.getReactant ().iterator (); reactantIterator.hasNext (); )
               addSpeciesStoichiometry (model, speciesList, reactionIndex, (SpeciesReference) reactantIterator.next (), false);
         } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException ("Reaction " + reaction.getName () + e.getMessage ());
         }
      }
      stoichiometryMatrixFilled = true;
   }

   private void fillAugmentedStoichiometryMatrix (Model model, List speciesList, List reactions, int speciesCount, int reactionCount, int entryCount) {
      if ( !stoichiometryMatrixFilled ) throw new RuntimeException("Must call fillStoichiometryMatrix before calling fillAugmentedStoichiometryMatrix.");
      for (int speciesIndex = 0; speciesIndex < speciesCount; speciesIndex++) {
         double speciesRow [] = new double [entryCount + 1];
         speciesRow [speciesIndex + reactionCount] = 1.0;
         speciesRow [entryCount] = speciesIndex;
         augmentedStoichiometryMatrix [speciesIndex] = speciesRow;
         Species species = (Species) speciesList.get (speciesIndex);
         speciesIds.add (species.getId ());
         illegalSpecies [speciesIndex] = species.isSetByRule (model) || species.isSetByEvent (model);
      }
      for ( int i = 0; i < speciesCount; i++ ) {
         System.arraycopy(stoichiometryMatrix[i],0,augmentedStoichiometryMatrix[i],0,reactionCount);
      }
   }
}
