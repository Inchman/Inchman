package jigcell.sbml2.math;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;

/**
 * SymbolTable keeps track of all names and strings in a model and their
 * translated names in Fortran.  The names are translated by making the strings
 * Fortran safe.  See makeFortranIdSafe.  The symbol table can make strings safe
 * as Fortran identifiers or safe as Fortran character strings.
 * <p>
 * This code is licensed under the DARPA BioCOMP Open Source License.  See
 * LICENSE for more details.
 *
 * @author Jason Zwolak
 * @author Nicholas Allen
 */

public class SymbolTable {

    /*
     * The following constants deal with mutating strings when there already
     * exists a string of the same name.  See the routine findAFortranName for
     * more documentation.
     */
    public static final int ID_START = 1;
    public static final int ID_END = 9999;
    public static final int ID_WIDTH = 4;

    /**
     * These constants describe the type of symbol table.  A FORTRAN_IDENTIFIER
     * symbol table will create symbols suitable as Fortran identifiers.  A
     * FORTRAN_CHARACTER_STRING symbol table will create symbols suitable as
     * Fortran character strings.  IDENTITY will not modify the name at all.  It
     * will simply copy it to "fortranName".  It will not modify the name even
     * if it is identical to another fortranName.
     */
    public static final class Mode {
        public static final Mode FORTRAN_IDENTIFIER_UNDERSCORED =
            new Mode( "FORTRAN_IDENTIFIER_UNDERSCORED", 31, true, false)
        ;
        public static final Mode FORTRAN_IDENTIFIER       =
            new Mode( "FORTRAN_IDENTIFIER", 31, false, false)
        ;
        /**
         * The mode for Fortran character strings.  This mode limits characters
         * to the length of 65000.  This is an empirically determined number for
         * the the Lahey-Fujitsu Fortran 95 compiler v6.1 Express for Linux.
         */
        public static final Mode FORTRAN_CHARACTER_STRING =
            new Mode( "FORTRAN_CHARACTER_STRING", 65000, false, false )
        ;
        public static final Mode IDENTITY =
            new Mode( "IDENTITY", Integer.MAX_VALUE, false, false )
            //new Mode( "IDENTITY", (1<<31) - 1 )
        ;
        /**
         * Allows multiple symbols to map to the different ids and names.  If
         * the id 'k' maps to 'k_1' (k-&gt;k_1) then another id, say 'a', can
         * also map to 'k_1' (a-&gt;k_1).  This does not allow duplicate names
         * and ids.
         */
        public static final Mode IDENTITY_WITH_DUPLICATES =
            new Mode( "IDENTITY", Integer.MAX_VALUE, false, true )
            //new Mode( "IDENTITY", (1<<31) - 1 )
        ;
        private String name;
        private int maxLength;
     /*
     * appendUnderscore specifies whether to append an '_' to the end of symbols
     * to ensure the symbol does not conflict with the target environments
     * predefined symbols (e.g., Fortran or XPP).
     *
     */
        private boolean appendUnderscore;
        private boolean allowDuplicates;
        public Mode(
            String name,
            int maxLength,
            boolean appendUnderscore,
            boolean allowDuplicates
        ) {
            this.name             = name;
            this.maxLength        = maxLength;
            this.appendUnderscore = appendUnderscore;
            this.allowDuplicates  = allowDuplicates;
        }
        public int maxLength() { return maxLength; }
        public boolean appendUnderscore() { return appendUnderscore; }
        public boolean allowDuplicates() { return allowDuplicates; }
        public String toString() { return name; }
    }
    /*
     * forbiddenSymbols ares symbols already defined in the Fortran Code and
     * cannot be used by the modeller.
     *
     * Note this isn't used yet, it is just a place holder and reminder.
     */
    private static Hashtable forbiddenSymbols;

    /*
     * nameCol - list of names in the symbol table.
     * fortranNameCol - list of fortran names in the symbol table.
     * idCol - list of ids in the symbol table
     * nameColRev - list of names in the symbol table.
     * fortranNameColRev - list of fortran names in the symbol table.
     *
     * nameCol, fortranNameCol, and idCol form the symbol table with nameCol
     * containing the names that are keys.
     *
     * nameColRev and fortranNameColRev provide reverse lookup of names based on
     * fortran names.
     *
     * nameCol and fortranNameColRev are alphabetically sorted.  fortranNameCol,
     * idCol, and nameColRev are ordered so they match up with their
     * counterpart names in nameCol and fortranNameColRev.
     *
     * fortranNameColRev has all its entries capitolized.  Fortran is case
     * insensitive.  Therefore, the names must be unique regardless of case.
     * Capitolizing all entries in fortranNameColRev ensures this.
     *
     * fortranNameCol contains the names that will be written to the Fortran
     * file.
     *
     * idNameMap contains a hash of id-name pairs.  This is for quickly looking
     * up names given an id.  Note that the ids are assumed to be unique (an
     * exception is thrown otherwise, see code below to ensure an exception is
     * really thrown).
     */
    private ArrayList nameCol;
    private ArrayList fortranNameCol;
    private ArrayList idCol;
    private ArrayList nameColRev;
    private ArrayList fortranNameColRev;
    private Hashtable idNameMap;

    /*
     * maxLength specifies the maximum length of symbols in this table.  After
     * a symbol is mangled it is guaranteed to be this length or smaller.
     *
     * mode specifies the mode of this symbol table.  The mod can be
     * TYPE_IDENTIFIER or TYPE_CHARACTER_STRING.  See TYPE_IDENTIFIER and
     * TYPE_CHARACTER_STRING for more info.
     *
     * If padToLength is true then FORTRAN_CHARACTER_STRINGs will be padded with
     * the padCharacter to the maxLength.
     */
    private int maxLength;
    private Mode mode = Mode.FORTRAN_IDENTIFIER;
    private boolean padToLength = false;
    private char padCharacter = ' ';

    /**
     * Create a new empty symbol table.
     */
    public SymbolTable() {
        maxLength         = mode.maxLength();
        nameCol           = new ArrayList();
        fortranNameCol    = new ArrayList();
        idCol             = new ArrayList();
        nameColRev        = new ArrayList();
        fortranNameColRev = new ArrayList();
        idNameMap         = new Hashtable();
    }

    /**
     * Create a new SymbolTable object with specified maximum ID length.  The
     * user may wish to limit the maximum ID length to smaller than Fortran
     * 95's limit on the length.  This is useful if the user wishes to prepend
     * or append strings with the '__' (double underscore) separator for a new
     * unique ID based on an old ID.  This is also useful if Fortran is not the
     * target for code generation.  The maxLength cannot be greater than the
     * default for the given mode.
     */
    public SymbolTable( int maxLength ) {
        this();
        this.maxLength =
            maxLength > mode.maxLength() ?
            mode.maxLength() : maxLength
        ;
    }

    /**
     * Creates a new empty SymbolTable object with the specified mode.  The mode
     * will control the nature in which the symbol table mangles symbols.  See
     * the {@link Mode} documentation for more information.
     */
    public SymbolTable( Mode mode ) {
        this();
        this.mode = mode;
        this.maxLength = mode.maxLength();
    }

    /**
     * Creates a new empty SymbolTable with the specified mode and maximum id
     * length.
     */
    public SymbolTable( Mode mode, int maxLength ) {
        this();
        this.mode = mode;
        this.maxLength = maxLength;
    }

    /**
     * Creates a new empty SymbolTable with the specified mode and maximum id
     * length and request to pad.
     */
    public SymbolTable( Mode mode, int maxLength, boolean padToLength ) {
        this();
        this.mode = mode;
        this.maxLength = maxLength;
        this.padToLength = padToLength;
    }

    /**
     * Write the symbol table to out.
     */
    public void dumpSymbolTable( PrintWriter out ) {
        Iterator n = nameCol.iterator();
        Iterator f = fortranNameCol.iterator();
        while ( n.hasNext() && f.hasNext() ) {
            String name = (String)n.next();
            out.println(name+" ("+lookupId(name)+")"+" => "+f.next());
        }
        out.flush();
    }
    public void dumpSymbolTable( ) {
        dumpSymbolTable(new PrintWriter(System.out));
    }

    public void emptySymbolTable(){
        nameCol.clear();
        fortranNameCol.clear();
        idCol.clear();
        nameColRev.clear();
        fortranNameColRev.clear();
        idNameMap.clear();
    }

    public Map getSymbolToIdMap() {
        Map map = new HashMap();
        for (int i=0; i<fortranNameCol.size(); i++) {
            map.put(fortranNameCol.get(i),idCol.get(i));
        }
        return map;
    }

    /**
     * Find the fortranName for name.  If name doesn't exist in the current
     * symbol table then return null otherwise return the fortran name.
     *
     * @return Fortran name on success, null on error.
     */
    public String lookupSymbol( String name ) {
        if ( name == null ) return null;
        int index = lookupSymbolIndex( name, nameCol );
        if ( index < 0 ) return null;
        return (String)fortranNameCol.get(index);
    }
    /**
     * Returns the Fortran name for id.
     */
    public String lookupSymbolById( String id ) {
        String name = (String)idNameMap.get(id);
        if ( name == null )
            return null;
        String ret = lookupSymbol(name);
        return ret;
    }
    /**
     * Returns the id for the JigCell name.
     */
    public String lookupId( String name ) {
        if ( name == null ) return null;
        int index = lookupSymbolIndex( name, nameCol );
        if ( index < 0 ) return null;
        return (String)idCol.get(index);
    }
    /**
     * Returns the JigCell name for id.
     */
    public String lookupName( String id ) {
        return (String)idNameMap.get(id);
    }
    /**
     * Returns the JigCell id for fortranName.  Throws an exception if
     * <i>allowDuplicates</i> is true from the <i>Mode</i>.
     */
    public String reverseLookupId( String fortranName ) throws Exception {
        if ( mode.allowDuplicates() ) {
            throw new Exception("reverseLookupId does not work when "+
                "duplicates are allowed.");
        }
        int index = lookupFortranNameIndex( fortranName );
        if ( index < 0 ) return null;
        return lookupId( (String)nameColRev.get(index) );
    }


    /**
     * Finds fname in fortranNameColRev.  Converts fname to upper case first
     * because all entries in fortranNameColRev are upper case.  Calls
     * lookupSymbolIndex to perform the search.
     */
    private int lookupFortranNameIndex( String fname ) {
        if ( mode == Mode.IDENTITY )
            return lookupSymbolIndex( fname, fortranNameColRev );
        return lookupSymbolIndex( fname.toUpperCase(), fortranNameColRev );
    }

    /**
     * Finds the index of name in the internal symbol table.  If name doesn't
     * exist in the table then the index returned is the negative of the place
     * in the table where this symbol should be placed -1; the minus one makes
     * the code work in the case where the symbol should be inserted at the 0th
     * position.
     *
     * @return Index of name if name is in the symbol table, negative position
     * in the table otherwise.
     */
    private int lookupSymbolIndex( String name, ArrayList array ) {
        int left = 0, right = array.size();
        int next = ( left + right ) / 2;
        if ( right == left ) return -1;
        int compare = 0;
        if ( next < array.size() )
            compare = name.compareTo( (String)array.get(next) );
        while ( right > left && compare != 0 ) {
            if ( compare > 0 ) left  = next + 1;
            else               right = next;
            next = ( left + right ) / 2;
            if ( next < array.size() )
                compare = name.compareTo( (String)array.get(next) );
        }
        if ( compare == 0 ) return next;
        return -next-1;
    }

    /**
     * Add <i>id</i> to the symbol table.  If name already exists then throw an
     * exception.  This routine adds a symbol with the id <i>id</i> and the
     * name <i>id</i>.  In other words, the id is copied to the name since the
     * name isn't present and symbols are added by name.
     *
     * @param id The ID of the symbol table entry.
     * @return The Fortran name of name.
     */
    public String addSymbol( String id ) throws Exception {
        return addSymbol( id, id );
    }
    /**
     * Add a symbol with name and id.  All ids entered here must be unique or
     * an exception is thrown.
     */
    public String addSymbol( String name, String id ) throws Exception {
        String fortranName = makeFortranSafe( name, mode, maxLength );
        int fortranIndex = lookupFortranNameIndex( fortranName );
        if ( fortranIndex >= 0 && mode != Mode.IDENTITY )
            fortranName = findAFortranName( fortranName );
        return addSymbol(name,id,fortranName);
    }
    /**
     * Add a symbol with name, id, and specific symbol.
     * This allows caller to specify symbol, but the symbol must be unique!
     */
    public String addSymbol( String name, String id, String symbol )
        throws Exception
    {
        if ( name == null ) {
            System.err.println(
                "WARNING: symbol name is null, setting name = id"
            );
            name = id;
        }
        if ( id   == null ) throw new Exception("id is null");
        if ( nameCol.size() == Integer.MAX_VALUE )
            throw new Exception( "Symbol table is full." );
        int index        = lookupSymbolIndex( name, nameCol );
        if ( index >= 0 )
            throw new Exception(
                "'"+name+"' already exists in the symbol table as a name." );
        int forIndex = lookupFortranNameIndex( symbol );
        if ( forIndex >= 0 && !mode.allowDuplicates() ) {
            throw new Exception(
                "'"+symbol+"' already exists in the symbol table as a "+
                "symbol." );
        } else if ( forIndex >= 0 ) {
            forIndex = -forIndex-1; // just so the code below looks clean.  The
                                    // code below is expecting the negative of
                                    // the position the symbol will be inserted
                                    // into.
        }
        nameCol.add(           -index-1,    name                      );
        fortranNameCol.add(    -index-1,    symbol                    );
        idCol.add(             -index-1,    id                        );
        nameColRev.add(        -forIndex-1, name                      );
        if ( mode == Mode.IDENTITY )
            fortranNameColRev.add( -forIndex-1, symbol );
        else
            fortranNameColRev.add( -forIndex-1, symbol.toUpperCase() );
        if ( idNameMap.get(id) != null )
            throw new Exception("Id "+id+" isn't unique.");
        idNameMap.put(    id, name );
        return symbol;
    }

    /**
     * Takes a valid Fortran identifier name or string and finds a similar
     * identifier name
     * that is not already used in fortranNameCol.  If name is already in
     * fortranNameCol then findAFortranName mutates name in a recognizable way
     * until it is unique among all the fortranNameCol rows or findAFortranName
     * has exhausted its mutations.  If the mutations are exhausted then an
     * exception is thrown.
     */
    public String findAFortranName( String name ) throws Exception {
        int index = 0;
        int id = ID_START;
        String suffix = "";
        String newName = "";
        if ( lookupFortranNameIndex( name ) < 0 ) return name;
        // Append or replace end with _id.
        suffix = Integer.toString(id);
        while ( suffix.length() < ID_WIDTH ) { suffix = "0" + suffix; }
        suffix = "_" + suffix;
        if ( name.length() < maxLength - suffix.length() ) {
            newName = name + suffix;
        } else {
            newName =
                name.substring( 0, maxLength - suffix.length() )
                + suffix
            ;
        }
        // If name still in fortranNameCol then increment id and replace old id
        // in name with new.  Only repeat until id = ID_END.
        while (
            ( index = lookupFortranNameIndex( newName ) ) >= 0
            && id < ID_END
        ) {
            id++;
            suffix = Integer.toString(id);
            while ( suffix.length() < ID_WIDTH ) { suffix = "0" + suffix; }
            suffix = "_" + suffix;
            if ( name.length() < maxLength - suffix.length() ) {
                newName = name + suffix;
            } else {
                newName =
                    name.substring( 0, maxLength - suffix.length() )
                    + suffix
                ;
            }
        }
        if ( index >= 0 )
            throw new Exception( "A unique Fortran name couldn't be found." );
        return newName;
    }

    /**
     * Returns a Fortran safe string according to <i>mode</i>
     */
    public String makeFortranSafe( String string )
        throws Exception
    {
        return makeFortranSafe(
            string, mode, maxLength, padToLength, padCharacter
        );
    }

    /**
     * Makes string Fortran safe according to <i>mode</i>.
     */
    public String makeFortranSafe( String string, Mode mode )
        throws Exception
    {
        return makeFortranSafe(
            string, mode, maxLength, padToLength, padCharacter
        );
    }

    public String makeFortranSafe(
        String string, Mode mode, int maxLength
    ) throws Exception {
        return makeFortranSafe(
            string, mode, maxLength, padToLength, padCharacter
        );
    }

    /**
     * Makes string Fortran safe according to <i>mode</i> and <i>maxLength</i>.
     */
    public String makeFortranSafe(
        String string, Mode mode, int maxLength, boolean padToLength,
        char padCharacter
    ) throws Exception {
        if        ( mode == Mode.IDENTITY ) {
            return string;
        } else if (
            mode == Mode.FORTRAN_IDENTIFIER ||
            mode == Mode.FORTRAN_IDENTIFIER_UNDERSCORED
        ) {
            return makeFortranIdSafe_( string, maxLength );
        } else if ( mode == Mode.FORTRAN_CHARACTER_STRING ) {
            return makeFortranCharacterSafe_(
                string, maxLength, padToLength, padCharacter )
            ;
        } else throw new Exception("Unrecognized mode.");
    }

    /**
     * Makes a string safe as a Fortran character string.  Uses the same
     * algorithm as {@link #makeFortranIdSafe_(String,int)} with the following
     * differences:
     * <ul>
     * <li>No string is ever prepended.
     * <li>The characters [0-9A-Za-z~`!#$%&amp;*()_+-=[]\{};':",./&lt;&gt;? ]
     * are permitted anywhere in the string.
     * <li>The characters [@^|] are converted to _u??, u??_, or _u??_ depending
     * on whether the character is at the end, beginning, or middle of the
     * string, respectively.
     * <li>All other characters are converted to _u??, u??_, or _u??_ depending
     * on whether the character is at the end, beginning, or middle of the
     * string, respectively.
     * </ul>
     */
    public String makeFortranCharacterSafe_(
        String string, int maxLength, boolean padToLength, char padCharacter
    ) throws Exception {
        StringBuffer out = new StringBuffer();
        StringBuffer in  = new StringBuffer(string);
        while ( in.length() > 0 && out.length() < maxLength ) {
            if ( in.toString().matches(
                "^[0-9A-Za-z~`!#$%&*()_+\\-=\\[\\]\\\\{};':\",./<>? ].*")
            ) {
                out.append(in.charAt(0));
            } else {
                if ( out.length() > 1 ) out.append("_");
                out.append("u"+Integer.toHexString((int)in.charAt(0)));
                if ( in.length() > 1 ) out.append("_");
            }
            in.deleteCharAt(0);
        }
        if ( out.length() > maxLength ) {
            return out.substring(0,maxLength);
        }
        if ( padToLength )
            while ( out.length() < maxLength )
                out.append(padCharacter);
        return out.toString();
    }

    /**
     * Convert a string to a Fortran safe identifier.
     * <p><b>Documentation out of date.</b><p>
     * Fortran 95 only accepts
     * [A-Za-z0-9_] in identifier names.  Addtionally an identifier must begin
     * with a letter and be less than 31 characters.  Identifiers can be
     * Fortran keywords.
     * <p>
     * The algorithm works as follows:
     <pre>
     If the first character is not a letter prepend PREFIX.
     Convert series of ['"]* to _ppp_ if in the middle of the string and _ppp if
         at the end, where the number of ps is the number of
         apostrophes (each ' counts as one p and each " counts as 2).
     Convert all other non alphanumeric characters to a number as they are
         defined in Java.  Surround the number with "_" if in the middle of the
         string and prepend with "_" if at the end.  Underscores are also
         converted to numbers.  This allows symbols from the symbol table to
         take suffixes of the form "__MYUSERSUFFIX" and still guarantee
         uniqueness.  The user must take care not to use suffixes that might be
         generated by this routine (i.e. _p, and _u34).  Also, the user must
         take care not to exceed the Fortran limit on variable names.  This
         routine allows the user to restrict var name length to ensure suffixes
         can be appended.
     </pre>
     */
    public String makeFortranIdSafe( String string ) throws Exception {
        return makeFortranIdSafe_( string, maxLength );
    }

   private String makeFortranIdSafe_ (String in, int maxLength) {
      int length = in.length ();
      if (length > maxLength)
         length = maxLength;
      if (length == 0)
         return "";
      StringBuffer out = new StringBuffer ();
      char c = in.charAt (0);
      int count = 0;
      int pos = 1;

      // Avoid identifiers starting with non-alphabetic characters
      if ((c < 'a' || c > 'z') && (c < 'A' || c > 'Z')) {
         out.append ("A_");
         count += 2;
      }

outer:
      while (count < maxLength) {
         if (c == '_' || c >= 'a' && c <= 'z' || c >= 'A' && c <= 'Z' || c >= '0' && c <= '9') {
            out.append (c);
            count++;
         } else if (c == '\'' || c == '"') {
            out.append ("_");
            count++;
            while (true) {
               if (c == '\'') {
                  out.append ("p");
                  count++;
               } else if (c == '"') {
                  out.append ("pp");
                  count += 2;
               } else
                  continue outer;
               if (pos == length)
                  break outer;
               c = in.charAt (pos++);
            }
         } else {
            out.append ("_u");
            String hexCode = Integer.toHexString ((int) c);
            out.append (hexCode);
            count += 2 + hexCode.length ();
         }
         if (pos == length)
            break;
         c = in.charAt (pos++);
      }
      String temp = 
         count > maxLength-1 ? out.substring (0, maxLength-1) : out.toString ();
      return temp+(mode.appendUnderscore()?"_":"");
   }
}
