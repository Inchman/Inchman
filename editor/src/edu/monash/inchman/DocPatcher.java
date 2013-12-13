package edu.monash.inchman;

import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Iterator;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.parsers.FactoryConfigurationError;
import javax.xml.stream.XMLStreamException;
import javax.xml.xpath.XPathExpressionException;

import jigcell.sbml2.Annotations;
import jigcell.sbml2.Compartment;
import jigcell.sbml2.KineticLaw;
import jigcell.sbml2.Model;
import jigcell.sbml2.Parameter;
import jigcell.sbml2.Reaction;
import jigcell.sbml2.SBMLLevel2Document;
import jigcell.sbml2.SBase;
import jigcell.sbml2.SBaseId;
import jigcell.sbml2.Species;
import jigcell.sbml2.SpeciesReference;

import com.ociweb.xml.ElementWAX;
import com.ociweb.xml.WAX;

import edu.monash.inchman.ImModel.ImCompartment;
import edu.monash.inchman.ImModel.ImMath;
import edu.monash.inchman.ImModel.ImParameter;
import edu.monash.inchman.ImModel.ImReaction;
import edu.monash.inchman.ImModel.ImReactionCompartmentIsAllowed;
import edu.monash.inchman.ImModel.ImReactionSpeciesStoichiometry;
import edu.monash.inchman.ImModel.ImSpecies;
import edu.monash.inchman.ImModel.ImSpeciesCompartmentInitialAmount;

import com.google.gson.Gson; // todo: only for debug purpose .. remove

public abstract class DocPatcher {
	
	public static final String SolverKey         = "solver"; 
    public static final String RunsKey           = "runs";
    //public static final String StepsKey          = "steps";
    public static final String TimeKey           = "time";
    public static final String OutputIntervalKey = "outputInterval";
    
    public static final String GridWidthKey      = "gridWidth";
    public static final String GridHeightKey     = "gridHeight";
    public static final String PhysicalWidthKey  = "physicalWidth";
    public static final String PhysicalHeightKey = "physicalHeight";
    
    public static final String InchmanPrefix = "im";
    public static final String InchmanNS     = "http://www.csse.monash.edu.au/~berndm/inchman";
    
    /*
     *  java.io.FileWriter is not supported by Google App Engine's Java runtime environment
    public static void patchSbmlDocumentToFile(ImModel data, String fileName) throws XMLStreamException, FactoryConfigurationError, XPathExpressionException, IOException
    {
    	Writer w = new BufferedWriter(new FileWriter(fileName));
    	w.write(patchSbmlDocumentToString(data));
    	w.close(); // this proved critical in practice!!! (files were incomplete)
    }
    */
    
    public static String patchSbmlDocumentToString(ImModel data) throws XMLStreamException, FactoryConfigurationError, XPathExpressionException, IOException
    {
    	SBMLLevel2Document sbml = patchSbmlDocumentToDocument(data);
    	
    	// Try not to clobber the original sbml tag - keep the original
    	// Thus we'll preserve extra attributes like 'xmlns:celldesigner="http://www.sbml.org/2001/ns/celldesigner"'
		StringWriter writer = new StringWriter();
		sbml.writeDocument(writer); // write entire document, including xml header and sbml element, annotations, notes, etc...
		writer.flush();
		
		String newDoc = writer.toString();
		
		if (data.origSbml != null) {
		    Pattern p = Pattern.compile("(<[sS][bB][mM][lL][^>]*?>)"); // note that "*?" is a "reluctant" search (which prevents it from matching too much)
    		Matcher origDocMatcher = p.matcher(data.origSbml);
    		if (origDocMatcher.find()) {
    			String quotedReplacement = Matcher.quoteReplacement(origDocMatcher.group());
    			newDoc = p.matcher(newDoc).replaceFirst(quotedReplacement);
    		}
		}
		
		return newDoc;
    }
    
    @SuppressWarnings("unchecked")
    private static SBMLLevel2Document patchSbmlDocumentToDocument(ImModel data) throws XMLStreamException, FactoryConfigurationError, XPathExpressionException, IOException
    {
        SBMLLevel2Document sbml = ((data.origSbml == null
                                    || data.origSbml.isEmpty())
                                   ? new SBMLLevel2Document(new Model(data.options.name))
                                   : SBMLLevel2Document.readDocument(new StringReader(data.origSbml)));
        Model model = sbml.getModel();
        
        // Update the MODEL name
        model.setName(data.options.name);
        
        // SBML root
        InchmanWax sbmlWax = new InchmanWax(sbml, "sbml");
        sbmlWax.child(SolverKey, data.options.solver);
        sbmlWax.unescapedChild(RunsKey, data.options.runs.cMML);
        //sbmlWax.unescapedChild(StepsKey, data.options.steps.cMML);
        sbmlWax.unescapedChild(TimeKey, data.options.time.cMML);
        sbmlWax.unescapedChild(OutputIntervalKey, data.options.outputInterval.cMML);
        {
            // Initialization script
            sbmlWax.start("init");
            sbmlWax.attr("type", data.init.type);
            {
                sbmlWax.start("script");
                sbmlWax.wax.cdata(data.init.script);
                sbmlWax.end();
            }
            sbmlWax.end();

            // Events script and time-stamps
            sbmlWax.start("events");
            sbmlWax.attr("type", data.events.type);
            {
                sbmlWax.start("script");
                sbmlWax.wax.cdata(data.events.script);
                sbmlWax.end();
                
                sbmlWax.start("listOfTimeStamps");
                for (ImMath timeStamp : data.events.timeStamps) {
                    sbmlWax.unescapedChild("timeStamp", timeStamp.cMML);
                }
                sbmlWax.end();
            }
            sbmlWax.end();
            
            // Drift Diffusivity
            sbmlWax.start("driftDiffusivity");
            sbmlWax.attr("nonLinear", data.driftDiffusivity.nonLinear);
            sbmlWax.attr("computeMoments", data.driftDiffusivity.computeMoments);
            {
                sbmlWax.start("method");
                sbmlWax.wax.cdata(data.driftDiffusivity.method);
                sbmlWax.end();
            }
            sbmlWax.end();
            
            // Update Fields method
            sbmlWax.start("setField");
            {
            	sbmlWax.start("method");
            	sbmlWax.wax.cdata(data.updateFields.method);
            	sbmlWax.end();
            }
            sbmlWax.end();
            
            // New individuals method
            sbmlWax.start("newIndividualsMethod");
            {
            	sbmlWax.start("method");
            	sbmlWax.wax.cdata(data.newIndividualsMethod.method);
            	sbmlWax.end();
            }
            sbmlWax.end();
            
        }
        sbmlWax.close();
        
        // Model root
        InchmanWax modelWax = new InchmanWax(model, "model");
        modelWax.unescapedChild(GridWidthKey, data.options.gridWidth.cMML);
        modelWax.unescapedChild(GridHeightKey, data.options.gridHeight.cMML);
        modelWax.unescapedChild(PhysicalWidthKey, data.options.physicalWidth.cMML);
        modelWax.unescapedChild(PhysicalHeightKey, data.options.physicalHeight.cMML);
        modelWax.close();
        
        //
        // Parameters
        //
        Map<String, Parameter> sbmlParametersRemaining = collectionListToHash(model.getParameters());
        // Add and update
        for (ImParameter imParam : data.parameters)
        {
            Parameter sbmlParam = sbmlParametersRemaining.get(imParam.sbmlId);
            if (sbmlParam == null) {
                sbmlParam = new Parameter(imParam.sbmlId, imParam.name);
                sbmlParam.setConstant(true); // constant for the life-time of the simulation - although it vary from one simulation to another
                model.addParameter(sbmlParam);
            } else {
                // sbmlParam.setId(imParam.sbmlId); -- if they were different then sbmlSpecies would be null and this code wouldn't be called anyway
                sbmlParam.setName(imParam.name);
            }
            
            // Add some sanity!
            if (!imParam.data.from.isEmpty()) {
                sbmlParam.setValue(Double.parseDouble(imParam.data.from));
            }
            
            InchmanWax paramWax = new InchmanWax(sbmlParam, "parameter");
            paramWax.attr("type", imParam.data.type);
            paramWax.attr("domain", imParam.data.domain);
            paramWax.attr("from", imParam.data.from);
            paramWax.attr("to", imParam.data.to);
            paramWax.attr("step", imParam.data.step);
            paramWax.attr("points", imParam.data.points);
            paramWax.close();
            
            sbmlParametersRemaining.remove(imParam.sbmlId);
        }
        // Remove parameters that no longer exist in the Inchman version
        model.getParameters().removeAll(sbmlParametersRemaining.values()); // remove the Parameter objects (the values) from the List (non-keyed)

        
        //
        // Compartments
        //
        Map<String, Compartment> sbmlCompartmentsRemaining = collectionListToHash(model.getCompartments());
        // Add and update/patch
        for (ImCompartment imCompart : data.compartments)
        {
            Compartment sbmlCompart = sbmlCompartmentsRemaining.get(imCompart.sbmlId);
            if (sbmlCompart == null) {
                sbmlCompart = new Compartment(imCompart.sbmlId, imCompart.name);
                model.addCompartment(sbmlCompart);
            } else {
                // sbmlCompart.setId(imCompart.sbmlId); -- if they were different then sbmlSpecies would be null and this code wouldn't be called anyway
                sbmlCompart.setName(imCompart.name);
            }
            
            InchmanWax compartWax = new InchmanWax(sbmlCompart, "compartment");
            compartWax.unescapedChild("x", imCompart.x.cMML);
            compartWax.unescapedChild("y", imCompart.y.cMML);
            compartWax.unescapedChild("width", imCompart.width.cMML);
            compartWax.unescapedChild("height", imCompart.height.cMML);
            compartWax.close();
            
            sbmlCompartmentsRemaining.remove(imCompart.sbmlId);
        }
        // Remove compartments that no longer exist in the Inchman version
        model.getCompartments().removeAll(sbmlCompartmentsRemaining.values()); // remove the Compartment objects (the values) from the List (non-keyed)
        
        
        //
        // Species
        //
        Map<String, Species> sbmlSpeciesRemaining = collectionListToHash(model.getSpecies());
        // Add and update/patch
        for (ImSpecies imSpecies : data.species)
        {
            Species sbmlSpecies = sbmlSpeciesRemaining.get(imSpecies.sbmlId);
            if (sbmlSpecies == null) {
                sbmlSpecies = new Species(imSpecies.sbmlId, imSpecies.name);
                // we need to set "World" as the standard compartment to satisfy SBML validity requirement
                // (every species needs to be associated with a compartment)
                // Inchman does not support species to only live in certain compartments so they all live in "World"
                sbmlSpecies.setCompartment("World");
                
                model.addSpecies(sbmlSpecies);
            } else {
                // sbmlSpecies.setId(imSpecies.sbmlId); -- if they were different then sbmlSpecies would be null and this code wouldn't be called anyway
                sbmlSpecies.setName(imSpecies.name);
            }

            // Generate the new inchman annotation
            InchmanWax speciesWax = new InchmanWax(sbmlSpecies, "species");
            speciesWax.unescapedChild("diffusionConstant", imSpecies.diffusionConstant.cMML);
            // check if it's an individual species
            if (imSpecies.individual) {
            	speciesWax.start("individual");
            	speciesWax.end();
            }
            speciesWax.start("listOfCompartmentParameters");
            {
                for (ImSpeciesCompartmentInitialAmount compartmentInitialAmount : imSpecies.compartmentInitialAmounts) {
                    speciesWax.start("compartmentParameters");
                    speciesWax.attr("compartment", compartmentInitialAmount.compartmentSbmlId); // simple value- fine for it to be an attribute
                    speciesWax.unescapedChild("initialAmount", compartmentInitialAmount.initialAmount.cMML);
                    speciesWax.end(); // compartmentParameters
                }
            }
            speciesWax.close();
            
            sbmlSpeciesRemaining.remove(imSpecies.sbmlId);
        }
        // Remove species that no longer exist in the Inchman version
        model.getSpecies().removeAll(sbmlSpeciesRemaining.values()); // remove the Species objects (the values) from the List (non-keyed)
        
        
        //
        // Reactions
        //
        Map<String, Reaction> sbmlReactionsRemaining = collectionListToHash(model.getReactions());
        // Add and update/patch
        for (ImReaction imReaction : data.reactions)
        {
            Reaction sbmlReaction = sbmlReactionsRemaining.get(imReaction.sbmlId);
            if (sbmlReaction == null) {
                sbmlReaction = new Reaction(imReaction.sbmlId, imReaction.name);
                model.addReaction(sbmlReaction);
            } else {
                // sbmlReaction.setId(imReaction.sbmlId); -- if they were different then sbmlReaction would be null and this code wouldn't be called anyway
                sbmlReaction.setName(imReaction.name);
            }
            
            
            // Kinetic Law
            KineticLaw kl = new KineticLaw();
            kl.setMath(imReaction.kineticLaw.cMML);
            sbmlReaction.setKineticLaw(kl);
            
            // Reactants
            // future: factorize this with products
            // future: support sbml's stoichiometry math
            
            // these ones are the reactants which were defined in the original (unpatched) SBML
            // we now need to compare to the actual model
            Map<String, SpeciesReference> sbmlReactantsRemaining = speciesReferenceListToHash(sbmlReaction.getReactant());
            // Add and update/patch
            for (ImReactionSpeciesStoichiometry imReactant : imReaction.reactants)
            {            	
                SpeciesReference sbmlReactant = sbmlReactantsRemaining.get(imReactant.speciesSbmlId);
                if (sbmlReactant == null) {
                    sbmlReactant = new SpeciesReference();
                    sbmlReactant.setSpecies(imReactant.speciesSbmlId);
                    sbmlReaction.addReactant(sbmlReactant);
                }
                
                // we set the stoichiometry
                sbmlReactant.setStoichiometry(Double.parseDouble(imReactant.stoichiometry));

                // and remove it from the existing list..
                // we can do that since we now that species in the inchman model only exist once (they're clean)
                sbmlReactantsRemaining.remove(imReactant.speciesSbmlId);
            }
            // Remove reactants that no longer exist in the Inchman version
            //System.out.println(new Gson().toJson(imReaction.reactants));
            //System.out.println(new Gson().toJson(sbmlReaction.getReactant()));
            sbmlReaction.getReactant().removeAll(sbmlReactantsRemaining.values()); // remove the SpeciesReference objects (the values) from the List (non-keyed)
            
            // Products
            // future: factorize this with reactants
            // future: support sbml's stoichiometry math
            Map<String, SpeciesReference> sbmlProductsRemaining = speciesReferenceListToHash(sbmlReaction.getProduct());
            // Add and update/patch
            for (ImReactionSpeciesStoichiometry imProduct : imReaction.products)
            {
                SpeciesReference sbmlProduct = sbmlProductsRemaining.get(imProduct.speciesSbmlId);
                if (sbmlProduct == null) {
                    sbmlProduct = new SpeciesReference();
                    sbmlProduct.setSpecies(imProduct.speciesSbmlId);
                    sbmlReaction.addProduct(sbmlProduct);
                }
                sbmlProduct.setStoichiometry(Double.parseDouble(imProduct.stoichiometry));
                
                sbmlProductsRemaining.remove(imProduct.speciesSbmlId);
            }
            // Remove species that no longer exist in the Inchman version
            sbmlReaction.getProduct().removeAll(sbmlProductsRemaining.values()); // remove the SpeciesReference objects (the values) from the List (non-keyed)
            

            // Generate the new inchman annotation
            InchmanWax reactionWax = new InchmanWax(sbmlReaction, "reaction");
            reactionWax.start("listOfCompartmentParameters");
            {
                for (ImReactionCompartmentIsAllowed compartmentIsAllowed : imReaction.compartmentIsAlloweds) {
                    reactionWax.start("compartmentParameters");
                    reactionWax.attr("compartment", compartmentIsAllowed.compartmentSbmlId); // simple value- fine for it to be an attribute
                    reactionWax.unescapedChild("isAllowed", compartmentIsAllowed.isAllowed.cMML);
                    reactionWax.end();
                }
            }
            reactionWax.close();
            
            sbmlReactionsRemaining.remove(imReaction.sbmlId);
        }
        // Remove reactions that no longer exist in the Inchman version
        model.getReactions().removeAll(sbmlReactionsRemaining.values()); // remove the Reaction objects (the values) from the List (non-keyed)

        return sbml;
    }
    
    
    // Note: this class is designed to be non-destructive - preserving any other system's annotations
    private static class InchmanWax {
        SBase obj = null;
        String objClassName = null;
        StringWriter stringWriter = new StringWriter();
        WAX wax = new WAX(stringWriter);
        
        InchmanWax(SBase obj, String objClassName) {
            this.obj = obj;
            this.objClassName = objClassName;
            this.start(objClassName);
            this.wax.namespace(InchmanPrefix, InchmanNS); // add this to the first element only
        }
        
        ElementWAX attr(String key, String value) {
            return this.wax.attr(key, value);
        }
        ElementWAX start(String key) {
            return this.wax.start(InchmanPrefix, key);
        }
        ElementWAX child(String key, String value) {
            return this.wax.child(InchmanPrefix, key, value);
        }
        ElementWAX unescapedChild(String key, String value) { // Handy for inserting existing XML data, such as MathML
            return this.wax.start(InchmanPrefix, key).unescapedText(value).end();
        }
        ElementWAX end() {
            return this.wax.end();
        }
        
        void close() throws XPathExpressionException {
            this.wax.close(); // end objClassName
            patchAnnotations(this.obj,
                             InchmanPrefix + ":" + this.objClassName,
                             this.stringWriter.toString());
        }
    };
    
    // Find and patch or add the annotation
    private static void patchAnnotations(SBase obj, String tag, String imAnnotation) throws XPathExpressionException
    {
        imAnnotation = "\n" + imAnnotation + "\n";
        
		Annotations as = obj.getAnnotations();
        if (as.size() <= 0) {
            // Add annotation
            as.add(imAnnotation);
        }
        else {
            if (as.size() > 1) {
                System.err.println("WARNING: Multiple annotation elements found for \"" + tag + "\"."
                                 + "         According to reference: L2V4 Section 4.1\n"
                                 + "         Only one <annotation> element is permitted inside a particular containing element."
                                 + "         Therefore, only the first annotation elemenet will be used.");
            }
            
        	// Don't even try to parse foreign annotations with xpath!
        	// Especially since software like CellDesigner (celldesigner.org) put the xmlns prefix
        	// in the root SBML element and hence parsing the annotation alone will be problematic.
            String existingAnnotation = as.get(0);
        	Pattern imPattern = Pattern.compile("(<" + tag + /* DON'T CLOSE THE ELEMENT with '>" - read through the xmlns prefix */ ".+?</" + tag + "\\w*?>)"); // reluctant match of entire element
        	Matcher imMatcher = imPattern.matcher(existingAnnotation);
        	String patched = null;
        	if (imMatcher.find()) {
        	    // Replace existing inchman annotation
        	    patched = imMatcher.replaceFirst(Matcher.quoteReplacement(imAnnotation));
        	}
        	else {
        	    // Append to existing annotation element
        	    Pattern endPattern = Pattern.compile("</annotation\\w*?>");
        	    String newEnd = Matcher.quoteReplacement(imAnnotation + "</annotation>");
        	    patched = endPattern.matcher(existingAnnotation).replaceFirst(newEnd);
        	}
        	as.set(0, patched);
        }
    }
    
    
    @SuppressWarnings({ "unchecked", "rawtypes" })
    private static <T> Map<String,T> collectionListToHash(List list)
    {
        HashMap<String,T> map = new HashMap<String,T>();
        for (Object o : list) {
            T item = (T)o;
            map.put(((SBaseId)item).getId(), item);
        }
        return map;
    }

    private static Map<String,SpeciesReference> speciesReferenceListToHash(List<SpeciesReference> list)
    {
        HashMap<String,SpeciesReference> map = new HashMap<String,SpeciesReference>();
        //for (SpeciesReference item : list) {
        // need to use iterators so we can safely remove
        for (Iterator<SpeciesReference> it = list.iterator(); it.hasNext();) {
        	SpeciesReference item = it.next();
        	
        	// we need to check if the species is already in there .. in which case we just increase the stoichiometry
        	if (map.containsKey(item.getSpecies())) {
        		SpeciesReference sbmlReactant = map.get(item.getSpecies());
        		// need to increase by the stoichiometry of the new item
        		sbmlReactant.setStoichiometry(sbmlReactant.getStoichiometry()+item.getStoichiometry());
        		// and remove the old one to avoid duplicates
        		it.remove();
        	} else
        		map.put(item.getSpecies(), item);
        }
        return map;
    }
    

}
