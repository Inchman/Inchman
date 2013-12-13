// Created by Aidan Lane, Thu Jun 23, 2011

package edu.monash.inchman;

import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.UUID;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.transform.stream.StreamSource;

import jigcell.sbml2.Annotations;
import jigcell.sbml2.Compartment;
import jigcell.sbml2.Model;
import jigcell.sbml2.Parameter;
import jigcell.sbml2.Reaction;
import jigcell.sbml2.SBMLLevel2Document;
import jigcell.sbml2.SBase;
import jigcell.sbml2.Species;
import jigcell.sbml2.SpeciesReference;
import net.sf.saxon.s9api.Processor;
import net.sf.saxon.s9api.QName;
import net.sf.saxon.s9api.SaxonApiException;
import net.sf.saxon.s9api.XPathCompiler;
import net.sf.saxon.s9api.XPathSelector;
import net.sf.saxon.s9api.XdmItem;
import net.sf.saxon.s9api.XdmNode;
import net.sf.saxon.s9api.XdmSequenceIterator;
import net.sf.saxon.s9api.XdmValue;
import edu.monash.inchman.ImModel;
import edu.monash.inchman.ImModel.ImCompartment;
import edu.monash.inchman.ImModel.ImDriftDiffusivity;
import edu.monash.inchman.ImModel.ImEvents;
import edu.monash.inchman.ImModel.ImInit;
import edu.monash.inchman.ImModel.ImMath;
import edu.monash.inchman.ImModel.ImOptions;
import edu.monash.inchman.ImModel.ImParameter;
import edu.monash.inchman.ImModel.ImParameterData;
import edu.monash.inchman.ImModel.ImReaction;
import edu.monash.inchman.ImModel.ImReactionCompartmentIsAllowed;
import edu.monash.inchman.ImModel.ImReactionSpeciesStoichiometry;
import edu.monash.inchman.ImModel.ImSpecies;
import edu.monash.inchman.ImModel.ImSpeciesCompartmentInitialAmount;
import edu.monash.inchman.ImModel.ImUpdateFields;
import edu.monash.inchman.DocPatcher;

public class DocLoader {
    
    //static final ImMath math_value_neg1  = encodeMathFromDouble(-1.0);
    static final ImMath math_value_0     = encodeMathFromDouble(0.0);
    static final ImMath math_value_1     = encodeMathFromDouble(1.0);
    static final ImMath math_value_10    = encodeMathFromDouble(10.0);
    static final ImMath math_value_32    = encodeMathFromDouble(32.0);
    static final ImMath math_value_40    = encodeMathFromDouble(40.0);
    static final ImMath math_value_100   = encodeMathFromDouble(100.0);
    static final ImMath math_value_10000 = encodeMathFromDouble(10000.0);
    static final ImMath math_value_false = encodeMathFromIdentifier("false");
    static final ImMath math_value_true  = encodeMathFromIdentifier("true");
    
    
    static Processor globalProcessor = new Processor(false);
    
    
	public static ImModel loadSbmlDocument(String sbmlContents) throws Exception
    {
		ImModel imModel = new ImModel();
	    
	    // TODO: remove the need for this WORKAROUND...
	    // XPath is choking on complex python code on the CDATA section
	    //System.out.println("before:" + sbmlContents);
	    HashMap<String, String> cdataHash = new HashMap<String, String>();
	    StringBuilder sbmlStringBuilder = new StringBuilder();
	    Matcher m = Pattern.compile("\\Q<![CDATA[\\E"
	                                + "(.*?)" // reluctant quantifiers - prevent grabbing multiple CDATA sections as a single big one
	                                + "\\Q]]>\\E",
	                                Pattern.DOTALL).matcher(sbmlContents);
	    int prevIndex = 0;
	    while (m.find()) {
	        sbmlStringBuilder.append(sbmlContents.substring(prevIndex, m.start())); // last char = endIndex-1 = m.start()-1
	        
	        String tag = UUID.randomUUID().toString();
            cdataHash.put(tag, m.group(1));
            sbmlStringBuilder.append(tag); // a placeholder
	        
	        prevIndex = m.end();
	    }
	    if (prevIndex < sbmlContents.length()-1) {
	        sbmlStringBuilder.append(sbmlContents.substring(prevIndex, sbmlContents.length())); // last char = endIndex-1 =  sbmlContents.length()-1
	    }
	    sbmlContents = sbmlStringBuilder.toString();
	    
        
        SBMLLevel2Document sbml = SBMLLevel2Document.readDocument(new StringReader(sbmlContents));
        Model model = sbml.getModel();
        
        imModel.origSbml         = sbmlContents;
        imModel.options          = readOptions(sbml);
        imModel.parameters       = readParameters(model);
        imModel.compartments     = readCompartments(model);
        imModel.species          = readSpecies(model);
        imModel.reactions        = readReactions(model);
        imModel.init             = readInit(sbml, cdataHash);
        imModel.newIndividualsMethod = readNewIndividualsMethod(sbml, cdataHash);
        imModel.events           = readEvents(sbml, cdataHash);
        imModel.driftDiffusivity = readDriftDiffusivity(sbml, cdataHash);
        imModel.updateFields	 = readUpdateFields(sbml, cdataHash);
        
        return imModel;
    }
	
	private static ImOptions readOptions(SBMLLevel2Document sbml) throws SaxonApiException, IOException {
        ImOptions options      = new ImOptions();
        Model     model        = sbml.getModel();  
        XdmItem   modelNode    = getXdmItemFromAnnotation(model, "//im:model");
        XdmItem   sbmlNode     = getXdmItemFromAnnotation(sbml,  "//im:sbml");
        options.name           = model.getName().isEmpty() ? "Untitled" : model.getName();
        options.gridWidth      = getMathFromXPath(modelNode, DocPatcher.GridWidthKey,      math_value_32);
        options.gridHeight     = getMathFromXPath(modelNode, DocPatcher.GridHeightKey,     math_value_32);
        options.physicalWidth  = getMathFromXPath(modelNode, DocPatcher.PhysicalWidthKey,  math_value_40);
        options.physicalHeight = getMathFromXPath(modelNode, DocPatcher.PhysicalHeightKey, math_value_40);
        options.solver         = getStringFromXPath(sbmlNode,DocPatcher.SolverKey,         "stochastic_homogeneous");
        options.runs           = getMathFromXPath(sbmlNode,  DocPatcher.RunsKey,           math_value_10);
        //options.steps          = getMathFromXPath(sbmlNode,  DocPatcher.StepsKey,          math_value_neg1); // -1 steps -> infinity
        options.time           = getMathFromXPath(sbmlNode,  DocPatcher.TimeKey,           math_value_100);
        options.outputInterval = getMathFromXPath(sbmlNode,  DocPatcher.OutputIntervalKey, math_value_10);
        return options;
	}
	
	private static ArrayList<ImParameter> readParameters(Model model) throws SaxonApiException {
	    ArrayList<ImParameter> parameters = new ArrayList<ImParameter>(model.getParameters().size());
        for (Object o : model.getParameters())
        {
            Parameter sbmlParam = (Parameter)o;
            XdmNode n = getXdmNodeFromAnnotation(sbmlParam, "//annotation/im:parameter");
            
            ImParameter imParam = new ImParameter();
            imParam.sbmlId      = sbmlParam.getId();
            imParam.name        = sbmlParam.getName();
            imParam.data        = new ImParameterData();
            imParam.data.type   = getNodeAttribute(n, "type",   "float"); // SBML uses a floating-point encoded numbers - so default to "float"
            imParam.data.domain = getNodeAttribute(n, "domain", "single");
            imParam.data.from   = getNodeAttribute(n, "from",   Double.toString(sbmlParam.getValue()));  // fallback to the SBML param value if "from" unavailable
            imParam.data.to     = getNodeAttribute(n, "to",     "0");
            imParam.data.step   = getNodeAttribute(n, "step",   "1");
            imParam.data.points = getNodeAttribute(n, "points", "1");
            
            parameters.add(imParam);
        }
        return parameters;
	}
	
	private static ArrayList<ImCompartment> readCompartments(Model model) throws SaxonApiException, IOException {
	    ArrayList<ImCompartment> compartments = new ArrayList<ImCompartment>(model.getCompartments().size());
        for (Object o : model.getCompartments())
        {
            Compartment sbmlCompart = (Compartment)o;
            XdmNode n = getXdmNodeFromAnnotation(sbmlCompart, "//im:compartment");
            
            ImCompartment imCompart = new ImCompartment();
            imCompart.sbmlId = sbmlCompart.getId();
            imCompart.name   = sbmlCompart.getName();
            imCompart.x      = getMathFromXPath(n, "x",      math_value_0);
            imCompart.y      = getMathFromXPath(n, "y",      math_value_0);
            imCompart.width  = getMathFromXPath(n, "width",  math_value_1);
            imCompart.height = getMathFromXPath(n, "height", math_value_1);
            
            compartments.add(imCompart);
        }
        return compartments;
	}
	
	private static ArrayList<ImSpecies> readSpecies(Model model) throws Exception {
	    ArrayList<ImSpecies> species = new ArrayList<ImSpecies>(model.getSpecies().size());
        XPathSelector speciesCompartmentParamsXPathSelector = compileInchmanXPathSelector("im:listOfCompartmentParameters/im:compartmentParameters");
        for (Object o : model.getSpecies())
        {
            Species sbmlSpecies = (Species)o;
            
            ImSpecies imSpecies = new ImSpecies();
            imSpecies.sbmlId            = sbmlSpecies.getId();
            imSpecies.name              = sbmlSpecies.getName();
            imSpecies.diffusionConstant = math_value_0; // unless told otherwise...
            
            // Inchman Extensions: read diffusionConstant and the list of compartment parameters from the inchman annotations
            XdmItem n = getXdmItemFromAnnotation(sbmlSpecies, "//im:species");
            if (n != null)
            {
                imSpecies.diffusionConstant = getMathFromXPath(n, "diffusionConstant", imSpecies.diffusionConstant);
                
                // are they individuals?
                XdmValue individuals = (XdmValue) getXdmFromAnnotation(sbmlSpecies, "//im:species/im:individual", true);
                if (individuals != null) {
                	imSpecies.individual = true;
                }
                else
                {
                	imSpecies.individual = false;
                }
                
                speciesCompartmentParamsXPathSelector.setContextItem(n);
                XdmValue list = speciesCompartmentParamsXPathSelector.evaluate();
                if (list != null) {
                    imSpecies.compartmentInitialAmounts = new ArrayList<ImSpeciesCompartmentInitialAmount>(list.size());
                    XdmSequenceIterator i = list.iterator();
                    while (i.hasNext()) {
                    	XdmItem cpItem = i.next();
                    	if (cpItem.isAtomicValue())
                    		continue;
                    	XdmNode cp = (XdmNode)cpItem;
                    	String compartmentIdStr = cp.getAttributeValue(new QName("compartment"));
                        if (compartmentIdStr != null
                            /*&& compartments.containsKey(compartmentIdStr)*/) // TODO: re-enable this check?
                        {
                            ImSpeciesCompartmentInitialAmount ia = new  ImSpeciesCompartmentInitialAmount();
                            ia.compartmentSbmlId = compartmentIdStr;
                            ia.initialAmount     = getMathFromXPath(cp, "initialAmount", math_value_0);
                            imSpecies.compartmentInitialAmounts.add(ia);
                        }
                    }
                }
            }
            
            species.add(imSpecies);
        }
        return species;
	}
	
	private static ArrayList<ImReaction> readReactions(Model model) throws Exception {
        ArrayList<ImReaction> reactions = new ArrayList<ImReaction>(model.getReactions().size());
        for (Object o : model.getReactions())
        {
            Reaction sbmlReaction = (Reaction)o;
            
            ImReaction imReaction = new ImReaction();
            imReaction.sbmlId     = sbmlReaction.getId();
            imReaction.name       = sbmlReaction.getName();
            imReaction.kineticLaw = encodeMathFromCMML(sbmlReaction.getKineticLaw().getMath().replaceAll("math:", "").replace("<math>", "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">")); // TODO: remove the need for this hack!
            
            // Load Reactants
            // future: support sbml's stoichiometry math
            imReaction.reactants = new ArrayList<ImReactionSpeciesStoichiometry>(sbmlReaction.getReactant().size());
            for (Object reactant : sbmlReaction.getReactant()) {
                ImReactionSpeciesStoichiometry ss = new ImReactionSpeciesStoichiometry();
                ss.speciesSbmlId = ((SpeciesReference)reactant).getSpecies();
                ss.stoichiometry = Double.toString(((SpeciesReference)reactant).getStoichiometry());
                imReaction.reactants.add(ss);
            }
            
            // Load Products
            // future: support sbml's stoichiometry math
            imReaction.products = new ArrayList<ImReactionSpeciesStoichiometry>(sbmlReaction.getProduct().size());
            for (Object product : sbmlReaction.getProduct()) {
                ImReactionSpeciesStoichiometry ss = new ImReactionSpeciesStoichiometry();
                ss.speciesSbmlId = ((SpeciesReference)product).getSpecies();
                ss.stoichiometry = Double.toString(((SpeciesReference)product).getStoichiometry());
                imReaction.products.add(ss);
            }
            
            // Inchman Extensions: read diffusionConstant and the list of compartment parameters from the inchman annotations
            XdmValue list = getXdmValueFromAnnotation(sbmlReaction, "//im:reaction/im:listOfCompartmentParameters/im:compartmentParameters");
            if (list != null) {
                imReaction.compartmentIsAlloweds = new ArrayList<ImReactionCompartmentIsAllowed>(list.size());
                XdmSequenceIterator i = list.iterator();
                while (i.hasNext()) {
                    XdmItem cpItem = i.next();
                    if (cpItem.isAtomicValue())
                    	continue;
                	XdmNode cp = (XdmNode)cpItem;
                    String compartmentIdStr = cp.getAttributeValue(new QName("compartment"));
                    if (compartmentIdStr != null
                        /*&& compartments.containsKey(compartmentIdStr)*/) { // TODO: re-enable this check?
                        ImReactionCompartmentIsAllowed ca = new ImReactionCompartmentIsAllowed();
                        ca.compartmentSbmlId = compartmentIdStr;
                        ca.isAllowed = getMathFromXPath(cp, "isAllowed", math_value_0);
                        imReaction.compartmentIsAlloweds.add(ca);
                    }
                }
            }
            
            reactions.add(imReaction);
        }
        
        return reactions;
	}
	
	private static ImInit readInit(SBMLLevel2Document sbml, HashMap<String, String> cdataHash) throws SaxonApiException, IOException {
	    ImInit init = new ImInit();
        XdmNode initNode = getXdmNodeFromAnnotation(sbml, "//im:sbml/im:init");
        init.type   = getNodeAttribute(initNode, "type", "python");
        init.script = cdataHash.get(getStringFromXPath(initNode, "script", ""));
        if (init.script == null) init.script = ""; // if the hash table lookup failed then set to empty string (preventing GSON from returning "undefined") 
        return init;
    }

	private static ImEvents readEvents(SBMLLevel2Document sbml, HashMap<String, String> cdataHash) throws SaxonApiException, IOException {
	    ImEvents events = new ImEvents();
	    XdmNode eventsNode = getXdmNodeFromAnnotation(sbml, "//im:sbml/im:events");
        
        events.type   = getNodeAttribute(eventsNode, "type", "python");
        events.script = cdataHash.get(getStringFromXPath(eventsNode, "script", ""));
        
        // Read the timestamps
        if (eventsNode != null) {
            XdmValue list = evalXPathList(eventsNode, "//im:listOfTimeStamps/im:timeStamp");
            if (list != null) {
                events.timeStamps = new ArrayList<ImMath>(list.size());
                XdmSequenceIterator i = list.iterator();
                while (i.hasNext()) {
                	XdmItem item = i.next();
                	if (!item.isAtomicValue())
                		events.timeStamps.add(getMathFromNode((XdmNode)i.next(), math_value_0));
                }
            }
        }
        
        return events;
    }
	
	private static ImDriftDiffusivity readDriftDiffusivity(SBMLLevel2Document sbml, HashMap<String, String> cdataHash) throws SaxonApiException {
	    ImDriftDiffusivity driftDiffusivity = new ImDriftDiffusivity();
	    XdmNode driftDiffusivityNode    = getXdmNodeFromAnnotation(sbml, "//im:sbml/im:driftDiffusivity");
        driftDiffusivity.nonLinear      = getNodeAttribute(driftDiffusivityNode, "nonLinear",      "true");
        driftDiffusivity.computeMoments = getNodeAttribute(driftDiffusivityNode, "computeMoments", "false");
        driftDiffusivity.method         = cdataHash.get(getStringFromXPath(driftDiffusivityNode, "method", ""));
        if (driftDiffusivity.method == null) driftDiffusivity.method = ""; // if the hash table lookup failed then set to empty string (preventing GSON from returning "undefined") 
        return driftDiffusivity;
    }
	
	private static ImUpdateFields readUpdateFields(SBMLLevel2Document sbml, HashMap<String, String> cdataHash) throws SaxonApiException
	{
		ImUpdateFields updateFields = new ImUpdateFields();
		XdmNode updateFieldsNode    = getXdmNodeFromAnnotation(sbml, "//im:sbml/im:setField");
		
		if (updateFieldsNode != null) {
			updateFields.method = cdataHash.get(getStringFromXPath(updateFieldsNode, "method", ""));
			if (updateFields.method == null)
				updateFields.method = "";
		} else {
			updateFields.method = "";
		}
		
		return updateFields;
	}

	private static ImModel.ImNewIndividualsMethod readNewIndividualsMethod(SBMLLevel2Document sbml, HashMap<String, String> cdataHash) throws SaxonApiException
	{
		ImModel.ImNewIndividualsMethod newIndividualsMethod = new ImModel.ImNewIndividualsMethod();
		XdmNode newIndividualsMethodNode    = getXdmNodeFromAnnotation(sbml, "//im:sbml/im:newIndividualsMethod");
		
		if (newIndividualsMethodNode != null) {
			newIndividualsMethod.method = cdataHash.get(getStringFromXPath(newIndividualsMethodNode, "method", ""));
			if (newIndividualsMethod.method == null)
				newIndividualsMethod.method = "";
		} else {
			newIndividualsMethod.method = "";
		}
		
		return newIndividualsMethod;
	}
	
	private static XdmValue getXdmValueFromAnnotation(SBase obj, String xpath) throws SaxonApiException {
		return (XdmValue) getXdmFromAnnotation(obj, xpath, false);
	}
	
	private static XdmItem getXdmItemFromAnnotation(SBase obj, String xpath) throws SaxonApiException {
		return (XdmItem) getXdmFromAnnotation(obj, xpath, true);
	}
	
	private static XdmNode getXdmNodeFromAnnotation(SBase obj, String xpath) throws SaxonApiException {
		XdmItem item = getXdmItemFromAnnotation(obj, xpath);
		return (item != null && !item.isAtomicValue()) ? (XdmNode)item : null;
	}
	
	private static Object getXdmFromAnnotation(SBase obj, String xpath, boolean evaluateSingle) throws SaxonApiException
    {
		// FUTURE: introduce some caching here
		Annotations as = obj.getAnnotations();
		if (obj.getAnnotations().size() > 0)
		{
			XPathSelector xpathSelector = compileInchmanXPathSelector(xpath);
			Matcher pathRootMatcher = Pattern.compile("([^/]+)").matcher(xpath); // greedy match of non-forward slashes
			String pathRoot =  pathRootMatcher.find() ? pathRootMatcher.group() : "";
			Pattern pattern = Pattern.compile("(<" + pathRoot + /* DON'T CLOSE THE ELEMENT with '>" - read through the xmlns prefix */ ".+?</" + pathRoot + ">)"); // reluctant match of entire element
	        
			for (int i=0; i < as.size(); ++i) {
		    	String annotationItem = as.get(i);
		    	// Don't even try to parse foreign annotations!
		    	// Especially since software like CellDesigner (celldesigner.org) put the xmlns prefix
		    	// in the root xml element and hence parsing the annotation alone will be problematic.
		    	// Note: there may be multiple "annotations" per "i"
		    	Matcher m = pattern.matcher(annotationItem);
		    	while (m.find()) {
		    		xpathSelector.setContextItem(buildNodeFromString(m.group()));
		    		Object o = (evaluateSingle
		    					? xpathSelector.evaluateSingle()
		    					: xpathSelector.evaluate());
		    		if (o != null) // otherwise keep looking...
		    			return o;
		    	}
		    }
		}
        return null;
    }
    
    private static String getNodeAttribute(XdmNode node, String key, String defaultValue) {
    	if (node != null) {
    		String v = node.getAttributeValue(new QName(key));
    		if (v != null)
    			return v;
    	}
    	return defaultValue;
    }
    
    private static ImMath encodeMathFromDouble(Double value) {
        ImMath math = new ImMath();
        math.cMML = MathMLUtils.encodeCMMLReal(value);
        math.pMML = MathMLUtils.encodePMMLReal(value);
        math.text = value.toString();
        return math;
    }
    
    private static ImMath encodeMathFromCMML(String cMML) throws IOException {
        ImMath math = new ImMath();
        math.cMML = cMML;
        math.pMML = MathMLUtils.convertContentToPresentationML(math.cMML);
        math.text = MathMLUtils.convertContentMLToMathText(math.cMML);
        return math;
    }
    
    private static ImMath encodeMathFromIdentifier(String id) {
        ImMath math = new ImMath();
        math.cMML = MathMLUtils.encodeCMMLIdentifier(id);
        math.pMML = MathMLUtils.encodePMMLIdentifier(id);
        math.text = id;
        return math;
    }
    
    private static ImMath getMathFromXPath(XdmItem contextNode, String key, ImMath defaultValue) throws SaxonApiException, IOException {
    	String cMML = serializeXPathToString(contextNode, key + "/math:math");
        return (cMML == null || cMML.isEmpty()
                ? defaultValue
                : encodeMathFromCMML(cMML));
    }
    private static ImMath getMathFromNode(XdmNode node, ImMath defaultValue) throws SaxonApiException, IOException {
    	XdmNode mathNode = (node != null && node.size() >= 1 && !node.itemAt(0).isAtomicValue()) ? (XdmNode)node.itemAt(0) : null;  
        String cMML = serializeNodeToString(mathNode);
        return (cMML == null || cMML.isEmpty()
                ? defaultValue
                : encodeMathFromCMML(cMML));
    }
    
    
    public static XPathSelector compileInchmanXPathSelector(String exp) throws SaxonApiException {
    	XPathCompiler comp = globalProcessor.newXPathCompiler();
    	comp.declareNamespace(DocPatcher.InchmanPrefix, DocPatcher.InchmanNS);
    	comp.declareNamespace(MathMLUtils.MathPrefix, MathMLUtils.MathNS);
    	return comp.compile(exp).load();
    }
    
    private static String getStringFromXPath(XdmItem xpathNode, String key, String defaultValue) throws SaxonApiException {
        XdmItem item = evalXPathSingleItem(xpathNode, DocPatcher.InchmanPrefix + ":" + key);
        if (item != null) {
        	String str = item.getStringValue();
        	if (!str.isEmpty())
                return str;
        }
        return defaultValue;
    }
    
    private static XdmNode evalXPathSingleNode(XdmItem contextItem, String subPath) throws SaxonApiException {
    	XdmItem item = evalXPathSingleItem(contextItem, subPath);
    	return (item != null && !item.isAtomicValue()
    			? (XdmNode)item
    			: null);
    }
    private static XdmItem evalXPathSingleItem(XdmItem contextItem, String subPath) throws SaxonApiException {
    	if (contextItem == null) return null;
    	XPathSelector s = compileInchmanXPathSelector(subPath);
    	s.setContextItem(contextItem);
    	return s.evaluateSingle();
    }
    private static XdmValue evalXPathList(XdmItem contextItem, String subPath) throws SaxonApiException {
    	if (contextItem == null) return null;
    	XPathSelector s = compileInchmanXPathSelector(subPath);
    	s.setContextItem(contextItem);
    	return s.evaluate();
    }
    
    public static XdmNode buildNodeFromString(String contents) throws SaxonApiException {
    	StreamSource src = new StreamSource(new StringReader(contents));
    	return globalProcessor.newDocumentBuilder().build(src);
    }
    
    private static String serializeXPathToString(XdmItem contextItem, String key) throws SaxonApiException {
        XdmNode n = evalXPathSingleNode(contextItem, DocPatcher.InchmanPrefix + ":" + key);
        return (n != null
        		? serializeNodeToString((XdmNode)n)
        		: null);
    }
    public static String serializeNodeToString(XdmNode node) throws SaxonApiException {
        return node != null
        		? globalProcessor.newSerializer().serializeNodeToString(node)
        				// chop the xml header off - as it's going to be embedded into an XML document
        				// FUTURE: use a regex with a more flexible matching pattern
        				.replaceFirst("<\\?xml version=\"1.0\" encoding=\"UTF-8\"\\?>", "")
        		: null;
    }
}
