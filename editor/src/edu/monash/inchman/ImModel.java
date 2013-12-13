// Created by Aidan Lane, Wed May 18, 2011

package edu.monash.inchman;

import java.util.ArrayList;
import java.util.List;

public class ImModel {

    public static class ImMath {
        public String cMML;
        public String pMML;
        public String text;
    };
    
    public static class ImOptions {
        public String name;
        public ImMath gridWidth;
        public ImMath gridHeight;
        public ImMath physicalWidth;
        public ImMath physicalHeight;
        public String solver;
        public ImMath runs;
        //public ImMath steps;
        public ImMath time;
        public ImMath outputInterval;
    };
    
    public static class ImBase {
        // no 'id' here -> force ExtJs to generate them as needed
        public String sbmlId;
        public String name;
    };
    
    public static class ImParameterData {
        public String type;
        public String domain;
        public String from;
        public String to;
        public String step;
        public String points;
    }
    public static class ImParameter extends ImBase {
        public ImParameterData data;
    };
    
    public static class ImCompartment extends ImBase {
        public ImMath x;
        public ImMath y;
        public ImMath width;
        public ImMath height;
    };
    
    public static class ImSpeciesCompartmentInitialAmount {
        public String compartmentSbmlId;
        public ImMath initialAmount;
    };
    public static class ImSpecies extends ImBase {
        public ImMath diffusionConstant;
        public boolean individual;
        public List<ImSpeciesCompartmentInitialAmount> compartmentInitialAmounts;
    };
    
    public static class ImReactionSpeciesStoichiometry {
        public String speciesSbmlId;
        public String stoichiometry;
    };
    
    public static class ImReactionCompartmentIsAllowed {
        public String compartmentSbmlId;
        public ImMath isAllowed; // programmable MathML value, like (virtually) everything else
    };
    public static class ImReaction extends ImBase {
        public ImMath kineticLaw;
        public List<ImReactionSpeciesStoichiometry> reactants;
        public List<ImReactionSpeciesStoichiometry> products;
        public List<ImReactionCompartmentIsAllowed> compartmentIsAlloweds;
    };
    
    public static class ImInit {
        public String type;
        public String script;
    };
    
    public static class ImEvents {
        public String type;
        public String script;
        public List<ImMath> timeStamps;
    };
    
    public static class ImDriftDiffusivity {
        public String method;
        public String nonLinear;
        public String computeMoments;
    };
    
    public static class ImUpdateFields {
    	public String method;
    };
    
    public static class ImNewIndividualsMethod {
    	public String method;
    };    
    
    // members
    public String                   origSbml;
    public ImOptions                options;
    public ArrayList<ImParameter>   parameters;
    public ArrayList<ImCompartment> compartments;
    public ArrayList<ImSpecies>     species;
    public ArrayList<ImReaction>    reactions;
    public ImInit                   init;
    public ImUpdateFields			updateFields;
    public ImNewIndividualsMethod	newIndividualsMethod;
    public ImEvents                 events;
    public ImDriftDiffusivity       driftDiffusivity;
    
    // goes through all reactions and checks if reactants/products have duplicates 
    // duplicates will be incorporated into the stoichiometry
    public void cleanReactions() {
    	// go through reaction list    	
    	for (ImReaction imReaction : reactions) {
    		// final reactant list
    		final ArrayList<ImReactionSpeciesStoichiometry> reactants = new ArrayList<ImReactionSpeciesStoichiometry>();
    		
    		// go through reactants
    		for (ImReactionSpeciesStoichiometry stoich : imReaction.reactants) {
    			// check if it's in the final list..
    			boolean found=false;
    			for (ImReactionSpeciesStoichiometry s: reactants) {
    				if (s.speciesSbmlId.equals(stoich.speciesSbmlId)){
    					found = true;
    					s.stoichiometry = String.valueOf(Double.parseDouble(s.stoichiometry) + Double.parseDouble(stoich.stoichiometry));
    				}
    			}
    			// if we haven't found it yet we need to add it to the list
    			if (!found)
    				reactants.add(stoich);
    		}
    		
    		// and set them
    		imReaction.reactants = reactants;
    		
    		// final product list
    		final ArrayList<ImReactionSpeciesStoichiometry> products = new ArrayList<ImReactionSpeciesStoichiometry>();
    		
    		// go through products
    		for (ImReactionSpeciesStoichiometry stoich : imReaction.products) {
    			// check if it's in the final list..
    			boolean found=false;
    			for (ImReactionSpeciesStoichiometry s: products) {
    				if (s.speciesSbmlId.equals(stoich.speciesSbmlId)){
    					found = true;
    					s.stoichiometry = String.valueOf(Double.parseDouble(s.stoichiometry) + Double.parseDouble(stoich.stoichiometry));
    				}
    			}
    			// if we haven't found it yet we need to add it to the list
    			if (!found)
    				products.add(stoich);
    		}
    		
    		// and set them
    		imReaction.products = products;
    		
    	}
    }
}
