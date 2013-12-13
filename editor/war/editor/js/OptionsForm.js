/**
 * OptionsForm.js
 * 
 * Created by Aidan Lane, Wed Dec 8 2010
 */

var optionsForm = null;
var optionsStore = null;
var simulateModelWindow = null;

var localStorageOptionsKey = 'options';
var getOptionsInLocalStorage = function() { return getJsonInLocalStorage(localStorageOptionsKey); /* JSON, not string */ };
var setOptionsInLocalStorage = function(options) { return setJsonInLocalStorage(localStorageOptionsKey, options); /* JSON, not string */ };


var createOptionsForm = function() {
	
	var localStorageSavesDisabled = false;
    
    saveOptions = function() {
    	if (!localStorageSavesDisabled)
    		setOptionsInLocalStorage(optionsForm.getForm().getValues());
    };
    
    var getStochasticOrDeterministicGroupValue = function() {
    	var value = optionsForm.down('#stochasticOrDeterministicGroup').getValue();
    	value = Ext.isArray(value) ? value[0] : value;
    	return value.stochasticOrDeterministicGroup;
    };
    var getInhomogeneousOrHomogeneousGroupValue = function() {
    	var value = optionsForm.down('#inhomogeneousOrHomogeneousGroup').getValue();
    	value = Ext.isArray(value) ? value[0] : value;
    	return value.inhomogeneousOrHomogeneousGroup;
    };
    
    updateSolverFieldFromRadioGroups = function() {
    	var solver = getStochasticOrDeterministicGroupValue() + '_' + getInhomogeneousOrHomogeneousGroupValue();
    	var solverField = optionsForm.down('#solver');
    	if (solverField.getValue() != solver) {
    		solverField.setValue(solver);
        	saveOptions();
    	}
    };
    
    updateSolverRadioGroupsFromField = function() {
    	var solver = optionsForm.down('#solver').getValue();
    	var stochasticOrDeterministic =
    		solver.search('deterministic') != -1
    		? 'deterministic'
    		: 'stochastic'; // fall-back to a sane default
    	var inhomogeneousOrHomogeneous =
    		solver.search('inhomogeneous') != -1 // don't search for just 'homogeneous', as it's common to both
    		? 'inhomogeneous'
    		: 'homogeneous'; // fall-back to a sane default
    	// batch so that doens't fire an event before inhomogeneousOrHomogeneousGroup has a chance to not be 'undefined'
    	optionsForm.down('#stochasticOrDeterministicGroup').batchChanges(function() {
	    	if (getStochasticOrDeterministicGroupValue() != stochasticOrDeterministic) {
	    		optionsForm.down('#stochasticOrDeterministicGroup').setValue({stochasticOrDeterministicGroup: stochasticOrDeterministic});
	    	}
	    	if (getInhomogeneousOrHomogeneousGroupValue() != inhomogeneousOrHomogeneous) {
	    		optionsForm.down('#inhomogeneousOrHomogeneousGroup').setValue({inhomogeneousOrHomogeneousGroup: inhomogeneousOrHomogeneous});
	    	}
    	});
    };
    
    updateSolverFeaturesFromRadioGroups = function() {
    	var isStochasticText = (getStochasticOrDeterministicGroupValue() == 'stochastic' ? 'YES' : 'NO');
    	optionsForm.down('#solverFeatureMultipleRuns').setValue(isStochasticText);
    	var isInhomogeneousText = (getInhomogeneousOrHomogeneousGroupValue() == 'inhomogeneous' ? 'YES' : 'NO');
    	optionsForm.down('#solverFeatureDrift').setValue(isInhomogeneousText);
    	optionsForm.down('#solverFeatureCustomMethod').setValue(isInhomogeneousText);
    	optionsForm.down('#solverFeatureAnisotropic').setValue(isInhomogeneousText);
    };

    optionsForm = Ext.create('Ext.form.Panel', {
        api: {
            load: function() {
                var options = getOptionsInLocalStorage();
                if (options == null) {
                    options = {
                        name:           'Untitled',
                        gridWidth:      createMathObjectFromReal(32),
                        gridHeight:     createMathObjectFromReal(32),
                        physicalWidth:  createMathObjectFromReal(40.0),
                        physicalHeight: createMathObjectFromReal(40.0),
                        runs:           createMathObjectFromReal(10),
                        time:           createMathObjectFromReal(100),
                        outputInterval: createMathObjectFromReal(100),
                        solver:         'stochastic_homogeneous'
                    };
                    setOptionsInLocalStorage(options);
                }
                {
                	// IMPORTANT:
                	// localStorageSavesDisabled is not just an optimization!
                	// Given that the entire form is re-serialized from the widgets
                	// on each change, it's possible that when changing one widget's
                	// value, another one yet to be changed is reset back to it's old
                	// value from its widget, instead of the new value in the localStorage.
                	// This mechanism prevents this problem from occurring.
                	localStorageSavesDisabled = true;
                	optionsForm.getForm().setValues(options);
                    localStorageSavesDisabled = false;
                }
                updateSolverRadioGroupsFromField();
                updateSolverFeaturesFromRadioGroups();
            }
        },
        title: 'Options',
        frame: false,
        border: 0,
        bodyCls: 'left-panel',
        layout: 'anchor',
        bodyPadding: 6,
        defaults: {
            anchor: '95%'
        },
        fieldDefaults: {
            labelWidth: 90,
            labelAlign: 'right'
        },
        items: [{
            xtype: 'displayfield',
            cls: 'options-section-header',
            fieldLabel: '<b>Model</b>',
            labelSeparator: ''
        },{
            xtype: 'textfield',
            fieldLabel: 'Name',
            name: 'name',
            itemId: 'name',
            listeners: {
                blur: function() { saveOptions(); }
            }
        },{
            xtype: 'mathmlformfield',
            fieldLabel: 'Grid Width',
            name: 'gridWidth',
            itemId: 'gridWidth',
            listeners: {
                change: function() {
                    optionsForm.down('#gridHeight').setValue(this.getValue());
                    saveOptions();
                }
            }
        },{
            xtype: 'mathmlformfield',
            fieldLabel: 'Grid Height',
            name: 'gridHeight',
            itemId: 'gridHeight',
            // FUTURE: re-enable this when the GPU code allows for non-squared regions
            disabled: true
        },{
            xtype: 'fieldcontainer',
            fieldLabel: 'Physical Width',
            layout: 'hbox',
            items: [{
                width: 110,
                xtype: 'mathmlformfield',
                name: 'physicalWidth',
                itemId: 'physicalWidth',
                listeners: {
                    change: function() {
                        optionsForm.down('#physicalHeight').setValue(this.getValue());
                        saveOptions();
                    }
                }
            },
            {
                xtype: 'label',
                html: '&nbsp;&micro;m',
                style: 'line-height: 22px; text-align: center'
            }]
        },{
            xtype: 'fieldcontainer',
            fieldLabel: 'Physical Height',
            layout: 'hbox',
            disabled: true,
            items: [{
                width: 110,
                xtype: 'mathmlformfield',
                name: 'physicalHeight',
                itemId: 'physicalHeight'
                // FUTURE: re-enable this when the GPU code allows for non-squared regions
            },
            {
                xtype: 'label',
                html: '&nbsp;&micro;m',
                style: 'line-height: 22px; text-align: center'
            }]
        },
        {
            xtype: 'displayfield',
            cls: 'options-section-header',
            fieldLabel: '<b>Simulation</b>',
            labelSeparator: ''
        },{
        	xtype: 'fieldcontainer',
        	fieldLabel: 'Solver',
        	layout: {
                type: 'vbox',
                align: 'stretch',
                pack: 'start'
            },
        	height: 90,
        	items: [{
	        	xtype: 'radiogroup',
	        	itemId: 'stochasticOrDeterministicGroup', // for lookup by the radiogroup change event functions
	        	columns: 1,
	        	vertical: true,
	        	width: 120,
	        	defaults: { name: 'stochasticOrDeterministicGroup', submitValue: false /* submit via solver hiddenfield */ },
	        	items: [
	        	        { boxLabel: 'Stochastic',    inputValue: 'stochastic' },
	        	        { boxLabel: 'Deterministic', inputValue: 'deterministic' }
	        	],
	        	listeners: {
	        		change: function() {
	        			updateSolverFieldFromRadioGroups();
	        			updateSolverFeaturesFromRadioGroups();
	        		}
	        	}
        	},{
        		xtype: 'image',
        		maxWidth: 100,
        		height: 2,
        		src: 'editor/images/horiz_divider.png'
        	},{
	        	xtype: 'radiogroup',
	        	itemId: 'inhomogeneousOrHomogeneousGroup', // for lookup by updateSolverField
	        	columns: 1,
	        	vertical: true,
	        	width: 120,
	        	defaults: { name: 'inhomogeneousOrHomogeneousGroup', submitValue: false /* submit via solver hiddenfield */ },
	        	items: [{ boxLabel: 'Homogeneous',   inputValue: 'homogeneous' },
	        	        { boxLabel: 'Inhomogeneous', inputValue: 'inhomogeneous' }
	        	],
	        	listeners: {
	        		change: function() {
	        			updateSolverFieldFromRadioGroups();
	        			updateSolverFeaturesFromRadioGroups();
	        		}
	        	}
	        }]
        },{
            xtype: 'mathmlformfield',
            fieldLabel: 'Runs',
            name: 'runs',
            itemId: 'runs',
            listeners: {
                change: function() { saveOptions(); }
            }
        },{
            xtype: 'fieldcontainer',
            fieldLabel: 'Max Time',
            layout: 'hbox',
            items: [{
                xtype: 'mathmlformfield',
                name: 'time',
                itemId: 'time',
                listeners: {
                    change: function() { saveOptions(); }
                },
                width: 110
            },{
                xtype: 'label',
                html: '&nbsp;s',
                style: 'line-height: 22px; text-align: center'
            }]
        },{
            xtype: 'fieldcontainer',
            fieldLabel: 'Output Interval',
            layout: 'hbox',
            items: [{
                xtype: 'mathmlformfield',
                name: 'outputInterval',
                itemId: 'outputInterval',
                listeners: {
                    change: function() { saveOptions(); }
                },
                width: 110
            },{
                xtype: 'label',
                html: '&nbsp;s',
                style: 'line-height: 22px; text-align: center'
            }]
        },{
        	xtype: 'hiddenfield',
        	name: 'solver',  // for form save/load
        	itemId: 'solver', // for lookup by updateSolverField
        	listeners: {
        		change: function() {
        			updateSolverRadioGroupsFromField();
        			updateSolverFeaturesFromRadioGroups();
        		}
        	}
        },{
        	xtype: 'fieldset',
        	title: 'Solver Features',
        	margin: '6 0 6 12',
        	collapsible: true,
        	collapsed: true,
        	defaults: { xtype: 'displayfield', labelWidth: 160 },
        	items: [{ itemId: 'solverFeatureMultipleRuns', fieldLabel: 'Multiple Runs', value: '-' },
        	        { itemId: 'solverFeatureDrift', fieldLabel: 'Drift', value: '-' },
        	        { itemId: 'solverFeatureCustomMethod', fieldLabel: 'Custom Drift Diffusivity Method', value: '-' },
        	        { itemId: 'solverFeatureAnisotropic', fieldLabel: 'Anisotropic Diffusion', value: '-' }
        	]
        }]
    });


}; // createOptionsForm