/**
 * Editor.js
 * 
 * Created by Aidan Lane, July 2010.
 */

var eventTimeStampsStore = null;
var eventTimeStampsGrid  = null;

// code mirror editors
var initCodeMirrorTextAreaId = 'initCodeMirrorTextArea';
var eventsCodeMirrorTextAreaId = 'eventsCodeMirrorTextArea';
var DriftDiffusivityCodeMirrorTextAreaId = 'DriftDiffusivityCodeMirrorTextArea';
var updateFieldsCodeMirrorTextAreaId = 'updateFieldsCodeMirrorTextArea';
var newIndividualsMethodCodeMirrorTextAreaId = 'newIndividualsMethodCodeMirrorTextArea';

var initCodeMirrorEditor = null;
var eventsCodeMirrorEditor = null;
var DriftDiffusivityCodeMirrorEditor = null;
var updateFieldsCodeMirrorEditor = null;
var newIndividualsMethodCodeMirrorEditor = null;

// local storages
var localStorageUpdateFieldsKey = 'updateFields.script';
var getUpdateFieldsInLocalStorage = function() { var v=localStorage[localStorageUpdateFieldsKey]; return v ? v : ""; /* string, not JSON */ };
var setUpdateFieldsInLocalStorage = function(updateFields) { return localStorage[localStorageUpdateFieldsKey] = updateFields; /* string, not JSON */ };

var localStorageInitScriptKey = 'init.script';
var getInitScriptInLocalStorage = function() { var v=localStorage[localStorageInitScriptKey]; return v ? v : ""; /* string, not JSON */ };
var setInitScriptInLocalStorage = function(initScript) { return localStorage[localStorageInitScriptKey] = initScript; /* string, not JSON */ };

var localStorageNewIndividualsMethodKey = 'newIndividualsMethod.script';
var getNewIndividualsMethodInLocalStorage = function () { var v=localStorage[localStorageNewIndividualsMethodKey]; return v ? v : "";}
var setNewIndividualsMethodInLocalStorage = function (newIndividualsMethod) {return localStorage[localStorageNewIndividualsMethodKey] = newIndividualsMethod;}

var localStorageEventsScriptKey = 'events.script';
var getEventsScriptInLocalStorage = function() { var v=localStorage[localStorageEventsScriptKey]; return v ? v : ""; /* string, not JSON */ };
var setEventsScriptInLocalStorage = function(eventsScript) { return localStorage[localStorageEventsScriptKey] = eventsScript; /* string, not JSON */ };

var localStorageDriftDiffusivityMethodKey = 'DriftDiffusivity.method';
var getDriftDiffusivityMethodInLocalStorage = function() { var v=localStorage[localStorageDriftDiffusivityMethodKey]; return v ? v : ""; /* string, not JSON */ };
var setDriftDiffusivityMethodInLocalStorage = function(DriftDiffusivityMethod) { return localStorage[localStorageDriftDiffusivityMethodKey] = DriftDiffusivityMethod; /* string, not JSON */ };

var localStorageDriftDiffusivityNonLinearKey = 'DriftDiffusivity.nonLinear';
var localStorageDriftDiffusivityComputeMomentsKey = 'DriftDiffusivity.computeMoments';

// setters
var safeSetCodeMirrorValue = function(id, cm, value) {
	Ext.getCmp(id).setValue(value);
	if (cm) {
		if (cm.getValue() != value) {
			cm.setValue(value);
		}
		// always refresh
		cm.refresh();
	}
};
var safeSetInitCodeMirrorValue = function(value) {
	safeSetCodeMirrorValue(initCodeMirrorTextAreaId, initCodeMirrorEditor, value);
};
var safeSetEventsCodeMirrorValue = function(value) {
	safeSetCodeMirrorValue(eventsCodeMirrorTextAreaId, eventsCodeMirrorEditor, value);
};
var safeSetNewIndividualsMethodCodeMirrorValue = function(value) {
	safeSetCodeMirrorValue(newIndividualsMethodCodeMirrorTextAreaId, newIndividualsMethodCodeMirrorEditor, value);
};
var safeSetDriftDiffusivityCodeMirrorValue = function(value) {
	safeSetCodeMirrorValue(DriftDiffusivityCodeMirrorTextAreaId, DriftDiffusivityCodeMirrorEditor, value);
};
var safeSetUpdateFieldsCodeMirrorValue = function(value) {
	safeSetCodeMirrorValue(updateFieldsCodeMirrorTextAreaId, updateFieldsCodeMirrorEditor, value);	
};

Ext.Loader.setConfig({
    enabled: true
});

Ext.Loader.setPath('Ext.ux', 'http://cdn.sencha.io/ext-4.0.7-gpl/examples/ux');

Ext.require('Ext.ux.CheckColumn');

Ext.onReady( function() {
    
	Ext.tip.QuickTipManager.init();
//    Ext.History.init();
        


    createFileInfoForm();
    createOptionsForm();
    
    createParametersStoreAndGrid();
    createCompartmentsStoreAndGrid();
    createSpeciesStoreAndGrid();
    createReactionsStoreAndGrid();
    
		
    Ext.define('EventTimeStamp', {
        extend: 'Ext.data.Model',
        fields: [
            {name: 'id',        type: 'int'},
            {name: 'timeStamp', type: Ext.data.Types.MATHOBJECT}
        ],
        proxy: {
            type: 'localstorage',
            id: 'eventTimeStamps'
        }
    });
    
    eventTimeStampsStore = Ext.create('Ext.data.Store', {
        model: 'EventTimeStamp',
        storeId: 'eventTimeStampsStore',
        batchUpdateMode: 'complete',
        autoLoad: true,
        autoSync: true
    });
    
    eventTimeStampsGrid = Ext.create('Ext.grid.Panel', {
        xtype: 'grid',
        width: 140,
        title: 'Time Stamps',
        border: 0,
        style: { 'border-right': '1px solid #84878C' },
        store: eventTimeStampsStore,
        plugins: [
             Ext.create('Ext.grid.plugin.CellEditing', {
                 clicksToEdit: 1
             })
        ],
        multiSelect: true,
        hideHeaders: true,
        columns: [{
            dataIndex: 'timeStamp',
            flex: 1,
            field: {
                xtype: 'mathmlgridfield'
            },
            renderer: mathJaxPresentationMathMLColumnRenderer
        }],
        tools: [{
            type: 'plus',
            handler: function(event, toolEl, panel) {
                eventTimeStampsStore.add({timeStamp: createMathObjectFromReal(0)}); // TODO: increment time?
                eventTimeStampsStore.sync();
            }
        },{
            type: 'minus',
            handler: function(event, toolEl, panel) {
                eventTimeStampsStore.remove(eventTimeStampsGrid.getSelectionModel().getSelection());
                eventTimeStampsStore.sync();
            }
        }]
    });
    
    // TODO: on change record to local storage!
    var convertTextAreaToCodeMirrorEditor = function(id, setEditorCallback, onValueChangeCallback, mode) {
	    var textarea = Ext.getCmp(id);
	    textarea.addListener('afterrender', function() {
	    	var codemirror = CodeMirror.fromTextArea(textarea.inputEl.dom, {
	    		height: textarea.getHeight() + 'px',
		        width: textarea.getWidth(),
		    	content: textarea.getValue(),
		    	mode: mode,
		        lineNumbers: true,
		        textWrapping: false,
		        onChange: function(cm) {
		        	onValueChangeCallback(cm.getValue());
		        }
	    	});
	    	setEditorCallback(codemirror);
	    },
	    textarea,
	    {single: true}
		);
    };
    
    var insertIntoDriftDiffusivityCodeMirrorEditor = function(text) {
    	//var somethingWasAlreadySelected = DriftDiffusivityCodeMirrorEditor.somethingSelected();
    	DriftDiffusivityCodeMirrorEditor.replaceSelection(text); // will insert at cursor position and select if nothing currently selected
    	//if (!somethingWasAlreadySelected)
    	//	DriftDiffusivityCodeMirrorEditor.setSelection(null, null);
    };

    var viewport = Ext.create('Ext.container.Viewport', {
        layout: 'border',
        renderTo: Ext.getBody(),
        defaults: {
            collapsible: false,
            split: false,
            margins: '0 0 0 0'
        },
        items: [{
            region: 'north',
            xtype: 'box',
            contentEl: 'nav_section'
        },{
            id: 'west-panel',
            region:'west',
            //collapsible: true,
            style: { 'border-right': '1px solid #84878C'/*, opacity: '0.0'*/ },
            border: false,
            width: 250,
            layout: {
                type: 'vbox',
                align: 'stretch',
                pack: 'start'
            },
            items: [fileInfoForm, optionsForm, parametersGrid]
        },{
            id: 'center-panel',
            region:'center',
            xtype: 'tabpanel',
            style: { border: 'none', background: '#bbbbbd'/*, opacity: '0.0'*/ },
            frame: false,
            frameHeader: false,
            bodyStyle: { border: 'none', 'border-top': '1px solid #BCC0CC' },
            plain: true,
            padding: '1 0 0 0',
            activeTab: 0,
            items: [{
            	// Model Tab
                id: 'model',
                title: 'Model',
                layout: {
                    type: 'vbox',
                    align: 'stretch',
                    pack: 'start'
                },
                items: [{
                	xtype: 'box',
                	html: '<h1>Model Definition</h1>'
                	 	+ 'Language: <i>SBML (with Inchman extensions)</i>',
                	cls: 'editor-info-box'
                },
                compartmentsGrid,
                speciesGrid,
                reactionsGrid
                ]
            },{
            	// Initialization Tab
            	title: 'Initialization',
            	layout: {
                    type: 'vbox',
                    align: 'stretch',
                    pack: 'start'
                },
            	items: [{
            		xtype: 'box',
            		html: '<h1>Initialization Script</h1>'
            		 	+ 'Language: <i>Python</i>',
            		cls: 'editor-info-box'
            	},{
                    flex: 1,
            		layout: 'fit',
                    border: false,
                	// future: only enable each of these when each operation can actually be made
                	tbar: [{ tooltip: 'Undo', iconCls: 'x-btn-undo', handler: function() {initCodeMirrorEditor.undo();} },
                	       { tooltip: 'Redo', iconCls: 'x-btn-redo', handler: function() {initCodeMirrorEditor.redo();} }],
            		items: [{
		            	// Initialization Tab / Init Code Mirror Panel
		                id: initCodeMirrorTextAreaId,
		                xtype: 'textarea',
		                value: getInitScriptInLocalStorage(),
		                layout: 'fit'
            		}]
            	}]
            },{
            	// Update fields tab
                id: 'updateField',
                title: 'Update Fields',
                layout: {
                    type: 'vbox',
                    align: 'stretch',
                    pack: 'start'
                },
            	items: [{
            		xtype: 'box',
            		html: '<h1>Update Fields Method</h1>'
            		 	+ 'Language: <i>Open-CL</i>',
            		cls: 'editor-info-box'
            	},{
                    flex: 1,
                    border: false,
	                layout: {
	                    type: 'hbox',
	                    align: 'stretch',
	                    pack: 'start'
	                },
	                items: [
	                eventTimeStampsGrid,
	                {
	                	xtype: 'panel',
	                	layout: 'fit',
	                	flex: 1,
	                    border: false,
	                	// future: only enable each of these when each operation can actually be made
	                	tbar: [{ toolTip: 'Undo', iconCls: 'x-btn-undo', handler: function() {updateFieldsCodeMirrorEditor.undo();} },
	                	       { toolTip: 'Redo', iconCls: 'x-btn-redo', handler: function() {updateFieldsCodeMirrorEditor.redo();} }],
	                	items: [{
		                	id: updateFieldsCodeMirrorTextAreaId,
		                	xtype: 'textarea',
		                    value: getUpdateFieldsInLocalStorage(),
		                    layout: 'fit',
		                    border: false
	                	}]
	                }]
            	}]
            },{
            	// New Individuals Tab
                id: 'newIndividualsMethod',
                title: 'New Individuals Method',
                layout: {
                    type: 'vbox',
                    align: 'stretch',
                    pack: 'start'
                },
            	items: [{
            		xtype: 'box',
            		html: '<h1>New Individuals Method</h1>'
            		 	+ 'Language: <i>Open-CL</i>',
            		cls: 'editor-info-box'
            	},{
                    flex: 1,
                    border: false,
	                layout: {
	                    type: 'hbox',
	                    align: 'stretch',
	                    pack: 'start'
	                },
	                items: [
	                eventTimeStampsGrid,
	                {
	                	xtype: 'panel',
	                	layout: 'fit',
	                	flex: 1,
	                    border: false,
	                	// future: only enable each of these when each operation can actually be made
	                	tbar: [{ toolTip: 'Undo', iconCls: 'x-btn-undo', handler: function() {newIndividualsMethodCodeMirrorEditor.undo();} },
	                	       { toolTip: 'Redo', iconCls: 'x-btn-redo', handler: function() {newIndividualsMethodCodeMirrorEditor.redo();} }],
	                	items: [{
		                	id: newIndividualsMethodCodeMirrorTextAreaId,
		                	xtype: 'textarea',
		                    value: getNewIndividualsMethodInLocalStorage(),
		                    layout: 'fit',
		                    border: false
	                	}]
	                }]
            	}]
            }
	    ,{
            	// Event
                id: 'events',
                title: 'Events',
                layout: {
                    type: 'vbox',
                    align: 'stretch',
                    pack: 'start'
                },
            	items: [{
            		xtype: 'box',
            		html: '<h1>Events Script</h1>'
            		 	+ 'Language: <i>Python</i>',
            		cls: 'editor-info-box'
            	},{
                    flex: 1,
                    border: false,
	                layout: {
	                    type: 'hbox',
	                    align: 'stretch',
	                    pack: 'start'
	                },
	                items: [
	                eventTimeStampsGrid,
	                {
	                	xtype: 'panel',
	                	layout: 'fit',
	                	flex: 1,
	                    border: false,
	                	// future: only enable each of these when each operation can actually be made
	                	tbar: [{ toolTip: 'Undo', iconCls: 'x-btn-undo', handler: function() {eventsCodeMirrorEditor.undo();} },
	                	       { toolTip: 'Redo', iconCls: 'x-btn-redo', handler: function() {eventsCodeMirrorEditor.redo();} }],
	                	items: [{
		                	id: eventsCodeMirrorTextAreaId,
		                	xtype: 'textarea',
		                    value: getEventsScriptInLocalStorage(),
		                    layout: 'fit',
		                    border: false
	                	}]
	                }]
            	}]
            },
	    {
            	// Compute Drift Diffusivity Method
            	id: 'driftDiffusivityTab',
            	title: 'Drift Diffusivity',
            	xtype: 'panel',
            	frame: false,
                layout: {
                    type: 'vbox',
                    align: 'stretch',
                    pack: 'start'
                },
            	items: [{
            		id: 'driftDiffusivityInfoBox',
            		cls: 'editor-info-box',
            		layout: 'fit',
            		border: false,
            		items: [{
                		xtype: 'box',
                		html: '<h1>Compute Drift Diffusivity Method</h1>'
                			+ 'Language: <i>OpenCL</i><br>'
                			+ 'Note: <i>This feature is only supported by the Inhomogeneous solver.</i>',
            		},{
            			xtype: 'box',
            			height: 10
            		},{
            			id: 'driftDiffusivityAdvancedBox',
            			xtype: 'fieldset',
            			title: 'Advanced',
            			collapsible: true,
            			collapsed: true,
            			items: [{
	            			xtype: 'checkboxgroup',
	            			fieldLabel: 'After each timestep',
	            	        layout: 'fit',
	            			columns: 1,
	            	        vertical: true,
	                		items: [{
								id: 'driftDiffusivityNonLinearCheckBox',
								name: 'nonLinear',
								boxLabel: 'Re-compute Drift Diffusivity (Non-Linear)',
								tooltip: 'Enable non-linear to recompute diffusitivity and drift after each timestep.',
								handler: function(checkbox, checked) {
									localStorage[localStorageDriftDiffusivityNonLinearKey] = checked; 
								}
							},{
								id: 'driftDiffusivityComputeMomentsCheckBox',
								name: 'computeMoments',
								boxLabel: 'Compute Moments',
								tooltip: 'Enable to compute momnents after each timestep (currently limited to X average).',
								handler: function(checkbox, checked) {
									localStorage[localStorageDriftDiffusivityComputeMomentsKey] = checked; 
								}
							}]
            			}]
            		}]
            	},{
            		id: 'driftDiffusivityTabContents',
                    flex: 1,
                    border: false,
                    layout: {
	                    type: 'hbox',
	                    align: 'stretch',
	                    pack: 'start'
	                },
	                items: [{
	                	flex: 1,
	                    layout: 'fit',
	                    border: false,
		            	tbar: [// future: only enable each of these when each operation can actually be made
		            	       { toolTip: 'Undo', iconCls: 'x-btn-undo', handler: function() {DriftDiffusivityCodeMirrorEditor.undo();} },
		            	       { toolTip: 'Redo', iconCls: 'x-btn-redo', handler: function() {DriftDiffusivityCodeMirrorEditor.redo();} },
		            	{
		            	    xtype: 'tbseparator'
		            	},{
		            		xtype: 'tbtext',
							text: 'Insert:'
						},{
							text: 'Global',
							menu: {
								xtype: 'menu',
								listeners: { 'click': function(menu, item) {
									if (item && !item.isDisabled() && !item.menu /* valid leaf item (has no sub-menu) */)
										insertIntoDriftDiffusivityCodeMirrorEditor(item.text);
								}},
								items: [{
									text: 'Grid Coordinates',
									menu: {
										xtype: 'menu',
										bubbleEvents: [ 'click' ],
										items: [{ text: 'GridX' },
										        { text: 'GridY' },
										        { text: 'GridModelWidth' },
										        { text: 'GridModelHeight' }],
									}
						        },
						        {
						        	text: 'Physical Coordinates',
						        	menu: {
										xtype: 'menu',
										bubbleEvents: [ 'click' ],
										items: [{ text: 'PhysicalX' },
										        { text: 'PhysicalY' },
										        { text: 'PhysicalModelWidth' },
										        { text: 'PhysicalModelHeight' },
										        { text: 'PhysicalCellWidth' },
										        { text: 'PhysicalCellHeight' }]
						        	}
						        },
						        { text: 'NumReactions' },
						        { text: 'NumSpecies' },
						        { text: 'PhysicalSimTime' },
								]								
							}
						},{
							text: 'Species',
							menu: {
								xtype: 'menu',
								listeners: {
									'beforeshow': function(menu) {
										// future: introduce caching if this proves too slow
										speciesStore.load();
										menu.removeAll();
										var species = ['(all)'];
										species = species.concat(speciesStore.collect('sbmlId'));
										Ext.each(species, function(name, index) {
											menu.add({
												text: name,
												menu: {
													xtype: 'menu',
													bubbleEvents: [ 'click' ],
													items: [{ text: (index==0 ? '' : (name+'.')) + 'DiffusivityX' },
													        { text: (index==0 ? '' : (name+'.')) + 'DiffusivityY' },
													        { text: (index==0 ? '' : (name+'.')) + 'DriftX' },
													        { text: (index==0 ? '' : (name+'.')) + 'DriftY' },
													        { xtype: 'menuseparator', padding: 0 },
													        {
													        	text: 'Advanced',
													        	menu: {
																	xtype: 'menu',
																	bubbleEvents: [ 'click' ],
																	items: [{ text: (index==0 ? '' : (name+'.')) + 'MeanX' /* future: mark as read-only */}]
													        	}
													        }]
										        }
											});
										});
									},
									'click': function(menu, item) {
										if (item && !item.isDisabled() && !item.menu /* valid leaf item (has no sub-menu) */)
											insertIntoDriftDiffusivityCodeMirrorEditor(item.text);
									}
								}
							}
						},{
							text: 'Parameter',
							menu: {
								xtype: 'menu',
								listeners: {
									'beforeshow': function(menu) {
										// future: introduce caching if this proves too slow
										parametersStore.load();
										menu.removeAll();
										var params = parametersStore.collect('sbmlId');
										if (params.length == 0) {
											this.add({ text: '(none)', disabled: true });
										} else {
											Ext.each(params, function(name) {
												menu.add({ text: name });
											}, menu);
										}
									},
									'click': function(menu, item) {
										if (item && !item.isDisabled() && !item.menu /* valid leaf item (has no sub-menu) */)
											insertIntoDriftDiffusivityCodeMirrorEditor(item.text);
									}			
								}
							}
						},{
							xtype: 'tbfill'
						},{
							text: 'Load Preset',
							menu: {
								xtype: 'menu',
								// TODO: move these somewhere smarter!
								items: [{
									text: 'Homogeneous',
									method:
										'// sets the diffusivity and drift to the constants given by the host parameters\n'
							          + 'DiffusivityX = diffX;\n'
							          + 'DiffusivityY = diffY;\n'
							          + '\n'
							          + 'DriftX = driftX;\n'
							          + 'DriftY = driftY;'
								},{
									text: 'Multiplicative Noise',
									method:
							            'Real shiftedX = PhysicalX + PhysicalModelWidth;\n'
							          + 'Real shiftedY = PhysicalY + PhysicalModelHeight;\n'
							          + '\n'
							          + 'DiffusivityX = (sigmaX * sigmaX) * (shiftedX * shiftedX) / 2.0;\n'
							          + 'DiffusivityY = (sigmaY * sigmaY) * (shiftedY * shiftedY) / 2.0;'
							          + '\n'
							          + 'DriftX = muX * shiftedX;\n'
							          + 'DriftY = muY * shiftedY;'
								},{
									text: 'Ornstein-Uhlenbeck',
									method:
							            '// set up an Ornstein-Uhlenbeck process\n'
							          + 'DiffusivityX = diffX;\n'
							          + 'DiffusivityY = diffY;\n'
							          + '\n'
							          + 'DriftX = -(gammaXX*PhysicalX + gammaXY*PhysicalY);\n'
							          + 'DriftY = -(gammaYX*PhysicalX + gammaYY*PhysicalY);'
								},{
									text: 'Non-Linear',
									method:
							            '// diffusivity is constant\n'
							          + 'DiffusivityX = diffX;\n'
							          + 'DiffusivityY = diffY;\n'
							          + '\n'
							          + 'Real shiftedX = PhysicalX + origin;\n'
							          + 'Real shiftedMeanX = MeanX + origin;\n'
							          + '\n'
							          + 'DriftX = (-omega + shiftedX + theta*shiftedMeanX);\n'
							          + 'DriftY = 0.0;'
								}],
								listeners:{
							         scope: this,
							         'click': function(menu, item) {
							        	 if (item)
							        		 safeSetDriftDiffusivityCodeMirrorValue(item.method);
							         }
							    }
							}
						}],
		            	items: [{
		            		// TODO: swap this to using the C-style codemirror syntax, rathar than python one
		            		id: DriftDiffusivityCodeMirrorTextAreaId,
		            		xtype: 'textarea',
		            		value: getDriftDiffusivityMethodInLocalStorage(),
		            		layout: 'fit',
		            		flex: 1,
		            		border: false
		            	}]
	                }],
            	}
		
		
		
		
		],
	            listeners: {
	            	'tabchange': function(tabPanel, tab){
	                    //Ext.History.add(tab.id);
	                    Ext.Function.defer(function() { // TODO: remove the need for this defer workaround- it's so that the value is set once the editor is visible (doomed to fail though if the machine is slow or busy)
	                    	safeSetInitCodeMirrorValue(getInitScriptInLocalStorage());
				safeSetUpdateFieldsCodeMirrorValue(getUpdateFieldsInLocalStorage());
	                    	safeSetEventsCodeMirrorValue(getEventsScriptInLocalStorage());
	                    	safeSetDriftDiffusivityCodeMirrorValue(getDriftDiffusivityMethodInLocalStorage());
				safeSetNewIndividualsMethodCodeMirrorValue(getNewIndividualsMethodInLocalStorage());
	                    }, 10);
	                }
	            }
            }]
        },{
            region: 'south',
            xtype: 'box',
            height: 24,
            cls: 'editor-footer',
            //style: { background: '#AAA', 'border-top': '1px solid #84878C' },
            html: 'Copyright &copy; 2011-2013 Monash University'// &mdash; Version ' + BuildInfo.versionNumber + '.' + BuildInfo.revisionNumber
        }]
    });
    
/*
        // Handle this change event in order to restore the UI to the appropriate history state
        Ext.History.on('change', function(token) {
            var tab = token ? token : 0;
            Ext.getCmp('center-panel').setActiveTab(tab);
        });
*/
    convertTextAreaToCodeMirrorEditor(initCodeMirrorTextAreaId,
    								  function(cm) { initCodeMirrorEditor=cm; },
    								  setInitScriptInLocalStorage,
    								  'python');
    convertTextAreaToCodeMirrorEditor(updateFieldsCodeMirrorTextAreaId,
    								  function(cm) { updateFieldsCodeMirrorEditor=cm; },
    								  setUpdateFieldsInLocalStorage,
    								  'clike');
    convertTextAreaToCodeMirrorEditor(newIndividualsMethodCodeMirrorTextAreaId,
				      function(cm) { newIndividualsMethodCodeMirrorEditor = cm;},
				      setNewIndividualsMethodInLocalStorage,
				      'clike');
    convertTextAreaToCodeMirrorEditor(eventsCodeMirrorTextAreaId,
    								  function(cm) { eventsCodeMirrorEditor=cm; },
    								  setEventsScriptInLocalStorage,
    								  'python');
    convertTextAreaToCodeMirrorEditor(DriftDiffusivityCodeMirrorTextAreaId,
    								  function(cm) { DriftDiffusivityCodeMirrorEditor=cm; },
    								  setDriftDiffusivityMethodInLocalStorage,
    								  'clike');
    
    // Logic
    var onSolverChanged = function() {
    	var isStochastic = (optionsForm.down('#solver').getValue().search('stochastic') != -1);
    	optionsForm.down('#runs').setDisabled(!isStochastic);
    	var isInhomogeneous = (optionsForm.down('#solver').getValue().search('inhomogeneous') != -1);
    	viewport.down('#driftDiffusivityAdvancedBox').setDisabled(!isInhomogeneous);
    	viewport.down('#driftDiffusivityTabContents').setDisabled(!isInhomogeneous);
    };
    optionsForm.down('#solver').on('change', onSolverChanged);
    
    // Note: the Grids are configured to auto-load
    optionsForm.getForm().load();
    onSolverChanged(); // for init
    Ext.getCmp('driftDiffusivityNonLinearCheckBox').setValue(localStorage[localStorageDriftDiffusivityNonLinearKey]);
    Ext.getCmp('driftDiffusivityComputeMomentsCheckBox').setValue(localStorage[localStorageDriftDiffusivityComputeMomentsKey]);
    
    viewport.doLayout();
    
    //document.getElementById('loading-message').innerHTML = 'Loading (100%) Completed.';
    Ext.get('loading').remove();

    //setTimeout(function(){
        //Ext.Fx.syncFx();
        //Ext.get('west-panel').fadeIn();
        //Ext.get('center-panel').fadeIn();
    //}, 250);
    
}); // Ext.onReady
