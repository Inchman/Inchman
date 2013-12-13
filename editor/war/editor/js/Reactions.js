/**
 * Reactions.js
 * 
 * Created by Aidan Lane, Wed Dec 8 2010
 * 
 * Note: ExtJS model/record IDs must be integers, but we also need alphanumeric SBML IDS.
 *       Thus, each contains an 'id' and an 'sbmlId'. The 'id' is only needed on the client.
 */

var reactionCompartmentIsAllowedsStore = null;
var reactionsStore = null;
var reactionsGrid  = null;

var createReactionsStoreAndGrid = function() {
    
	// future: render whole expression as MathML instead of text??
    speciesStoichiometryColumnRenderer = function(data/*, metadata, record, rowIndex, colIndex*/) {
        var text = '';
        Ext.each(data, function(item, index) {
            if (item
                && item.speciesSbmlId
                && item.stoichiometry)
            {
                text +=
                	(item.stoichiometry != 1.0 ? item.stoichiometry : '') // if ==1 then leave out the number, as it's implicit
                	+ item.speciesSbmlId
                    + (index < data.length-1 ? ' + ' : '');
            }
        });
        if (text == '') text = '(none)'; // no items added? -> (none)
        return text;
    };
	
    // todo: render value as MathML instead of text
    isAllowedsColumnRenderer = function(data/*, metadata, record, rowIndex, colIndex*/) {
        var text = '';
        Ext.each(data, function(item, index) {
            if (item
                && item.compartmentSbmlId
                && item.isAllowed)
            {
                //compartmentsStore.load();
                //var compartment = compartmentsStore.findRecord('sbmlId', item.compartmentSbmlId);
                text +=
                	item.compartmentSbmlId//(compartment ? compartment.get('name') : '?')
                    + '='
                    + (item.isAllowed ? item.isAllowed.text : '?')
                    + (index < data.length-1 ? ', ' : '');
            }
        });
        if (text == '') text = '(none)'; // no items added? -> (none)
        return text;
    };
    
    Ext.define('Reaction', {
        extend: 'Ext.data.Model',
        fields: [
             {name: 'id',                     type: 'int'},
             {name: 'sbmlId',                 type: 'string'},             
             {name: 'name',                   type: 'string'},
             {name: 'kineticLaw',             type: Ext.data.Types.MATHOBJECT},
             {name: 'reactants',              type: 'auto', defaultValue: [] /* NOT empty string */}, /* allow javascript ARRAY object  */
             {name: 'products',               type: 'auto', defaultValue: [] /* NOT empty string */}, /* allow javascript ARRAY object  */
             {name: 'compartmentIsAlloweds',  type: 'auto', defaultValue: [] /* NOT empty string */}  /* allow javascript ARRAY object  */
        ],
        proxy: {
            type: 'localstorage',
            id: 'reactions'
        }
    });
    reactionsStore = Ext.create('Ext.data.Store', {
        model: 'Reaction',
        storeId: 'reactionsStore',
        batchUpdateMode: 'complete',
        autoLoad: true,
        autoSync: true
    });
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    Ext.define('ReactionSpeciesStoichiometry', {
        extend: 'Ext.data.Model',
        fields: [
             {name: 'id',            type: 'int'},
             {name: 'speciesSbmlId', type: 'string'},
             {name: 'stoichiometry', type: 'float'}
        ],
        proxy: {
            type: 'localstorage', // use localstorage so that we can use autoLoad and autoSync, which simplify our work
            id: 'reactionSpeciesStoichiometry'
        }
    });
    reactionSpeciesStoichiometryStore = Ext.create('Ext.data.Store', {
        model: 'ReactionSpeciesStoichiometry',
        batchUpdateMode: 'complete',
        autoLoad: true,
        autoSync: true
    });
    
    Ext.define('ReactionSpeciesListField', {
        extend: 'Ext.form.field.Picker',
        alias: 'widget.reactionspecieslistfield',
        
        editable: false, // prevent the user from typing text directly into the field
        matchFieldWidth: false,
        
        // private
        onFocus: function() {
            this.callParent();
            this.expand();
        },
        
        createPicker: function() {
            
            var stoichiometryGrid = Ext.create('Ext.grid.Panel', {
                title: 'Species Stoichiometry',
                hidden: true,
                floating: true,
                width: 220,
                height: 150,
                store: reactionSpeciesStoichiometryStore,
                plugins: [
                    Ext.create('Ext.grid.plugin.CellEditing', {
                        clicksToEdit: 1
                    })
                ],
                columns: [{
                    header: 'Species',
                    dataIndex: 'speciesSbmlId',
                    field: {
                        xtype: 'combo',
                        itemId: 'reactionSpeciesComboBox',
                        store: Ext.getStore('speciesStore'),
                        valueField: 'sbmlId',
                        displayField: 'name',
                        queryMode: 'local',
                        editable: false // need to select it (because display=name, value=sbmlId)
                    },
                    renderer: function(value) {
                        var r = speciesStore.findRecord('sbmlId', value);
                        return r ? r.get('name') : '(Unknown)';
                    }
                },{
                    header: 'Stoichiometry',
                    flex: 1,
                    dataIndex: 'stoichiometry',
                    field: {
                        xtype: 'textfield'
                    }
                }],
                tools: [{
                    type: 'plus',
                    handler: function(event, toolEl, panel, tc) {
                        // TODO: disable the whole tool if no species are available (and show that it is too)
                        var firstSpecies = Ext.getStore('speciesStore').first();
                        if (firstSpecies) {
                            reactionSpeciesStoichiometryStore.add({
                                speciesSbmlId: firstSpecies.get('sbmlId'),
                                stoichiometry: 1.0
                            });
                            reactionSpeciesStoichiometryStore.sync();
                        }
                    }
                },{
                    type: 'minus',
                    handler: function(event, toolEl, panel, tc) {
                        // Note: it's safe to use stoichiometryGrid here, as it works by reference
                        //       and it'll be set by the time this code is executed.
                        var selection = stoichiometryGrid.getSelectionModel().getSelection();
                        reactionSpeciesStoichiometryStore.remove(selection);
                        reactionSpeciesStoichiometryStore.sync();
                    }
                }]
            });
            
            return stoichiometryGrid;
        },
        
        onDataChanged: function() {
            var value = Ext.clone(getCleanRecordsArray(reactionSpeciesStoichiometryStore));
            this.setRawValue(value); // not setValue - don't reload the store
        },
        
        onExpand: function() {
            // VERY IMPORTANT:
            // We must suspend events while expanded so that this picker won't call cancelEdit()
            // (which would prevent data from being saved) if a sub-picker is focused.
            // And yes, we also need to queueSuspended so that editing finishes properly. 
            this.suspendEvents(true);
            
            reactionSpeciesStoichiometryStore.on('datachanged', this.onDataChanged, this);
        },
        
        collapse: function() {
            // Only collapse if the species combobox isn't visible
            // FUTURE: write a safer implementation
            if (Ext.ComponentQuery.query('#reactionSpeciesComboBox[isExpanded]').length == 0)
            {
                reactionSpeciesStoichiometryStore.un('datachanged', this.onDataChanged, this);
                
                // VERY IMPORTANT:
                // See onExpand()
                this.resumeEvents();
                
                this.callParent();
            }
        },
        
        getValue: function() {
            return this.value;
        },
        
        setValue: function(value) {
            reactionSpeciesStoichiometryStore.each(function(model) {model.destroy();}); // TODO: remove the need for this workaround for removeAll not having a lasting effect!
            reactionSpeciesStoichiometryStore.removeAll();
            if (value) {
                //Ext.each(value, function(item) {
                //    reactionSpeciesStoichiometryStore.add(item);
                //});
                // clone -> don't allow it to mess with our value (so renderer can still work)
                reactionSpeciesStoichiometryStore.add(Ext.clone(value));
            }
            reactionSpeciesStoichiometryStore.sync();
            
            this.setRawValue(value);
            if (this.rendered) {
                this.validate();
            }
            
            return this;
        },
        
        setRawValue: function(value) {
            this.value = value;
            if (this.rendered) {
                this.inputEl.dom.value = this.value ? speciesStoichiometryColumnRenderer(this.value) : '';
            }
            return this.value;
        }
    });
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    

    Ext.define('ReactionCompartmentIsAllowed', {
        extend: 'Ext.data.Model',
        fields: [
             {name: 'id',                type: 'int'},
             {name: 'compartmentSbmlId', type: 'string'},
             {name: 'isAllowed',         type: Ext.data.Types.MATHOBJECT} // yes, we allow it to be programmable
        ],
        proxy: {
            type: 'localstorage', // use localstorage so that we can use autoLoad and autoSync, which simplify our work
            id: 'reactionCompartmentIsAllowed'
        }
    });
    reactionCompartmentIsAllowedsStore = Ext.create('Ext.data.Store', {
        model: 'ReactionCompartmentIsAllowed',
        batchUpdateMode: 'complete',
        autoLoad: true,
        autoSync: true
    });
    
    Ext.define('ReactionCompartmentListField', {
        extend: 'Ext.form.field.Picker',
        alias: 'widget.reactioncompartmentlistfield',
        
        editable: false, // prevent the user from typing text directly into the field
        matchFieldWidth: false,
        
        // private
        onFocus: function() {
            this.callParent();
            this.expand();
        },
        
        createPicker: function() {
            
            var isAllowedsGrid = Ext.create('Ext.grid.Panel', {
                title: 'Allowed In',
                hidden: true,
                floating: true,
                width: 220,
                height: 150,
                store: reactionCompartmentIsAllowedsStore,
                plugins: [
                    Ext.create('Ext.grid.plugin.CellEditing', {
                        clicksToEdit: 1
                    })
                ],
                columns: [{
                    header: 'Compartment',
                    dataIndex: 'compartmentSbmlId',
                    field: {
                        xtype: 'combo',
                        itemId: 'reactionCompartmentComboBox',
                        store: Ext.getStore('compartmentsStore'),
                        valueField: 'sbmlId',
                        displayField: 'name',
                        queryMode: 'local',
                        editable: false // need to select it (because display=name, value=sbmlId)
                    },
                    renderer: function(value) {
                        var r = compartmentsStore.findRecord('sbmlId', value);
                        return r ? r.get('name') : '(Unknown)';
                    }
                },{
                    header: 'Is Allowed',
                    flex: 1,
                    dataIndex: 'isAllowed',
                    field: {
                        xtype: 'mathmlgridfield'
                    },
                    renderer: mathJaxPresentationMathMLColumnRenderer
                }],
                tools: [{
                    type: 'plus',
                    handler: function(event, toolEl, panel, tc) {
                        // TODO: disable the whole tool if no compartments are available (and show that it is too)
                        var firstCompartment = Ext.getStore('compartmentsStore').first();
                        if (firstCompartment) {
                            reactionCompartmentIsAllowedsStore.add({
                                compartmentSbmlId: firstCompartment.get('sbmlId'),
                                isAllowed:         createMathObjectFromIdentifier('true')
                            });
                            reactionCompartmentIsAllowedsStore.sync();
                        }
                    }
                },{
                    type: 'minus',
                    handler: function(event, toolEl, panel, tc) {
                        // Note: it's safe to use isAllowedsGrid here, as it works by reference
                        //       and it'll be set by the time this code is executed.
                        var selection = isAllowedsGrid.getSelectionModel().getSelection();
                        reactionCompartmentIsAllowedsStore.remove(selection);
                        reactionCompartmentIsAllowedsStore.sync();
                    }
                }]
            });
            
            return isAllowedsGrid;
        },
        
        onDataChanged: function() {
            var value = Ext.clone(getCleanRecordsArray(reactionCompartmentIsAllowedsStore));
            this.setRawValue(value); // not setValue - don't reload the store
        },
        
        onExpand: function() {
            // VERY IMPORTANT:
            // We must suspend events while expanded so that this picker won't call cancelEdit()
            // (which would prevent data from being saved) if a sub-picker is focused.
            // And yes, we also need to queueSuspended so that editing finishes properly. 
            this.suspendEvents(true);
            
            reactionCompartmentIsAllowedsStore.on('datachanged', this.onDataChanged, this);
        },
        
        collapse: function() {
            // Only collapse if the mathMLEditor and compartment combobox aren't visible
            // FUTURE: write a safer implementation
            if ((!window.mathMLEditor || !window.mathMLEditor.isVisible(true))
                && Ext.ComponentQuery.query('#reactionCompartmentComboBox[isExpanded]').length == 0)
            {
                reactionCompartmentIsAllowedsStore.un('datachanged', this.onDataChanged, this);
                
                // VERY IMPORTANT:
                // See onExpand()
                this.resumeEvents();
                
                this.callParent();
            }
        },
        
        getValue: function() {
            return this.value;
        },
        
        setValue: function(value) {
            reactionCompartmentIsAllowedsStore.each(function(model) {model.destroy();}); // TODO: remove the need for this workaround for removeAll not having a lasting effect!
            reactionCompartmentIsAllowedsStore.removeAll();
            if (value) {
                //Ext.each(value, function(item) {
                //    reactionCompartmentIsAllowedsStore.add(item);
                //});
                // clone -> don't allow it to mess with our value (so renderer can still work)
                reactionCompartmentIsAllowedsStore.add(Ext.clone(value));
            }
            reactionCompartmentIsAllowedsStore.sync();
            
            this.setRawValue(value);
            if (this.rendered) {
                this.validate();
            }
            
            return this;
        },
        
        setRawValue: function(value) {
            this.value = value;
            if (this.rendered) {
                this.inputEl.dom.value = this.value ? isAllowedsColumnRenderer(this.value) : '';
            }
            return this.value;
        }
    });
    
    reactionsGrid = new Ext.grid.Panel( {
        
        title: "Reactions",
        flex: 2,
        store: reactionsStore,
        plugins: [
            Ext.create('Ext.grid.plugin.CellEditing', {
                clicksToEdit: 1
            })
        ],
        multiSelect: true,
        border: 0,
        columns: [{
            header: 'SBML ID',
            dataIndex: 'sbmlId',
            field: {
                xtype: 'textfield',
                selectOnFocus: true
            }
        },{
            header: 'Name',
            dataIndex: 'name',
            field: {
                xtype: 'textfield',
                selectOnFocus: true
            }
        },{
            // NOTE: this currently experimental!!! - GPGMP can't support arbitrary laws yet either
            header: 'Kinetic Law',
            dataIndex: 'kineticLaw',
            field: {
                xtype: 'mathmlgridfield'
            },
            renderer: mathJaxPresentationMathMLColumnRenderer
        },{
            header: 'Reactants',
            dataIndex: 'reactants',
            field: {
                xtype: 'reactionspecieslistfield'
            },
            renderer: speciesStoichiometryColumnRenderer
        },{
            header: 'Products',
            dataIndex: 'products',
            field: {
                xtype: 'reactionspecieslistfield'
            },
            renderer: speciesStoichiometryColumnRenderer
        },{
            header: 'Compartments Allowed In',
            dataIndex: 'compartmentIsAlloweds',
            flex: 1,
            field: {
                xtype: 'reactioncompartmentlistfield'
            },
            renderer: isAllowedsColumnRenderer
        }],
        tools: [{
            type: 'plus',
            handler: function(event, toolEl, panel, tc) {
                var sbmlId = findUniqueSbmlId(reactionsStore, 'Reaction');
                reactionsStore.add({
                    sbmlId:     sbmlId,
                    name:       sbmlId,
                    kineticLaw: createEmptyMathObject(),
                    reactants:  [],
                    products:   [],
                    compartmentIsAlloweds: []
                });
                reactionsStore.sync();
            }
        },{
            type: 'minus',
            handler: function(event, toolEl, panel, tc) {
                reactionsStore.remove(reactionsGrid.getSelectionModel().getSelection());
                reactionsStore.sync();
            }
        }]
    });

}; // createReactionsStoreAndGrid