/**
 * Species.js
 * 
 * Created by Aidan Lane, Wed Dec 8 2010
 * 
 * Note: ExtJS model/record IDs must be integers, but we also need alphanumeric SBML IDS.
 *       Thus, each contains an 'id' and an 'sbmlId'. The 'id' is only needed on the client.
 *       
 * FUTURE: factorize code iwth Recations.js
 */

var speciesCompartmentInitialAmountStore = null;
var speciesStore = null;
var speciesGrid  = null;

var createSpeciesStoreAndGrid = function() {
    
    // TODO: render MathML instead of text
    initialAmountsColumnRenderer = function(data/*, metadata, record, rowIndex, colIndex*/) {
        var text = '';
        Ext.each(data, function(item, index) {
            if (item
                    && item.compartmentSbmlId
                    && item.initialAmount)
            {
                compartmentsStore.load();
                var compartment = compartmentsStore.findRecord('sbmlId', item.compartmentSbmlId);
                text +=
                    (compartment ? compartment.get('name') : '?')
                    + '='
                    + (item.initialAmount ? item.initialAmount.text : '?')
                    + (index < data.length-1 ? ', ' : '');
            }
        });
        if (text == '') text = '(none)'; // no items added? -> (none)
        return text;
    };
    
    Ext.define('Species', {
        extend: 'Ext.data.Model',
        fields: [
             {name: 'id',                        type: 'int'},
             {name: 'sbmlId',                    type: 'string'},
             {name: 'name',                      type: 'string'},
             {name: 'diffusionConstant',         type: Ext.data.Types.MATHOBJECT},
             {name: 'individual',                type: 'bool', defaultValue: true},
             {name: 'compartmentInitialAmounts', type: 'auto', defaultValue: [] /* NOT empty string */}   /* allow javascript ARRAY object */
        ],
        proxy: {
            type: 'localstorage',
            id: 'species'
        }
    });
    speciesStore = Ext.create('Ext.data.Store', {
        model: 'Species',
        storeId: 'speciesStore',
        batchUpdateMode: 'complete',
        autoLoad: true,
        autoSync: true
    });
    
    
    Ext.define('SpeciesCompartmentInitialAmount', {
        extend: 'Ext.data.Model',
        fields: [
             {name: 'id',                type: 'int'},
             {name: 'compartmentSbmlId', type: 'string'},
             {name: 'initialAmount',     type: Ext.data.Types.MATHOBJECT}
        ],
        proxy: {
            type: 'localstorage', // use localstorage so that we can use autoLoad and autoSync, which simplify our work
            id: 'speciesCompartmentInitialAmount'
        }
    });
    speciesCompartmentInitialAmountStore = Ext.create('Ext.data.Store', {
        model: 'SpeciesCompartmentInitialAmount',
        batchUpdateMode: 'complete',
        autoLoad: true,
        autoSync: true
    });


    Ext.define('SpeciesCompartmentListField', {
        extend: 'Ext.form.field.Picker',
        alias: 'widget.speciescompartmentlistfield',
        
        editable: false, // prevent the user from typing text directly into the field
        matchFieldWidth: false,
        
        // private
        onFocus: function() {
            this.callParent();
            this.expand();
        },
        
        createPicker: function() {
            
            var initialCountsGrid = Ext.create('Ext.grid.Panel', {
                title: 'Initial Counts',
                hidden: true,
                floating: true,
                width: 220,
                height: 150,
                store: speciesCompartmentInitialAmountStore,
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
                        itemId: 'speciesCompartmentComboBox',
                        store: compartmentsStore,
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
                    header: 'Initial Count', // Amount
                    flex: 1,
                    dataIndex: 'initialAmount',
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
                            speciesCompartmentInitialAmountStore.add({
                                compartmentSbmlId: firstCompartment.get('sbmlId'),
                                initialAmount:     createMathObjectFromReal(0)
                            });
                            speciesCompartmentInitialAmountStore.sync();
                        }
                    }
                },{
                    type: 'minus',
                    handler: function(event, toolEl, panel, tc) {
                        // Note: it's safe to use initialCountsGrid here, as it works by reference
                        //       and it'll be set by the time this code is executed.
                        var selection = initialCountsGrid.getSelectionModel().getSelection();
                        speciesCompartmentInitialAmountStore.remove(selection);
                        speciesCompartmentInitialAmountStore.sync();
                    }
                }]
            });
            
            return initialCountsGrid;
        },
        
        onDataChanged: function() {
            var value = Ext.clone(getCleanRecordsArray(speciesCompartmentInitialAmountStore));
            this.setRawValue(value); // not setValue - don't reload the store
        },
        
        onExpand: function() {
            // VERY IMPORTANT:
            // We must suspend events while expanded so that this picker won't call cancelEdit()
            // (which would prevent data from being saved) if a sub-picker is focused.
            // And yes, we also need to queueSuspended so that editing finishes properly. 
            this.suspendEvents(true);
            
            speciesCompartmentInitialAmountStore.on('datachanged', this.onDataChanged, this);
        },
        
        collapse: function() {
            // Only collapse if the mathMLEditor and compartment combobox aren't visible
            // FUTURE: write a safer implementation
            if ((!window.mathMLEditor || !window.mathMLEditor.isVisible(true))
                && Ext.ComponentQuery.query('#speciesCompartmentComboBox[isExpanded]').length == 0)
            {
                speciesCompartmentInitialAmountStore.un('datachanged', this.onDataChanged, this);
                
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
            speciesCompartmentInitialAmountStore.each(function(model) {model.destroy();}); // TODO: remove the need for this workaround for removeAll not having a lasting effect!
            speciesCompartmentInitialAmountStore.removeAll();
            if (value) {
                //Ext.each(value, function(item) {
                //    speciesCompartmentInitialAmountStore.add(item);
                //});
                // clone -> don't allow it to mess with our value (so renderer can still work)
                speciesCompartmentInitialAmountStore.add(Ext.clone(value));
            }
            speciesCompartmentInitialAmountStore.sync();
            
            this.setRawValue(value);
            if (this.rendered) {
                this.validate();
            }
            
            return this;
        },
        
        setRawValue: function(value) {
            this.value = value;
            if (this.rendered) {
                this.inputEl.dom.value = this.value ? initialAmountsColumnRenderer(this.value) : '';
            }
            return this.value;
        }
    });

    speciesGrid = new Ext.grid.Panel({
        title: 'Species',
        flex: 1,
        store: speciesStore,
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
            header: 'Diffusion Constant',
            dataIndex: 'diffusionConstant',
            width: 120,
            field: {
                xtype: 'mathmlgridfield'
            },
            renderer: mathJaxPresentationMathMLColumnRenderer 
        },{
         xtype: 'checkcolumn', text : 'Individual', dataIndex : 'individual' 
        },{
            header: 'Compartment Initial Counts',
            dataIndex: 'compartmentInitialAmounts',
            flex: 1,
            field: {
                xtype: 'speciescompartmentlistfield'
            },
            renderer: initialAmountsColumnRenderer
        }],
        tools: [{
            type: 'plus',
            handler: function(event, toolEl, panel, tc) {
                var sbmlId = findUniqueSbmlId(speciesStore, 'Species');
                speciesStore.add({
                    sbmlId: sbmlId,
                    name:   sbmlId,
                    diffusionConstant: createMathObjectFromReal(0),
                    compartmentInitialAmounts: []
                });
                speciesStore.sync();
            }
        },{
            type: 'minus',
            handler: function(event, toolEl, panel, tc) {
                var selection = speciesGrid.getSelectionModel().getSelection();
                speciesStore.remove(selection);
                speciesStore.sync();
            }
        }]
    });    
}; // createSpeciesStoreAndGrid