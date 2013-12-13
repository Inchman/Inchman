/**
 * Compartments.js
 * 
 * Created by Aidan Lane, Wed Dec 8 2010
 *
 * Note: ExtJS model/record IDs must be integers, but we also need alphanumeric SBML IDS.
 *       Thus, each contains an 'id' and an 'sbmlId'. The 'id' is only needed on the client.
 */

var compartmentsStore = null;
var compartmentsGrid  = null;

var createCompartmentsStoreAndGrid = function() {
    
    var createMMLColumn = function(id) {
        return new Ext.grid.Column({
            header: id,
            dataIndex: id,
            width: 65,
            field: {
                xtype: 'mathmlgridfield'
            },
            renderer: mathJaxPresentationMathMLColumnRenderer
        });
    };
    
    Ext.define('Compartment', {
        extend: 'Ext.data.Model',
        fields: [
             {name: 'id',     type: 'int'},
             {name: 'sbmlId', type: 'string'},
             {name: 'name',   type: 'string'},
             {name: 'x',      type: Ext.data.Types.MATHOBJECT},
             {name: 'y',      type: Ext.data.Types.MATHOBJECT},
             {name: 'width',  type: Ext.data.Types.MATHOBJECT},
             {name: 'height', type: Ext.data.Types.MATHOBJECT}
        ],
        proxy: {
            type: 'localstorage',
            id: 'compartments'
        }
    });
    
    compartmentsStore = Ext.create('Ext.data.Store', {
        model: 'Compartment',
        storeId: 'compartmentsStore',
        batchUpdateMode: 'complete',
        autoLoad: true,
        autoSync: true
    });
    
    compartmentsGrid = Ext.create('Ext.grid.Panel', {
        title: 'Compartments',
        flex: 1,
        store: compartmentsStore,
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
            width: 80,
            field: {
                xtype: 'textfield',
                selectOnFocus: true
            }
        },{
            header: 'Name',
            dataIndex: 'name',
            width: 100,
            field: {
                xtype: 'textfield',
                selectOnFocus: true
            }
        },
        createMMLColumn('x'),
        createMMLColumn('y'),
        createMMLColumn('width'),
        createMMLColumn('height'),
        {
            // fill to edge
            flex: 1,
            sortable: false,
            hideable: false,
            menuDisabled: true,
            draggable: false
        }
        ],
        tools: [{
            type: 'plus',
            handler: function(event, toolEl, panel, tc) {
                var sbmlId = findUniqueSbmlId(compartmentsStore, 'Compart');
                compartmentsStore.add({
                    sbmlId: sbmlId,
                    name:   sbmlId,
                    x:      createMathObjectFromReal(0),
                    y:      createMathObjectFromReal(0),
                    width:  createMathObjectFromReal(1),
                    height: createMathObjectFromReal(1)
                });
                compartmentsStore.sync();
            }
        },{
            type: 'minus',
            handler: function(event, toolEl, panel, tc) {
                compartmentsStore.remove(compartmentsGrid.getSelectionModel().getSelection());
                compartmentsStore.sync();
            }
        }
        ]
    });    
}; // createCompartmentsStoreAndGrid