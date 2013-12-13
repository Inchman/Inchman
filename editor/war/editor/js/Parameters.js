/**
 * Parameters.js
 * 
 * Created by Aidan Lane, Wed Dec 8 2010
 * 
 * Note: ExtJS model/record IDs must be integers, but we also need alphanumeric SBML IDS.
 *       Thus, each contains an 'id' and an 'sbmlId'. The 'id' is only needed on the client.
 *       
 * TODO: handle the the Name attribute properly
 * TODO: support strings (with enumerations)
 */

var parametersStore = null;
var parametersGrid  = null;

var createParametersStoreAndGrid = function() {

    paramDataColumnRenderer = function(data/*, metadata, record, rowIndex, colIndex*/) {
        var text = '';
        var pointsText = data.points + (data.points == 1 ? ' point' : ' points');
        switch (data.domain) {
            case 'single':
                text = data.from;
                break;
            case 'range':
                text = '[' + data.from + ', ' + data.to + '], '
                       + (data.step != ''
                           ? ('step ' + data.step)
                           : pointsText);
                break;
            case 'random':
                text = '{'
                        + ((data.from == data.to)
                            ? data.from 
                            : (data.from + ', ... ' + data.to))
                        + '}, ' + pointsText;
                break;
            case 'field':
                text = '[]';
        }
        return text;
    };
    
    Ext.define('Parameter', {
        extend: 'Ext.data.Model',
        fields: [
            {name: 'id',     type: 'int'},
            {name: 'sbmlId', type: 'string'},
            {name: 'data',   type: 'auto', defaultValue: {}}
        ],
        proxy: {
            type: 'localstorage',
            id: 'parameters'
        }
    });    
    
    parametersStore = Ext.create('Ext.data.Store', {
        model: 'Parameter',
        storeId: 'parametersStore',
        batchUpdateMode: 'complete',
        autoLoad: true,
        autoSync: true
    });
    
    parametersGrid = new Ext.grid.Panel({
        title: 'Parameters',
        frame: false,
        border: 0,
        /*bodyStyle: { background: 'rgb(218,223,229)' },*/
        flex: 1,
        store: parametersStore,
        plugins: [
            Ext.create('Ext.grid.plugin.CellEditing', {
               clicksToEdit: 1
            })
        ],
        multiSelect: true,
        hideHeaders: true,
        columns: [
          {
              dataIndex: 'sbmlId',
              width: 50,
              field: {
                  xtype: 'textfield',
                  selectOnFocus: true
              }
          },{
              id: 'data', // used by parametersGrid for autoExpandColumn
              dataIndex: 'data', // use by parametersStore
              field: {
                  xtype: 'paramfield'
              },
              renderer: paramDataColumnRenderer,
              flex: 1
          }
        ],
        tools: [
                {
                    type: 'plus',
                    handler: function(event, toolEl, panel, tc) {
                        var sbmlId = findUniqueSbmlId(parametersStore, 'p');
                        parametersStore.add({sbmlId: sbmlId, name: sbmlId, data: {type: 'float', domain: 'single', from: '0.0'}});
                        parametersStore.sync();
                    }
                },{
                    type: 'minus',
                    handler: function(event, toolEl, panel, tc) {
                        parametersStore.remove(parametersGrid.getSelectionModel().getSelection());
                        parametersStore.sync();
                    }
                }
        ]
    });
}; // createParametersStoreAndGrid


Ext.define('ParamField', {
    extend: 'Ext.form.field.Picker',
    alias: 'widget.paramfield',
    
    editable: false, // prevent the user from typing text directly into the field
    matchFieldWidth: false,
    
    // private
    onFocus: function() {
        this.callParent();
        this.expand();
    },
    
    createPicker: function() {
        var me = this;

        return Ext.create('ParamEditor', {
            hidden: true,
            width: 350,
            frame: true,
            floating: true,
            renderTo: document.body,
            listeners: {
                scope: me,
                change: me.onChange,
                finished: me.collapse
            },
            keyNavConfig: {
                esc: function() {
                    me.collapse();
                }
            }
        });
    },
    
    onExpand: function() {
        this.picker.setValue(this.value);
    },
    
    onChange: function() {
        this.setValue(this.picker.getValue());
    },
    
    getValue: function() {
        return this.value;
    },
    
    setValue: function(value) {
        this.setRawValue(value);
        if (this.rendered) {
            this.validate();
        }
        return this;
    },
    
    setRawValue: function(value) {
        this.value = value;
        if (this.rendered) {
            this.inputEl.dom.value = this.value ? paramDataColumnRenderer(this.value) : '';
        }
        return this.value;
    }
});


Ext.define('ParamEditor', {
    extend: 'Ext.form.Panel',
    alias: 'widget.ParamEditor',
    
    // private
    initComponent : function() {
        Ext.apply(this, {
            fieldDefaults: {
                labelWidth: 50,
                labelAlign: 'right'
            },
            items: [{
                        fieldLabel: 'Type',
                        xtype: 'fieldcontainer',
                        flex: 1,
                        layout: 'hbox',
                        items: [{
                            xtype: 'buttonsegment',
                            flex: 1,
                            defaults: { xtype: 'button', enableToggle: true, allowDepress: false, toggleGroup: 'type', toggleHandler: this.typeToggleHandler, width: 50 },
                            items: [{itemId: 'typeIntegerButton', id: 'param-type-integer', text: 'integer', cls: 'segment-left'},
                                    {itemId: 'typeFloatButton', id: 'param-type-float', text: 'float', cls: 'segment-right'}]
                        }]
                    },
                    {
                        fieldLabel: 'Domain',
                        xtype: 'fieldcontainer',
                        flex: 1,
                        layout: 'hbox',
                        items: [{
                            xtype: 'buttonsegment',
                            flex: 1,
                            defaults: { xtype: 'button', enableToggle: true, allowDepress: false, toggleGroup: 'domain', toggleHandler: this.domainToggleHandler, width: 50 },
                            items: [{itemId: 'domainSingleButton', id: 'param-domain-single', text: 'single', cls: 'segment-left'},
                                    {itemId: 'domainRangeButton', id: 'param-domain-range', text: 'range', cls: 'segment-middle'},
                                    {itemId: 'domainRandomButton', id: 'param-domain-random', text: 'random', cls: 'segment-middle'},
                                    {itemId: 'domainFieldButton', id: 'param-domain-field', text: 'field', cls: 'segment-right'}]
                        }]
                    },
                    {
                        fieldLabel: 'Value',
                        xtype: 'fieldcontainer',
                        layout: 'hbox',
                        items: [{
                                    xtype: 'textfield',
                                    itemId: 'fromEdit',
                                    width: 60,
                                    allowDecimals: true,
                                    selectOnFocus: true,
                                    listeners: {
                                        specialkey: function(field, e) {
                                            if (e.getKey() == e.ENTER) {
                                                var editor = field.up('.ParamEditor');
                                                editor.onFromEditBlurred(field); // ensure that it occurs straight away
                                                editor.fireEvent('finished', editor); // better than editor.hide(), as it will submit the value, rather than waiting for the field to loose focus
                                            }
                                        }
                                    }
                                },{
                                    xtype: 'label',
                                    itemId: 'toLabel',
                                    text:  'to',
                                    style: 'line-height: 22px; text-align: center',
                                    width: 20
                                },{
                                    xtype: 'textfield',
                                    itemId: 'toEdit',
                                    width: 60,
                                    allowDecimals: true,
                                    selectOnFocus: true,
                                    listeners: {
                                        specialkey: function(field, e) {
                                            if (e.getKey() == e.ENTER) {
                                                var editor = field.up('.ParamEditor');
                                                editor.onToEditBlurred(field); // ensure that it occurs straight away
                                                editor.fireEvent('finished', editor); // better than editor.hide(), as it will submit the value, rather than waiting for the field to loose focus
                                            }
                                        }
                                    }
                                },{
                                    xtype: 'buttonsegment',
                                    itemId: 'stepOrPointsSwitch',
                                    margins: "0 5",
                                    width: 80,
                                    defaults: { enableToggle: true, allowDepress: false, toggleGroup: 'stepOrPoints', width: 40 },
                                    items: [{xtype: 'button', itemId: 'stepButton', text: 'step', cls: 'segment-left'},
                                            {xtype: 'button', itemId: 'pointsButton', text: 'points', cls: 'segment-right'}]
                                },{
                                    xtype: 'textfield',
                                    itemId: 'stepOrPointsEdit',
                                    flex: 1,
                                    allowDecimals: true,
                                    selectOnFocus: true,
                                    listeners: {
                                        specialkey: function(field, e) {
                                            if (e.getKey() == e.ENTER) {
                                                var editor = field.up('.ParamEditor');
                                                editor.onStepOrPointsEditBlurred(field); // ensure that it occurs straight away
                                                editor.fireEvent('finished', editor); // better than editor.hide(), as it will submit the value, rather than waiting for the field to loose focus
                                            }
                                        }
                                    }
                                }]
                    }
                    ]
        }
        );
        
        if (Ext.isArray(this.initialConfig)){
            Ext.apply(this, {items:this.initialConfig});
        }
        
        ParamEditor.superclass.initComponent.call(this);
        
        this.down('#fromEdit').on('blur', this.onFromEditBlurred, this);
        this.down('#toEdit').on('blur', this.onToEditBlurred, this);
        this.down('#stepButton').on('toggle', this.onStepButtonToggled, this); // this.pointsButton.on('toggle', this.onPointsButtonToggle, this); -- only need to connect to one of the buttons
        this.down('#stepOrPointsEdit').on('blur', this.onStepOrPointsEditBlurred, this);
        
        this.addEvents('change');
        this.addEvents('finished');
    },
    
    show: function() {
        this.callParent(arguments);
        this.down('#fromEdit').focus(true, 100);
    },
    
    getValue: function() {
        return this.value;
    },
    
    setValue: function(value) {
        this.value = Ext.clone(value); // VERY important - can't just use a reference - we need to detect changes in order to save them
        
        if (!this.value)
            return;
        
        /*
         * Set from, to, step and points BEFORE type and domain,
         * so that the latter don't try to 'repair' / 'fix' them,
         * resetting their values.
         */
        
        this.down('#fromEdit').setValue(this.value.from);
        this.down('#toEdit').setValue(this.value.to);

        var step = this.value.step;
        var points = this.value.points;
        if (step != '') {
            this.down('#stepOrPointsEdit').setValue(step); // set value BEFORE toggle, so that the record is updated with the correct value, not ''
            this.down('#stepButton').toggle(true);
        }
        else if (points != '') {
            this.down('#stepOrPointsEdit').setValue(points); // set value BEFORE toggle, so that the record is updated with the correct value, not ''
            this.down('#pointsButton').toggle(true);
        }

        var typeButton = Ext.getCmp('param-type-' + this.value.type);
        if (typeButton == null) typeButton = this.down('#typeIntegerButton');
        if (typeButton.pressed)
            this.typeToggleHandler(typeButton);
        else
            typeButton.toggle(true);

        var domainButton = Ext.getCmp('param-domain-' + this.value.domain);
        if (domainButton == null) domainButton = this.down('#domainSingleButton');
        if (domainButton.pressed)
            this.domainToggleHandler(domainButton);
        else
            domainButton.toggle(true);
    },
    
    // private
    doCallback: function() {
        var me = this;
        me.fireEvent('change', me, me.getValue());
    },
    
    // private
    typeToggleHandler: function(button) {
        var p = button.up('.ParamEditor');
        var regex = (p.down('#typeIntegerButton').pressed
                     ? /^([-+]?\d+([Ee][-+]?[0-2]?\d{1,2})?)$/i                      // integer, with or without exponent
                     : /^([-+]?(\d+\.?\d*|\d*\.?\d+)([Ee][-+]?[0-2]?\d{1,2})?)$/i);  // floating point, with or without exponent
        p.down('#fromEdit').regex = regex;
        p.down('#toEdit').regex = regex;
        p.down('#stepOrPointsEdit').regex = regex;
        p.down('#fromEdit').validate();
        p.down('#toEdit').validate();
        p.down('#stepOrPointsEdit').validate();
        p.value.type = p.down('#typeIntegerButton').pressed ? 'integer' : 'float';
        p.doCallback();
    },
    
    // private
    domainToggleHandler: function(button) {
        var p = button.up('.ParamEditor');
        var isRange  = p.down('#domainRangeButton').pressed;
        var isRandom = p.down('#domainRandomButton').pressed;
        var isField = p.down('#domainFieldButton').pressed;
        
        // Fix values as needed
        if (isRange || isRandom)
        {
            var toEdit = p.down('#toEdit');
            if (toEdit.getValue() == '') {
                toEdit.setValue(p.down('#fromEdit').getValue());
                p.onToEditBlurred(toEdit);
            }
            
            var stepOrPointsEdit = p.down('#stepOrPointsEdit');
            if (stepOrPointsEdit.getValue() == '') {
                stepOrPointsEdit.setValue('1');
                p.onStepOrPointsEditBlurred(stepOrPointsEdit);
            }
            
            var stepButton = p.down('#stepButton');
            var pointsButton = p.down('#pointsButton');
            if (isRange
                && !(stepButton.pressed
                     || pointsButton.pressed)) {
                stepButton.toggle(true);
            }
            else if (isRandom) {
                pointsButton.toggle(true);
            }
        }
        
        // Show, hide and resize as needed
        if (!(isRange || isRandom)) {
            //if (p.isVisible()) p.getResizeEl().scale(220, p.getHeight()); else p.setWidth(220);
        }
        p.down('#toLabel').setVisible(isRange | isRandom);
        p.down('#fromEdit').setVisible(!isField);
        p.down('#toEdit').setVisible(isRange | isRandom);
        p.down('#stepOrPointsSwitch').setVisible(isRange | isRandom);
        p.down('#stepButton').setDisabled(!isRange);
        p.down('#stepOrPointsEdit').setVisible(isRange | isRandom);
        if (isRange | isRandom) {
            //if (p.isVisible()) p.getResizeEl().scale(350, p.getHeight()); else p.setWidth(350);
        }
        
        p.value.domain = isRandom ? 'random' : (isRange ? 'range' : (isField ? 'field' : 'single'));
        p.doCallback();
    },
    
    // private
    onFromEditBlurred: function(field) {
        this.value.from = field.getValue();
        this.doCallback();
    },
    // private
    onToEditBlurred: function(field) {
        this.value.to = field.getValue();
        this.doCallback();
    },
    // private
    onStepButtonToggled: function(button) {
        this.onStepOrPointsEditBlurred(this.down('#stepOrPointsEdit'));
    },
    // private
    onStepOrPointsEditBlurred: function(field) {
        this.value.step = this.down('#stepButton').pressed ? field.getValue() : '';
        this.value.points = this.down('#pointsButton').pressed ? field.getValue() : '';
        this.doCallback();
    }
});