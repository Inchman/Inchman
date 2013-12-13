/**
 * MMLEditor.js
 * 
 * Created by Aidan Lane, Mon Jan 10 2011
 * 
 * FUTURE: handle focus better (e.g. tabbed-in) - automatically expand?
 *         or only do so if a value's been entered?
 */

Ext.define('MathMLBaseField', {
    extend: 'Ext.form.field.Picker',
    
    editable: false, // prevent the user from typing text directly into the field
    matchFieldWidth: false,
    
    createPicker: function() {
        // Share / re-use the heavy picker.
        // This also allows multiple fields to share the same state of text and CMML
        //  panel visibility, across both MathMLFormField and MathMLGridField.
        if (true) {// -- TODO: renable when it's not so quirky! e.g. so yo can click from one field to another explicitly collapsing the first and still no change event issues ... !window.mathMLEditor) {
            window.mathMLEditor = new MathMLEditor({
                hidden: true,
                width: 300,
                height: 200,
                floating: true
            });
        }
        else {
            window.mathMLEditor.clearListeners(); // just in case collapse() wasn't called
        }
        return window.mathMLEditor;
    },
    
    expand: function() {
        this.callParent();
        if (this.picker) {
            this.picker.clearListeners(); // just in case collapse() wasn't called
            this.picker.setValue(this.value); // set when no listeners attached
            this.picker.on('change', this.onChange, this);
            //this.picker.on('finished', this.collapse, this);
        }
    },
    
    collapse: function() {
        this.callParent();
        if (this.picker) {
            this.picker.clearListeners();
        }
    },
    
    onChange: function() {
        if (this.picker) // may not be set yet, as field may be focused but not expanded yet (e.g. tabbed to)
            this.setValue(this.picker.getValue());

        // TODO: remove the need for this workaround - checkChange() didn't always seem to work when needed...
        //this.checkChange();
        this.fireEvent('change', this.value, this.aidansOldValue);
        this.aidansOldValue = this.value;
    },
    
    // private
    onRender: function(container, position) {
        this.callParent(arguments);
        this.richOverlay = new Ext.Component({
            renderTo: this.bodyEl,
            style: {position:'absolute'},
            padding: 4
        });
        this.mon(this.richOverlay.el, {
            scope: this,
            click: function() { if (!this.isDisabled()) this.expand(); }
        });
    },

    // This is called by the form's getValues() method
    // Unlike the standard version, this will return the value even if it's disabled
    getSubmitData: function() {
        var data = {};
        data[this.name] = this.value;
        return data;
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
            var html = '';
            if (!Ext.isEmpty(value)) {
                var divId = this.id + '-mathjax'; // stable id, that will stay the same between updates 
                html = '<div id="' + divId + '">' + value.pMML + '</div>';
                
                // TODO: more checking on this PROTOTYPE caching code... - ALSO factorize it with mathJaxPresentationMathMLColumnRenderer
                if (mathJaxCache[divId]
                    && mathJaxCache[divId].pMML == value.pMML
                    && mathJaxCache[divId].mathJax) {
                    html = mathJaxCache[divId].mathJax;
                    //console.log('hit');
                } else {
                    mathJaxCache[divId] = {pMML: value.pMML};
                    Ext.defer(updateMathJax, 1, this, [divId]); // defer event - so that it's processed after the html result is returned
                }
                
                /*
                if (Ext.get(divId) == null
                    || Ext.get(divId).dom.value != html)
                    Ext.defer(updateMathJax, 1, this, [divId]); // defer event - so that it's processed after the html result is returned
                    */
            }
            if (this.richOverlay)
                this.richOverlay.el.dom.innerHTML = html;
        }
        return this.value;
    }
});


Ext.define('MathMLFormField', {
    extend: 'MathMLBaseField',
    alias: 'widget.mathmlformfield'
});


Ext.define('MathMLGridField', {
    extend: 'MathMLBaseField',
    alias: 'widget.mathmlgridfield',
    
    // private
    onFocus: function() {
        this.callParent();
        this.expand();
    }
});


Ext.define('MathMLEditor', {
    extend: 'Ext.panel.Panel',
    alias: 'widget.mathmleditor',
    
    // private
    textContents: null,
    // private
    cMMLContents: null,
    // private
    pMMLContents: null,
    
    // private
    processingEditorBlur: false,
    // private
    doCallbackWhenEditorBlurFinished: false,
    
    layout: {
    	type: 'vbox',
    	align: 'stretch'
    },
    
    // private
    initComponent : function() {
        Ext.apply(this, {
            defaults: {
                xtype: 'panel',
                //collapsible: true,
                titleCollapse: true,
                bodyBorder: false,
                layout: 'fit'
            },
            items: [{
                title: 'Equation',
                items: [{
                    itemId: 'mathTextEditor',
                    xtype: 'textfield',
                    fieldStyle: { border: '0px' },
                    listeners: {
                        specialkey: function(field, e) {
                            if (e.getKey() == e.ENTER) {
                                var editor = field.up('.mathmleditor');
                                editor.onMathTextEditorBlur(); // ensure that it occurs straight away
                                editor.hide(); // -- PROVED TO BE BAD: editor.fireEvent('finished', editor); // better than editor.hide(), as it will submit the value, rather than waiting for the field to loose focus
                            }
                        }
                    }
                }]
            },{
                title: 'Content MathML',
                //collapsed: true,
                flex: 1,
                items: [{
                    itemId: 'mathmlContentEditor',
                    xtype: 'textarea',
                    htmlEncode: true,
                    fieldStyle: { border: '0px' }
                }]
            }]
            }
        );
        
        if (Ext.isArray(this.initialConfig)){
            Ext.apply(this, {items:this.initialConfig});
        }
        
        MathMLEditor.superclass.initComponent.call(this);
        
        this.down('#mathTextEditor').on('blur', this.onMathTextEditorBlur, this);
        this.down('#mathmlContentEditor').on('blur', this.onMathmlContentEditorBlur, this);
        
        this.addEvents('change');
        //this.addEvents('finished');
    },
    
    show: function() {
        this.callParent(arguments);
        this.down('#mathTextEditor').focus(true, 100);
    },
    
    getValue: function() {
        return {cMML: this.cMMLContents, pMML: this.pMMLContents, text: this.textContents};
    },
    
    setValue: function(value) {
    	var local = Ext.clone(value); // VERY important - can't just use a reference - we need to detect changes in order to save them
        this.setContentEditorValue(local.cMML);
        this.pMMLContents = local.pMML;
        this.setMathTextEditorValue(local.text);
    },
    
    // private
    onMathTextEditorBlur: function() {
        var text = this.down('#mathTextEditor').getValue();
        this.processingEditorBlur = true;
        Ext.Ajax.request({
        	url: '/api/convertMathTextToContentAndPresentationML',
        	params: {mathText: text},
        	scope: this,
        	success: function(response) {
        		var data = Ext.decode(response.responseText);
        		if (data.success === 'true') {
	            	this.setContentEditorValue(data.contentML);
	                this.pMMLContents = data.presentationML;
	                this.textContents = text; // set this now, now that we know it's valid
	            }
	            else {
	                this.down('#mathTextEditor').markInvalid(data.message);
	            }
	            this.processingEditorBlur = false;
	            //if (this.doCallbackWhenEditorBlurFinished == true) {
	                this.doCallback();
	            //}
        	}
        });
    },
    
    // private
    onMathmlContentEditorBlur: function() {
        var cMML = this.down('#mathmlContentEditor').getValue();
        this.processingEditorBlur = true;
        Ext.Ajax.request({
        	url: '/api/convertContentMLToMathTextAndPresentationML',
        	params: {contentML: cMML},
        	scope: this,
        	success: function(response) {
        		var data = Ext.decode(response.responseText);
        		if (data.success === 'true') {
	                this.cMMLContents = cMML; // set this now, now that we know it's valid
	                this.pMMLContents = data.presentationML;
	                this.setMathTextEditorValue(data.mathText);
	            }
	            else {
	                this.down('#mathmlContentEditor').markInvalid(data.message);
	            }
	            this.processingEditorBlur = false;
	            //if (this.doCallbackWhenEditorBlurFinished == true) {
	                this.doCallback();
	            //}
        	}
        });
    },
    
    // private
    doCallback: function() {
        var me = this;
        me.fireEvent('change', me, me.getValue());
    },
    
    // private
    setContentEditorValue: function(value) {
        this.cMMLContents = value;
        this.down('#mathmlContentEditor').setValue(value);
    },
    
    // private
    setMathTextEditorValue: function(value) {
        this.textContents = value;
        this.down('#mathTextEditor').setValue(value);
    }
});
