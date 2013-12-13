/**
 * FileInfoForm.js
 * 
 * Created by Aidan Lane, Wed Dec 8 2010
 */

var fileInfoForm = null;

var localStorageOrigSbmlKey = 'origSbml';
var getOrigSbmlInLocalStorage = function() { return localStorage[localStorageOrigSbmlKey]; /* string, not JSON */ };
var setOrigSbmlInLocalStorage = function(origSbml) { return localStorage[localStorageOrigSbmlKey] = origSbml; /* string, not JSON */ };

var guiClearAllData = null;
var clearAllData = null;

var createFileInfoForm = function() {
	
	guiClearAllData = function() {
		Ext.Msg.confirm('Clear All?',
                'Are you sure that you would like to clear all data? This cannot be undone.',
                function(buttonId) {
            if (buttonId == 'yes') {
                clearAllData();
                optionsForm.getForm().load(); // reload default data 
            }
         });
	};
    
    clearAllData = function() {
        /*
         * Clear all existing data
         */
        localStorage.removeItem(localStorageOrigSbmlKey);
        
        localStorage.removeItem(localStorageOptionsKey);
        optionsForm.getForm().load();
        
        parametersStore.each(function(model) {model.destroy();}); // TODO: remove the need for this workaround for removeAll not having a lasting effect!
        parametersStore.removeAll();
        parametersStore.sync();
        
        compartmentsStore.each(function(model) {model.destroy();}); // TODO: remove the need for this workaround for removeAll not having a lasting effect!
        compartmentsStore.removeAll();
	
        speciesStore.each(function(model) {model.destroy();}); // TODO: remove the need for this workaround for removeAll not having a lasting effect!
        speciesStore.removeAll();
        speciesStore.sync();
        
        reactionsStore.each(function(model) {model.destroy();}); // TODO: remove the need for this workaround for removeAll not having a lasting effect!
        reactionsStore.removeAll();
        reactionsStore.sync();
        
        safeSetInitCodeMirrorValue('');
        setInitScriptInLocalStorage(''); // still need a valid value, unlike: localStorage.removeItem(localStorageInitScriptKey);
        
        safeSetEventsCodeMirrorValue('');
        setEventsScriptInLocalStorage(''); // still need a valid value, unlike: localStorage.removeItem(localStorageEventsScriptKey);
        
        safeSetDriftDiffusivityCodeMirrorValue('');
        setDriftDiffusivityMethodInLocalStorage(''); // still need a valid value, unlike: localStorage.removeItem(localStorageDriftDiffusivityMethodKey);
        
	safeSetUpdateFieldsCodeMirrorValue('');
	setUpdateFieldsInLocalStorage('');
	
	safeSetNewIndividualsMethodCodeMirrorValue('');
	setNewIndividualsMethodInLocalStorage('');
	
        eventTimeStampsStore.each(function(model) {model.destroy();}); // TODO: remove the need for this workaround for removeAll not having a lasting effect!
        eventTimeStampsStore.removeAll();
        eventTimeStampsStore.sync();
        
        localStorage[localStorageDriftDiffusivityNonLinearKey] = true;
        localStorage[localStorageDriftDiffusivityComputeMomentsKey] = false;
        Ext.getCmp('driftDiffusivityNonLinearCheckBox').setValue(localStorage[localStorageDriftDiffusivityNonLinearKey]);
        Ext.getCmp('driftDiffusivityComputeMomentsCheckBox').setValue(localStorage[localStorageDriftDiffusivityComputeMomentsKey]);
    };
    
    
    var loadModelData = function(imModel)
    {
    	clearAllData();
        
        setOrigSbmlInLocalStorage(imModel.origSbml);
        
        setOptionsInLocalStorage(imModel.options);
        optionsForm.getForm().load();
        
        parametersStore.add(imModel.parameters);
        parametersStore.sync();
        
        compartmentsStore.add(imModel.compartments);
        compartmentsStore.sync();
        
        speciesStore.add(imModel.species);
        speciesStore.sync();
        
        reactionsStore.add(imModel.reactions);
        reactionsStore.sync();
       
        var initScript = imModel.init.script ? imModel.init.script : ''; // prevent it from being set to "undefined"
        safeSetInitCodeMirrorValue(initScript); // 'should' indirectly call setInitScriptInLocalStorage
        setInitScriptInLocalStorage(initScript);
        
	var newIndividualsMethod = imModel.newIndividualsMethod.method ? imModel.newIndividualsMethod.method : '';
	safeSetNewIndividualsMethodCodeMirrorValue(newIndividualsMethod);
	setNewIndividualsMethodInLocalStorage(newIndividualsMethod);
	
        var eventsScript = imModel.events.script ? imModel.events.script : ''; // prevent it from being set to "undefined"
        safeSetEventsCodeMirrorValue(eventsScript); // 'should' indirectly call setEventsScriptInLocalStorage
        setEventsScriptInLocalStorage(eventsScript);
        
        if (imModel.events.timeStamps) {
            Ext.Array.each(imModel.events.timeStamps, function(item) {
                eventTimeStampsStore.add({timeStamp: item});
            });
            eventTimeStampsStore.sync();
        }
        
        var driftDiffusivityMethod = imModel.driftDiffusivity.method ? imModel.driftDiffusivity.method : ''; // prevent it from being set to "undefined"
        safeSetDriftDiffusivityCodeMirrorValue(driftDiffusivityMethod); // 'should' indirectly call setDriftDiffusivityMethodInLocalStorage
        setDriftDiffusivityMethodInLocalStorage(driftDiffusivityMethod);
        
        localStorage[localStorageDriftDiffusivityNonLinearKey] = imModel.driftDiffusivity.nonLinear;
        localStorage[localStorageDriftDiffusivityComputeMomentsKey] = imModel.driftDiffusivity.computeMoments;
        Ext.getCmp('driftDiffusivityNonLinearCheckBox').setValue(localStorage[localStorageDriftDiffusivityNonLinearKey]);
        Ext.getCmp('driftDiffusivityComputeMomentsCheckBox').setValue(localStorage[localStorageDriftDiffusivityComputeMomentsKey]);


	var updateFieldsMethod = imModel.updateFields.method ? imModel.updateFields.method : '';
	safeSetUpdateFieldsCodeMirrorValue(updateFieldsMethod);
	setUpdateFieldsInLocalStorage(updateFieldsMethod);
    };
    
    
    var gatherModelDataForUpload = function()
    {
        var timeStamps = [];
        eventTimeStampsStore.each(function(model) {
            timeStamps.push(model.get('timeStamp'));
        });
	        
        return {
            origSbml:     getOrigSbmlInLocalStorage(),
            options:      getOptionsInLocalStorage(),
            parameters:   getCleanRecordsArray(parametersStore),
            compartments: getCleanRecordsArray(compartmentsStore),
            species:      getCleanRecordsArray(speciesStore),
            reactions:    getCleanRecordsArray(reactionsStore),
            init: {
                type:        'python',
                script:      getInitScriptInLocalStorage()
            },
	    updateFields: {
		method:		getUpdateFieldsInLocalStorage()
	    },
	    newIndividualsMethod: {
		method:		getNewIndividualsMethodInLocalStorage()
	    },
            events: {
                type:        'python',
                script:      getEventsScriptInLocalStorage(),
                timeStamps:  timeStamps
            },
            driftDiffusivity: {
            	method:         getDriftDiffusivityMethodInLocalStorage(),
            	nonLinear:      localStorage[localStorageDriftDiffusivityNonLinearKey],
            	computeMoments: localStorage[localStorageDriftDiffusivityComputeMomentsKey]
            }
        };
    };
    
    

    var modelOpen = function(form) {
        if (!form.isValid())
            return;
        
        form.submit({
            url: '/api/modelOpen',
            waitMsg: 'Loading model...',
            success: function(form, action) {
                if (action.result.success != 'true') // != 'true' better than == 'false', as it covers the case of null and the alike
                {
                    if (action.result.redirectUrl
                    		&& action.result.redirectUrl != '')
                        window.location.href = action.result.redirectUrl;
                    else
                        Ext.Msg.alert('Open Failed', action.result.message);
                }
                else
                {
                	loadModelData(action.result.imModel);
                }
            },
            // TODO: sure that this is actually called upon failure -- see the submitExperiment code to see that a bit of extra work may be required
            failure: function(form, action) {
                Ext.Msg.alert('Open Failed', action.result.message);
            }
        });
    };


    var modelSave = function() {
        var waitMsgBox = Ext.Msg.wait('Preparing model for save, please wait...');
        
        console.log('Data sent:'+Ext.JSON.encode(getCleanRecordsArray(reactionsStore)));

        Ext.Ajax.request({
            url: '/api/modelSavePrepare',
            jsonData: gatherModelDataForUpload(),
            callback: function(options, success, response) {
            	waitMsgBox.close();
                // TODO: handle case where responseText is not valid JSON, that is - e.g. a struts2 error message
                var actualResponse = Ext.JSON.decode(response.responseText);
                if (actualResponse.success === 'true') {
                	window.location.href = '/api/modelSaveDownload?id=' + actualResponse.id + '&filenameHint=' + actualResponse.filenameHint;
                }
                else {
                    if (!actualResponse.message)
                        actualResponse.message = 'An unknown error occured.';
                
                    Ext.Msg.alert('Save Failed', actualResponse.message);
                }
            }
        });
    };
    
    
    var submitExperiment = function() {
    	var waitMsgBox = Ext.Msg.wait('Simulating experiment, please wait...',
                null,
                { increment: 20 } // may take a while, so double the number of increments before wrapping
        );
    	Ext.Ajax.request({
            url: '/api/modelSavePrepare',
            jsonData: gatherModelDataForUpload(),
            callback: function(options, success, response) {
            	waitMsgBox.close();
                // TODO: handle case where responseText is not valid JSON, that is - e.g. a struts2 error message
                var actualResponse = Ext.JSON.decode(response.responseText);
                if (actualResponse.success === 'true') {
                	window.location.href = 'inchman://simulate/' + actualResponse.id;
                }
                else {
                    if (!actualResponse.message)
                        actualResponse.message = 'An unknown error occured.';
                
                    Ext.Msg.alert('Simulation Failed', actualResponse.message);
                }
            }
        });
    };

    
   fileInfoForm = Ext.create('Ext.form.Panel', {
        frame: false,
        bodyCls: 'left-panel',
        border: 0,
        bodyBorder: false,
        bodyPadding: 4,
        layout: 'hbox',
        height: 42,
        minButtonWidth: null,
        items: [{
            // Note: this field must live inside a form with 'standardSubmit: true'
            name: 'upload',
            xtype: 'filefield',
            buttonOnly: true,
            width: 38, // TODO: remove this hack!
            buttonConfig: {
                cls: 'inchman-tool-button',
                text: '', // replace the default with emptiness
                icon: '/editor/images/open.png',
                scale: 'large',
                tooltip: 'Open&hellip;' // TODO: ensure that this is displayed
            },
            listeners: {
                change: function() {
                	var form = this.up('form').getForm();
                	if (form.isValid())
                		modelOpen(form);
                }
            }
        },{
            xtype: 'button',
            cls: 'inchman-tool-button',
            icon: '/editor/images/new.png',
            scale: 'large',
            tooltip: 'New',
            handler: function(button, event) {
                guiClearAllData();
            }
        },{
            xtype: 'button',
            cls: 'inchman-tool-button',
            icon: '/editor/images/save.png',
            scale: 'large',
            tooltip: 'Save',
            handler: function(button, event) {
                modelSave();
            }
        },
        { xtype: 'component', flex: 1 },
        {
            xtype: 'button',
            cls: 'inchman-tool-button',
            icon: '/editor/images/simulate.png',
            scale: 'large',
            tooltip: 'Submit Experiment&hellip;',
            handler: function(button, event) {
                Ext.Msg.confirm('Submit Experiment',
                        'Are you sure that you would like to submit this experiment?',
                        function(buttonId) {
                    if (buttonId == 'yes') {
                        submitExperiment();
                    }
                });
            }
        }]
    });

}; // createFileInfoForm