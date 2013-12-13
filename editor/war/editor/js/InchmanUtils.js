/**
 * InchmanUtils.js
 * 
 * Created by Aidan Lane, Mon Jan 10 2011
 */

var getCleanRecordsArray = function(store) {
    var array = [];
    store.each(function(record) {
        var clean = {};
        record.fields.each(function(field){
            clean[field.name] = record.get(field.name);
        });
        array.push(clean);
    });
    return array;
};

var getJsonInLocalStorage = function(key) { return Ext.JSON.decode(localStorage[key]); };
var setJsonInLocalStorage = function(key, value) { return localStorage[key] = Ext.JSON.encode(value); };

var findUniqueSbmlId = function(store, idBase) {
    var id = idBase;
    if (store.findRecord('sbmlId', id) != null) {
        var suffix = 2;
        while (store.findRecord('sbmlId', id + suffix) != null) {
            suffix++;
        }
        id += suffix;
    }
    return id;
};