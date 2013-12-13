/**
 * MathJaxUtils.js
 * 
 * Created by Aidan Lane, Wed Dec 8 2010
 */

Ext.data.Types.MATHOBJECT = {
    convert: function(v, data) {
        return v;
    },
    sortType: function(v) {
        return v.text;  // When sorting, order by text
    },
    type: 'MathObject'
};

var mathJaxCache = {};

var updateMathJax = function(el) {
    if (Ext.get(el)) { // make sure that it's ready
        MathJax.Hub.Queue(
            ["Typeset",MathJax.Hub,el],
            [function(el) {
                if (mathJaxCache[el]) {
                    mathJaxCache[el].mathJax = Ext.clone(Ext.get(el).dom.innerHTML);
                }
            }, el]
        );
    }
};

var mathJaxPresentationMathMLColumnRenderer = function(value, metadata, record, rowIndex, colIndex) {
    var divId = record.id + '-' + rowIndex + '-' + colIndex + '-mathjax'; // stable id, that will stay the same between updates 
    var html = '<div id="' + divId + '">' + value.pMML + '</div>';
    
    // TODO: more checking on this PROTOTYPE caching code... - ALSO, factorize it with MMLEditor
    if (mathJaxCache[divId]
        && mathJaxCache[divId].pMML == value.pMML
        && mathJaxCache[divId].mathJax) {
        html = mathJaxCache[divId].mathJax;
        //console.log('hit');
    } else {
        mathJaxCache[divId] = {pMML: value.pMML};
        Ext.defer(updateMathJax, 1, this, [divId]); // defer event - so that it's processed after the html result is returned
    }
    return html;
};

var createMathObjectFromReal = function(value) {
    return {
    	// Important: CellDesigner (celldesigner.org) doesn't like (and thus re-export) '<cn type="real">' (or any other type), it prefers it to be simply: '<cn>'
        cMML: '<math xmlns="http://www.w3.org/1998/Math/MathML"><cn>' + value + '</cn></math>',
        pMML: '<math xmlns="http://www.w3.org/1998/Math/MathML"><mn>' + value + '</mn></math>',
        text: '' + value // ensure that it's encoded as a string
    };
};

var createMathObjectFromIdentifier = function(id) {
	return {
        cMML: '<math xmlns="http://www.w3.org/1998/Math/MathML"><ci>' + id + '</ci></math>',
        pMML: '<math xmlns="http://www.w3.org/1998/Math/MathML"><mi>' + id + '</mi></math>',
        text: id
    };
};

var createEmptyMathObject = function(value) {
    return {
        cMML: '<math xmlns="http://www.w3.org/1998/Math/MathML"></math>',
        pMML: '<math xmlns="http://www.w3.org/1998/Math/MathML"></math>',
        text: ''
    };
};