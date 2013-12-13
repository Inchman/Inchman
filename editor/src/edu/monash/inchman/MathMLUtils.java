// Created by Aidan Lane, Thu Jun 23, 2011

package edu.monash.inchman;

import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.concurrent.ArrayBlockingQueue;

import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;

import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;
import org.jmathml.ASTNode;
import org.jmathml.ASTRootNode;
import org.jmathml.ASTToXMLElementVisitor;
import org.jmathml.FormulaFormatter;
import org.jmathml.MathMLReader;
import org.jmathml.TextToASTNodeMathParser2;


public class MathMLUtils {
	
    public static final String MathPrefix    = "math";
    public static final String MathNS        = "http://www.w3.org/1998/Math/MathML";
    
	
	public static String encodeCMMLReal(String value) {
		// Important: CellDesigner (celldesigner.org) doesn't like (and thus re-export) "<cn type='real'>" (or any other type), it prefers it to be simply: "<cn>"
        return "<math xmlns='http://www.w3.org/1998/Math/MathML'><cn>" + value + "</cn></math>";
    }
    public static String encodeCMMLReal(Double value) {
        return encodeCMMLReal(value.toString());
    }
    public static String encodePMMLReal(String value) {
        return "<math xmlns='http://www.w3.org/1998/Math/MathML'><mn>" + value + "</mn></math>";
    }
    public static String encodePMMLReal(Double value) {
        return encodePMMLReal(value.toString());
    }
    
    public static String encodeCMMLIdentifier(String id) {
        return "<math xmlns='http://www.w3.org/1998/Math/MathML'><ci>" + id + "</ci></math>";
    }
    public static String encodePMMLIdentifier(String id) {
        return "<math xmlns='http://www.w3.org/1998/Math/MathML'><mi>" + id + "</mi></math>";
    }
    
    private static ArrayBlockingQueue<Transformer> contentToPresentationMLTransformers = new ArrayBlockingQueue<Transformer>(10); // capacity of 10
    public static void transformContentToPresentationML(StreamSource xmlSource, StreamResult outputTarget)
    {
        Transformer t = null;
        
        try {
        	// "offer on demand"
        	// WARNING: don't try to pre-create all the transformers (e.g. all 10), as that will make the first call take a LONG time
        	if (contentToPresentationMLTransformers.peek() == null
        			&& contentToPresentationMLTransformers.remainingCapacity() > 0) {
        		try {
        			TransformerFactory tfactory = net.sf.saxon.TransformerFactoryImpl.newInstance(); // the standard one had issues with the large xslt file: TransformerFactory.newInstance();
                    Transformer newt = tfactory.newTransformer(new StreamSource(new StringReader(Mathmlc2p_xsl.contents)));
                    
                    contentToPresentationMLTransformers.offer(newt); // OFFER, NOT ADD, which makes this whole operation safe -> if we can't add it then that's ok, we'll just wait for one instead upon take()
                } catch (TransformerConfigurationException e) {
                    e.printStackTrace();
                }
        	}
        	
            t = contentToPresentationMLTransformers.take();    
        }
        catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Unexpected interruption");
        }
        
        try {
            t.transform(xmlSource, outputTarget);
        }
        catch (TransformerException e) {
            e.printStackTrace();
        }
        
        // Now put it back, even if the transform failed
        try {
            contentToPresentationMLTransformers.put(t);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Unexpected interruption");
        }
    }
	
	
    public static String convertMathTextToContentML(String mathText) throws IOException
    {
        // TODO: remove the need for this workaround- 
        
        TextToASTNodeMathParser2 parser = new TextToASTNodeMathParser2(); // Note: be sure to use JMathML from JMathML's svn head
        ASTRootNode root = new ASTRootNode();
        parser.parseString(mathText, root);

        ASTToXMLElementVisitor visitor = new ASTToXMLElementVisitor();
        root.accept(visitor);
        
        return getElementString(visitor.getElement())
        // fixes issue where JMathML interprets "^" and "pow" as simply "pow", not ImMathML's "power":
        .replaceAll(MathPrefix+":pow", MathPrefix+":power") // JMathML -> ImMathML
        .replaceAll(MathPrefix+":cn type=\"\\w+?\"", MathPrefix+":cn") // CellDesigner (celldesigner.org) doesn't like (and thus re-export) "<cn type='real'>" (or any other type), it prefers it to be simply: "<cn>"
        .replaceFirst("<\\?xml version=\"1.0\" encoding=\"UTF-8\"\\?>", "") // TODO: remove the need for this hack!
        .trim();
    }
    
	
	public static String convertContentToPresentationML(String contentML)
    {
        StringWriter writer = new StringWriter();

        transformContentToPresentationML(
                new StreamSource(new StringReader(contentML)),
                new StreamResult(writer));
        writer.flush();
        
        return writer.toString().replaceAll(MathPrefix+":", "");
    }

	
	// TODO: fix known bug - reaction CMML of k*A*b will be wrongly converted to '(k * B')
    //       cause: some reactions are encoded with more than 2 arguments per "apply" -> JMathML requires the use of nesting for that
    public static String convertContentMLToMathText(String contentML) throws IOException
    {
        // fixes issue where JMathML interprets "^" and "pow" as simply "pow", not ImMathML's "power":
        String fixedContentML = contentML.replaceAll(MathPrefix+":power", MathPrefix+":pow"); // ImMathML -> JMathML
        
        MathMLReader reader = new MathMLReader();
        ASTNode math = reader.parseMathMLFromString(fixedContentML);
        
        FormulaFormatter formatter = new FormulaFormatter();
        String text = formatter.formulaToString(math);
        
        // Try to simplify the math - e.g. "(3.0)" -> "3.0"
        // Note: can't just remove brackets simply when the text begins and ends with them: e.g. '(3*x)/(2*y)'
        // Check if the top level children have any children
        // TODO: support "(3.0 * x)" -> "3.0 * x"
        boolean grandChildExists = false;
        for (ASTNode n : math.getChildren()) {
            if (n.getNumChildren() > 0) {
                grandChildExists = true;
                break;
            }
        }
        if (!grandChildExists
            && text.startsWith("(")
            && text.endsWith(")")) {
            text = text.substring(1, text.length()-1).trim(); // yes, -1 (not -2) -- endIndex: the ending index, __exclusive__
        }
        
        return text;
    }

    
    private static String getElementString(org.jdom.Element node) {
        org.jdom.Document mathDoc = new org.jdom.Document();
        mathDoc.setRootElement(node);
        String xmlString = xmlToString(mathDoc);
        return xmlString;
    }
    private static String xmlToString(org.jdom.Document xmlDoc) {
        // pretty formatting/printing, as this is used (indirectly) by convertMathTextToContentML for user visible text
        XMLOutputter xmlOut = new XMLOutputter(Format.getPrettyFormat());
        return xmlOut.outputString(xmlDoc);
    }
}
