package edu.monash.inchman;

import java.io.IOException;

import javax.servlet.http.*;

@SuppressWarnings("serial")
public class ConvertContentMLToMathTextAndPresentationMLServlet extends HttpServlet {

	public static class ResponseData extends CommonResponseData {
        public String mathText = "";
        public String presentationML = "";
    }
	
	public void doPost(HttpServletRequest req, HttpServletResponse resp) {
		ResponseData d = new ResponseData();
		try {
			d.mathText = MathMLUtils.convertContentMLToMathText(req.getParameter("contentML"));
			d.presentationML = MathMLUtils.convertContentToPresentationML(req.getParameter("contentML"));
	        d.success = String.valueOf(true);
		    d.message = "Success";
		} catch (IOException exception) {
			d.success = String.valueOf(false);
		    d.message = "Invalid Content MathML";
		}
		d.tryRespondAsJson(resp);
	}
}
