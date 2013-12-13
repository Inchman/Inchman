package edu.monash.inchman;

import java.io.IOException;

import javax.servlet.http.*;

@SuppressWarnings("serial")
public class ConvertMathTextToContentAndPresentationMLServlet extends HttpServlet {

	public static class ResponseData extends CommonResponseData {
        public String contentML = "";
        public String presentationML = "";
    }
	
	public void doPost(HttpServletRequest req, HttpServletResponse resp) throws IOException {
		ResponseData d = new ResponseData();
		try {
			d.contentML = MathMLUtils.convertMathTextToContentML(req.getParameter("mathText"));
			d.presentationML = MathMLUtils.convertContentToPresentationML(d.contentML);
	        d.success = String.valueOf(true);
		    d.message = "Success";
		} catch (IOException exception) {
			d.success = String.valueOf(false);
		    d.message = "Invalid or unsupported equation";
		}
        d.tryRespondAsJson(resp);
	}
}
