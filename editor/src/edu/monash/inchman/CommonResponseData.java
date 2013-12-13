// Aidan Lane - 5 Oct, 2011

package edu.monash.inchman;

import java.io.IOException;

import javax.servlet.http.HttpServletResponse;

import com.google.gson.Gson;

public class CommonResponseData {
	
	public String success = String.valueOf(false);	
	public String message = "Unknown status";
    
    public void tryRespondAsJson(HttpServletResponse resp) {
    	try {
			respondAsJson(this, resp);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    
    public void respondAsJson(HttpServletResponse resp) throws IOException {
    	respondAsJson(this, resp);
    }
    
    public static void respondAsJson(Object responseData, HttpServletResponse resp) throws IOException {
		// From the Sencha ExtJS 4 documentation:
		// "The server response is parsed by the browser to create the document for the IFRAME.
		// If the server is using JSON to send the return object, then the Content-Type header
		// must be set to "text/html" in order to tell the browser to insert the text unchanged
		// into the document body."
		resp.setContentType("text/html"); // prevent wrapping the response up in <pre> tags
		resp.setCharacterEncoding("UTF-8"); // VERY important
        resp.setHeader("Cache-Control", "no-cache");
        
		new Gson().toJson(responseData, resp.getWriter());
		resp.getWriter().flush(); // ensure that it written completely
	}
}
