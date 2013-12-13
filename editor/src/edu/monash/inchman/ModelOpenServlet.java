package edu.monash.inchman;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;

import javax.servlet.http.*;

import org.apache.commons.fileupload.FileItemIterator;
import org.apache.commons.fileupload.FileUploadException;
import org.apache.commons.fileupload.servlet.ServletFileUpload;

@SuppressWarnings("serial")
public class ModelOpenServlet extends HttpServlet {
	
	public static class ResponseData extends CommonResponseData {
		public ImModel imModel = null;
	}
	
	public void doPost(HttpServletRequest req, HttpServletResponse resp) throws IOException {
		/*
		boolean isMultipart = ServletFileUpload.isMultipartContent(req);
		if (!isMultipart) {
			resp.getWriter().println("not isMultipart");
			// TODO: return error
			return;
		}
		*/

		// Create a new file upload handler
		ServletFileUpload upload = new ServletFileUpload();

		ResponseData d = new ResponseData();
		try {
			FileItemIterator iter = upload.getItemIterator(req);
			if (iter.hasNext()) { // not while, as we only support one file
				InputStream inputStream =  iter.next().openStream();
		        	
				// read the file
	            ByteArrayOutputStream buffer = new ByteArrayOutputStream();
			    byte[] temp = new byte[1024];
			    int read;
			    while ((read = inputStream.read(temp)) > 0){
			       buffer.write(temp, 0, read);
			    }
			    d.imModel = DocLoader.loadSbmlDocument(buffer.toString());
			    d.success = String.valueOf(true);
			    d.message = "Success";
			}
		}
		catch (FileUploadException e) {
			e.printStackTrace();
			d.success = String.valueOf(false);
			d.message = "Failed to upload file";
		}
		catch (IOException e) {
			e.printStackTrace();
			d.success = String.valueOf(false);
			// FUTURE: if caused by SAXException, print that message instead
			d.message = "Failed to open model: " + e.getLocalizedMessage();
		}
		catch (Exception e) {
			e.printStackTrace();
			d.success = String.valueOf(false);
			d.message = "Failed to open model."; // Generic Exception handler- don't give the message away - we don't know what me might otherwise disclose!
		}
		d.respondAsJson(resp);
	}
}
