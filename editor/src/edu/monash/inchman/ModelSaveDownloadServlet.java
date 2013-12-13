package edu.monash.inchman;
import java.io.IOException;
import java.util.logging.Level;

import javax.servlet.http.*;

import com.google.appengine.api.memcache.ErrorHandlers;
import com.google.appengine.api.memcache.MemcacheService;
import com.google.appengine.api.memcache.MemcacheServiceFactory;

@SuppressWarnings("serial")
public class ModelSaveDownloadServlet extends HttpServlet {
	
	public void doGet(HttpServletRequest req, HttpServletResponse resp) throws IOException {
		try {
			String id = req.getParameter("id");
			String filenameHint = req.getParameter("filenameHint"); // may not be present, in which case it'll be null
			
			// Try get it from the cache
			MemcacheService syncCache = MemcacheServiceFactory.getMemcacheService();
		    syncCache.setErrorHandler(ErrorHandlers.getConsistentLogAndContinue(Level.INFO));
		    byte[] compressedContents = (byte[]) syncCache.get(id);
		    if (compressedContents == null) {
		    	//throw new Exception("Cache entry not found for key: \"" + id + "\"");
		    	resp.setStatus(404); // Not Found
		    	return;
		    }
		    
		    // Uncompress the contents
		    String contents = Compressor.decompressString(compressedContents);
		    
		    // Write as UTF-8
		    byte[] contentsBytes = contents.getBytes("UTF-8");  
		    resp.setContentType("application/force-download");
		    if (filenameHint != null) resp.setHeader("Content-Disposition","attachment; filename=\"" + filenameHint + "\"");
		    resp.setCharacterEncoding("UTF-8");
		    resp.setContentLength(contentsBytes.length);  
		    resp.getOutputStream().write(contentsBytes);  
		    resp.getOutputStream().flush();
		}
		catch (Exception e) {
			e.printStackTrace();
			resp.setStatus(500);
		}
	}
}