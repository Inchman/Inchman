package edu.monash.inchman;
import java.io.IOException;
import java.util.UUID;
import java.util.logging.Level;

import javax.servlet.http.*;

import com.google.appengine.api.memcache.ErrorHandlers;
import com.google.appengine.api.memcache.Expiration;
import com.google.appengine.api.memcache.MemcacheService;
import com.google.appengine.api.memcache.MemcacheServiceFactory;
import com.google.gson.Gson;

@SuppressWarnings("serial")
public class ModelSavePrepareServlet extends HttpServlet {
	
	public static class ResponseData extends CommonResponseData {
		public String id = null;
		public String filenameHint = null;
	}
	
	public void doPost(HttpServletRequest req, HttpServletResponse resp) throws IOException {
		ResponseData d = new ResponseData();
		try {
			// Patch the SBML document with the latest changes
			ImModel model = new Gson().fromJson(req.getReader(), ImModel.class);
			
			// clean the reactions
			model.cleanReactions();
			
			String contents = DocPatcher.patchSbmlDocumentToString(model);
			
			// Compress the contents, as it's XML data that should compress up well
			byte[] compressedContents = Compressor.compressString(contents);
			
			//System.out.println(new Gson().toJson(model.reactions));
			//System.out.println(contents);
			
			// Try to add it to the cache - it's going to get read only once and very very soon... (TODO: check how reliable the cache is?)
			MemcacheService syncCache = MemcacheServiceFactory.getMemcacheService();
		    syncCache.setErrorHandler(ErrorHandlers.getConsistentLogAndContinue(Level.INFO));
		    for (int i=0; i<10; i++) { // UUID 10 times... if it doesn't get it by then, then something else is wrong!
		    	String tryId = UUID.randomUUID().toString().replaceAll("-", ""); // e.g "5ced1362fb634780b93ac05708d1ba02"
		    	if (!syncCache.contains(tryId)) {
		    		d.id = tryId;
		    		break;
		    	}
		    }
		    if (d.id == null) {
		    	throw new Exception("Failed to store item in cache");
		    }
		    syncCache.put(d.id, compressedContents, Expiration.byDeltaSeconds(60*10)); // they have 10 minutes to access it... should be instantly downloaded anyway
		    
		    // Finish response data, now that all looks okay
		    d.filenameHint = model.options.name + ".xml"; 
		    d.success = String.valueOf(true);
		    d.message = "Success";
		}
		catch (Exception e) {
			e.printStackTrace();
			d.success = String.valueOf(false);
			d.message = "Failed to prepare data for saving";
		}
		d.respondAsJson(resp);
	}
}