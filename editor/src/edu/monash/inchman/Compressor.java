// From: http://thiscouldbebetter.wordpress.com/2011/08/26/compressing-and-uncompressing-data-in-java-using-zlib/

package edu.monash.inchman;

import java.io.ByteArrayOutputStream;
import java.util.zip.DataFormatException;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

class Compressor {

	public static byte[] compressString(String str) {
		byte[] myByte = str.getBytes();

		Deflater compressor = new Deflater();
		compressor.setLevel(Deflater.BEST_COMPRESSION);
		compressor.setInput(myByte);
		compressor.finish();

		// Create an expandable byte array to hold the compressed data.
		// It is not necessary that the compressed data will be smaller than the
		// uncompressed data.
		ByteArrayOutputStream bos = new ByteArrayOutputStream(myByte.length);

		// Compress the data
		byte[] buf = new byte[1024];
		while (!compressor.finished()) {
			int count = compressor.deflate(buf);
			bos.write(buf, 0, count);
		}
		try {
			bos.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		return bos.toByteArray();
	}

	
	public static String decompressString(byte[] compressedBytes) {

		Inflater decompressor = new Inflater();
		byte[] buf = new byte[1024];
		decompressor.setInput(compressedBytes);

		// Create an expandable byte array to hold the decompressed data
		ByteArrayOutputStream bos = new ByteArrayOutputStream(compressedBytes.length);

		// Decompress the data
		buf = new byte[1024];
		while (!decompressor.finished()) {
			try {
				int count = decompressor.inflate(buf);
				bos.write(buf, 0, count);
			} catch (DataFormatException e) {
				e.printStackTrace();
			}
		}
		try {
			bos.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		return bos.toString();
	}
}