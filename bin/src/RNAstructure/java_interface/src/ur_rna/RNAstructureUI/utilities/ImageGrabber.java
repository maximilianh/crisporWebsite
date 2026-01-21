/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.utilities;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;

/**
 * A class that creates images in various contexts for the RNAstructure GUI.
 *
 * @author Jessica S. Reuter
 */
public class ImageGrabber
	implements Serializable {
	private static final long serialVersionUID = 20120802;

	/**
	 * Create an image.
	 *
	 * @param path           The path to the image.
	 * @return               The image.
	 * @throws IOException   If an image cannot be read correctly. Since all
	 *                       image paths are hard coded, this should never
	 *                       be thrown.
	 */
	public static BufferedImage getImage( String path )
		throws IOException {
		InputStream iconIn = Resources.get(path);
		return ImageIO.read( iconIn );
	}

	/**
	 * Create an image.
	 *
	 * @param path           The path to the image.
	 * @return               The image, or null if an image cannot be read due to an IO error.
	 */
	public static BufferedImage tryGetImage( String path ) {
		try {
			InputStream iconIn = Resources.tryGet(path);
			return ImageIO.read(iconIn);
		} catch (IOException ex) {
			return null;
		}
	}

	/**
	 * Create an image icon.
	 *
	 * @param path           The path to the image.
	 * @return               The image icon.
	 * @throws IOException   If an image cannot be read correctly. Since all
	 *                       image paths are hard coded, this should never
	 *                       be thrown.
	 */
	public static ImageIcon getImageIcon( String path )
		throws IOException {
		return new ImageIcon( getImage( path ) );
	}
	/**
	 * Create an image icon.
	 *
	 * @param path           The path to the image.
	 * @return               The image icon, nor null if it is missing or cannot be read.
	 */
	public static ImageIcon tryGetImageIcon( String path ) {
		try {
			return new ImageIcon(getImage(path));
		} catch (IOException ex) {
			return null;
		}
	}


	/**
	 * Create an image label.
	 *
	 * @param path           The path to the image.
	 * @return               The image label.
	 * @throws IOException   If an image cannot be read correctly. Since all
	 *                       image paths are hard coded, this should never
	 *                       be thrown.
	 */
	public static JLabel getImageLabel( String path )
		throws IOException {
		return new JLabel( getImageIcon( path ) );
	}
}
