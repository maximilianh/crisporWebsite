package ur_rna.Utilities.swing;

import java.awt.*;
import java.awt.color.ColorSpace;
import java.awt.image.*;
import java.util.Hashtable;

/**
 * Image Utilities
 */
public class ImageUtil {
    /**
     * This method returns a buffered image with the contents of an image.
     * @param image Image to be converted
     * @param forceCopy This parameter has no effect unless the source image is a BufferedImage.
     *                  If false the source image will be returned as-is, since it is already a BufferedImage.
     *                  If true, a new BufferedImage (which is compatible with the screen) will be created
     *                  even if the source image is already a buffered image.
     * @return a buffered image with the contents of the specified image
     */
    public static BufferedImage toBufferedImage(Image image, boolean forceCopy) {
        if (image instanceof BufferedImage && !forceCopy) return (BufferedImage) image;

        // If the image is already a buffered image, then its pixels are already loaded.
        // Otherwise, we need to ensure that the image is completely loaded before continuing.
        if (!(image instanceof BufferedImage)) {
            // This code ensures that all the pixels in the image are loaded
            image = new javax.swing.ImageIcon(image).getImage();         // TODO: load image without using ImageIcon
        }
        // Create a buffered image with a format that's compatible with the screen
        BufferedImage bimage = createCompatibleImage(image.getWidth(null), image.getHeight(null), hasAlpha(image));

        // Copy image to buffered image
        Graphics g = bimage.createGraphics();
        // Paint the image onto the buffered image
        g.drawImage(image, 0, 0, null);
        g.dispose();
        return bimage;
    }
    /**
     * This method returns a buffered image with the contents of an image.
     * If the source image is already a BufferedImage, it will be returned as-is.
     * @param image Image to be converted
     * @return a buffered image with the contents of the specified image
     */
    public static BufferedImage toBufferedImage(Image image) { return toBufferedImage(image, false); }

    /**
     * Converts an arbitrary image to a {@code BufferedImage}.
     * @param image Image that should be converted.
     * @return a buffered image containing the image pixels, or the original
     *         instance if the image already was of type {@code BufferedImage}.
     */
    public static BufferedImage toBufferedImage(RenderedImage image) {
        if (image instanceof BufferedImage) return (BufferedImage) image;
        ColorModel cm = image.getColorModel();
        WritableRaster raster = cm.createCompatibleWritableRaster(image.getWidth(), image.getHeight());
        boolean isRasterPremultiplied = cm.isAlphaPremultiplied();
        Hashtable<String, Object> properties = null;
        if (image.getPropertyNames() != null) {
            properties = new Hashtable<String, Object>();
            for (String key : image.getPropertyNames()) {
                properties.put(key, image.getProperty(key));
            }
        }
        BufferedImage bimage = new BufferedImage(cm, raster, isRasterPremultiplied, properties);
        image.copyData(raster);
        return bimage;
    }

    public static BufferedImage getAlphaImage(BufferedImage image) {
        WritableRaster alphaRaster = image.getAlphaRaster();
        int width = image.getWidth();
        int height = image.getHeight();

        ColorModel cm;
        WritableRaster raster;
        // TODO Handle bitmap masks (work on ImageDataStream is necessary)
		/*
		if (image.getTransparency() == BufferedImage.BITMASK) {
            byte[] arr = {(byte) 0, (byte) 255};
            cm = new IndexColorModel(1, 2, arr, arr, arr);
            raster = Raster.createPackedRaster(DataBuffer.TYPE_BYTE,
            		width, height, 1, 1, null);
		} else {*/
        ColorSpace colorSpace = ColorSpace.getInstance(ColorSpace.CS_GRAY);
        int[] bits = {8};
        cm = new ComponentColorModel(colorSpace, bits, false, true,
                Transparency.OPAQUE, DataBuffer.TYPE_BYTE);
        raster = cm.createCompatibleWritableRaster(width, height);
        //}

        BufferedImage alphaImage = new BufferedImage(cm, raster, false, null);

        int[] alphaValues = new int[image.getWidth()*alphaRaster.getNumBands()];
        for (int y = 0; y < image.getHeight(); y++) {
            alphaRaster.getPixels(0, y, image.getWidth(), 1, alphaValues);
            // FIXME Don't force 8-bit alpha channel (see TODO above)
            if (image.getTransparency() == BufferedImage.BITMASK) {
                for (int i = 0; i < alphaValues.length; i++) {
                    if (alphaValues[i] > 0) {
                        alphaValues[i] = 255;
                    }
                }
            }
            alphaImage.getRaster().setPixels(0, y, image.getWidth(), 1, alphaValues);
        }

        return alphaImage;
    }

    /**
     * This method returns {@code true} if the specified image has the
     * possibility to store transparent pixels.
     * Inspired by http://www.exampledepot.com/egs/java.awt.image/HasAlpha.html
     * @param image Image that should be checked for alpha channel.
     * @return {@code true} if the specified image can have transparent pixels,
     *         {@code false} otherwise
     */
    public static boolean hasAlpha(Image image) {
        ColorModel cm;
        // If buffered image, the color model is readily available
        if (image instanceof BufferedImage) {
            BufferedImage bimage = (BufferedImage) image;
            cm = bimage.getColorModel();
        } else {
            // Use a pixel grabber to retrieve the image's color model;
            // grabbing a single pixel is usually sufficient
            PixelGrabber pg = new PixelGrabber(image, 0, 0, 1, 1, false);
            try {
                pg.grabPixels();
            } catch (InterruptedException e) {
                return false;
            }
            // Get the image's color model
            cm = pg.getColorModel();
        }
        return cm.hasAlpha();
    }


    /**
     * This method returns {@code true} if the specified image has at least one
     * pixel that is not fully opaque.
     * @param image Image that should be checked for non-opaque pixels.
     * @return {@code true} if the specified image has transparent pixels,
     *         {@code false} otherwise
     */
    public static boolean usesAlpha(Image image) {
        if (image == null) return false;
        BufferedImage bimage = toBufferedImage(image);
        Raster alphaRaster = bimage.getAlphaRaster();
        if (alphaRaster == null) {
            return false;
        }
        DataBuffer dataBuffer = alphaRaster.getDataBuffer();
        for (int i = 0; i < dataBuffer.getSize(); i++) {
            int alpha = dataBuffer.getElem(i);
            if (alpha < 255) {
                return true;
            }
        }
        return false;
    }
    /**
     * Returns a <code>BufferedImage</code> that supports the specified
     * transparency and has a data layout and color model
     * compatible with this <code>GraphicsConfiguration</code>.
     *
     * The returned <code>BufferedImage</code> has a layout and
     * color model that can be optimally blitted to a device
     * with this <code>GraphicsConfiguration</code>.
     * */
    public static BufferedImage createCompatibleImage(int width, int height, boolean supportTransparency) {
        try {
            // Create the buffered image
            GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
            GraphicsDevice gs = ge.getDefaultScreenDevice();
            GraphicsConfiguration gc = gs.getDefaultConfiguration();
            return gc.createCompatibleImage(width, height, supportTransparency?Transparency.TRANSLUCENT:Transparency.OPAQUE);
        } catch (HeadlessException e) {
            // The system does not have a screen
            // Create a buffered image using the default color model
            return new BufferedImage(width, height, supportTransparency ? BufferedImage.TYPE_INT_ARGB : BufferedImage.TYPE_INT_RGB);
        }
    }
}