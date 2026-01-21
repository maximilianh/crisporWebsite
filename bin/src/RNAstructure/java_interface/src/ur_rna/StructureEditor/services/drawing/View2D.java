package ur_rna.StructureEditor.services.drawing;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.NoninvertibleTransformException;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

/**
 * Encapsulates the AffineTransform that converts from model coordinates to screen coordinates along with its
 * inverse transform and allows conversion of points between the two.
 * Also contains the bounds of the screen in both coordinate systems.
 */
public class View2D {
    public static final View2D IDENTITY = new View2D();
    public View2D() {
        trToScreen = new AffineTransform();
        trToModel = new AffineTransform();
        screenBounds = new Rectangle();
        modelBounds = new Rectangle2D.Float();
        screenMargin =  new Dimension(0,0);
    }
    public View2D(final AffineTransform transform, final Rectangle bounds) {
        this();
        setView(transform, bounds);
    }

    public void setView(View2D copyFrom) {
        trToScreen.setTransform(copyFrom.trToScreen);
        trToModel.setTransform(copyFrom.trToModel);
        screenBounds.setBounds(copyFrom.screenBounds);
        modelBounds.setRect(copyFrom.modelBounds);
        screenMargin.setSize(copyFrom.screenMargin);
    }
    public void setView(AffineTransform modelToScreenTransform, Rectangle screenBounds) {
        trToScreen.setTransform(modelToScreenTransform);
        trToModel.setTransform(modelToScreenTransform);
        try {
            trToModel.invert();
        } catch (NoninvertibleTransformException ex) {
            trToModel.setToTranslation(-trToScreen.getTranslateX(), -trToScreen.getTranslateY());
            trToModel.scale(1 / trToScreen.getScaleX(), 1 / trToScreen.getScaleY());
        }

        this.screenBounds.setRect(screenBounds);
        toModel(screenBounds, this.modelBounds);
    }

    /** Converts model coordinates to screen coordinates. */
    public final AffineTransform trToScreen;
    /** Converts screen coordinates to model coordinates. */
    public final AffineTransform trToModel;

    /** Represents the bounds of the screen (in screen) coordinates, so the location is always (0,0). */
    public final Rectangle screenBounds;
    /** Represents the bounds of the screen in model coordinates. */
    public final Rectangle2D modelBounds;

    public final Dimension screenMargin;

    public boolean isInScreenBounds(Point2D screenPoint) {
        return screenBounds.contains(screenPoint);
    }
    public boolean isInModelBounds(Point2D modelPoint) {
        return modelBounds.contains(modelPoint);
    }
    public boolean isInScreenBounds(Shape screenShape) {
        return screenBounds.contains(screenShape.getBounds2D());
    }
    public boolean isInModelBounds(Shape modelShape) {
        return modelBounds.contains(modelShape.getBounds2D());
    }

    public Point2D toScreen(Point2D p) { return trToScreen.transform(p, p); }
    public Point2D toModel(Point2D p) { return trToModel.transform(p, p); }
    public Point2D rayToScreen(Point2D p) { return trToScreen.deltaTransform(p, p); }
    public Point2D rayToModel(Point2D p) { return trToModel.deltaTransform(p, p); }

    public Point2D toScreen(Point2D p, Point2D pDest) { return trToScreen.transform(p, pDest); }
    public Point2D toModel(Point2D p, Point2D pDest) { return trToModel.transform(p, pDest); }
    public Point2D rayToScreen(Point2D p, Point2D pDest) { return trToScreen.deltaTransform(p, pDest); }
    public Point2D rayToModel(Point2D p, Point2D pDest) { return trToModel.deltaTransform(p, pDest); }

    private final Point2D.Double tmpPt1 = new Point2D.Double(), tmpPt2 = new Point2D.Double();
    public Rectangle2D toModel(Rectangle2D r) { return toModel(r, r); }
    public Rectangle2D toModel(Rectangle2D r, Rectangle2D rDest) {
        if (rDest == null) rDest = (Rectangle2D)r.clone();
        synchronized (tmpPt1) {
            tmpPt1.setLocation(r.getX(), r.getY());
            tmpPt2.setLocation(tmpPt1.x + r.getWidth(), tmpPt1.y + r.getHeight());
            rDest.setFrameFromDiagonal(toModel(tmpPt1), toModel(tmpPt2));
        }
        return rDest;
    }
    public Rectangle2D toScreen(Rectangle2D r) { return toScreen(r, r); }
    public Rectangle2D toScreen(Rectangle2D r, Rectangle2D rDest) {
        if (rDest == null) rDest = (Rectangle2D)r.clone();
        synchronized (tmpPt1) {
            tmpPt1.setLocation(r.getX(), r.getY());
            tmpPt2.setLocation(tmpPt1.x + r.getWidth(), tmpPt1.y + r.getHeight());
            rDest.setFrameFromDiagonal(toScreen(tmpPt1), toScreen(tmpPt2));
        }
        return rDest;
    }

    public int hashCode() {
        return trToScreen.hashCode() * 7247 ^ screenBounds.hashCode();
    }
}
