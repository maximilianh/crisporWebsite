package ur_rna.StructureEditor.services.drawing;

import ur_rna.StructureEditor.models.IScreenObject;

import java.awt.*;
import java.awt.geom.Point2D;

/**
 * Represents an object drawn onscreen that is meant to be interacted with by the user.
 */
public class ScreenShape implements IScreenObject {
    public ScreenShape(final Object model, final Shape shape, final int zorder) {
        this.shape = shape;
        this.zorder = zorder;
        this.model = model;
    }
    protected Shape shape;
    protected int zorder;
    protected Object model;

    /**
     * Called before any other methods to notify the screen object of the current screen coordinates.
     *
     * @param view
     */
    @Override
    public void updateView(final View2D view) {
        // no need. shape is already in screen coords.
    }
    /**
     * Determine if the object contains the given point (in screen coordinates)
     *
     * @param point
     */
    @Override
    public boolean hitTest(final Point2D point) {
        return shape.contains(point);

    }
    /**
     * Determine if the object intersects the given rectangle (in screen coordinates)
     *
     * @param rc
     */
    @Override
    public boolean hitTest(final Rectangle rc) {
        return zorder != IScreenObject.ZORDER_MIN && shape.intersects(rc);
    }
    /**
     * @return Returns the bounds of this screen object (in screen coordinates)
     *         The bounds should NOT include any additional selection border or buffer zone.
     */
    @Override
    public Rectangle getBounds() {
        return shape.getBounds();
    }
    /**
     * @return Returns the z-order of the screen object.
     * Objects with higher z-orders are given preference in hit-testing.
     */
    @Override
    public int getZOrder() {
        return zorder;
    }

    @Override
    public Object getModel() {
        return model;
    }
}
