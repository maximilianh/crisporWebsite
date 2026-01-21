package ur_rna.StructureEditor.models;

import ur_rna.StructureEditor.services.drawing.View2D;

import java.awt.*;
import java.awt.geom.Point2D;

/**
 * Represents a model-space object that can be
 */
public interface IScreenObject {
    int ZORDER_MIN = Integer.MIN_VALUE;
    int ZORDER_MAX = Integer.MAX_VALUE;

    /**
     * Called before any other methods to notify the screen object of the current screen coordinates.
     */
    void updateView(View2D view);

    /**
     * Determine if the object contains the given point (in screen coordinates)
     */
    boolean hitTest(Point2D point);

    /**
     * Determine if the object intersects the given rectangle (in screen coordinates)
     * @param rc
     */
    boolean hitTest(Rectangle rc);

    /**
     * @return Returns the model-space object that is represented by this screen object.
     */
    Object getModel();

    /**
     * @return Returns the bounds of this screen object (in screen coordinates)
     */
    Rectangle getBounds();

    /**
     * @return Returns the z-order of the screen object.
     * Objects with higher z-orders are given preference in hit-testing.
     */
    int getZOrder();
}
