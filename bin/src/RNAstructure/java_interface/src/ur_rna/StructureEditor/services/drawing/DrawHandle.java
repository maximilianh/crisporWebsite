package ur_rna.StructureEditor.services.drawing;

import ur_rna.StructureEditor.models.DrawSettings;
import ur_rna.StructureEditor.models.Nuc;
import ur_rna.StructureEditor.models.SceneUpdateInfo;
import ur_rna.StructureEditor.services.SceneController;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.util.Collection;

import static ur_rna.Utilities.geom.PointMath.*;

/**
 * Represents a handle that the user can click, drag or interact with on the screen to cause an effect in the scene.
 */
public abstract class DrawHandle {
    public static final float SNAP_DISTANCE = 6;
    private static int defaultBoxSize = 16;
    protected final static BasicStroke defaultOutline = new BasicStroke(1);
    protected final static Color defaultLineColor = new Color(0xcc004d99, true);
    protected final static Color defaultBgColor = new Color(0xcc3399ff, true);
    private static int DEFAULT_MARGIN = 4; // minimum pixels between a nucleoside and a draw-handle box

    private boolean isPerformingDrag = false;
    private boolean enabled = true;
    protected Point location = new Point();
    protected Point startLocation = new Point(); // location at start of drag.
    protected Rectangle cachedBounds;
    protected Image icon;
    protected Cursor customCursor;
    protected Nuc targetNuc;
    protected Collection<Nuc> allNucs;
    protected SceneController controller;
    protected final AffineTransform _tmpTr = new AffineTransform();
    protected final Point2D.Double _tmpPt = new Point2D.Double();
    protected View2D view = View2D.IDENTITY;
    protected int boxSize = defaultBoxSize;

    protected DrawHandle(SceneController controller) { this.controller = controller; }

    /** Get the size of the default DrawHandle box that the user can click. */
    public int getBoxSize() {
        return boxSize;
    }
    /** Set the size of the default DrawHandle box that the user can click. */
    public void setBoxSize(final int size) {
        boxSize = size;
    }

    public Rectangle getBounds() {
        if (cachedBounds == null)
            cachedBounds = new Rectangle(0, 0, boxSize, boxSize);
        cachedBounds.setLocation(location.x - boxSize / 2, location.y - boxSize /2 );
        return cachedBounds;
    }
    public void setEnabled(boolean enabled) { this.enabled = enabled; }
    public boolean isEnabled() {
        return enabled;
    }
    public void setIcon(BufferedImage icon) {
        this.icon = icon;
    }
    public void settingsChanged() { }

    /**
     * Called when the user has completed dragging (e.g. on mouse-button-release.)
     * It is appropriate to add an undo event.
     * @param success Whether the action was successful (and therefore should be stored in the UNDO history).
     */
    public void dragComplete(final boolean success) {
        //isDragging = false;
        if (success) {
            SceneUpdateInfo info = getCompletionEvent();
            if (info != null)
                controller.addUndo(info);
        }
//        if (customCursor != null)
//            controller.getCanvas().getComponent().setCursor(Cursor.getDefaultCursor());
    }

    public void dragStarted(Point startPos) {
        //isDragging = true;
        startLocation.setLocation(location);
//        if (customCursor != null)
//            controller.getCanvas().getComponent().setCursor(customCursor);
    }

    public final Shape draw(final Graphics2D g, View2D view, DrawSettings settings, boolean isDragging) {
        this.view = view;
        prepareDraw(isDragging);
//        int viewHash = view.hashCode();
//        if (viewHash != prevViewHash) {
//            this.view = view;
//            prevViewHash = viewHash;
//            viewChanged();
//        }
        return drawBox(g, view, settings);
    }

    protected Shape drawBox(final Graphics2D g, View2D view, DrawSettings settings) {
        Rectangle rc = getBounds();
        if (icon != null)
            g.drawImage(icon, rc.x, rc.y, rc.width, rc.height, null);
        else {
            g.setColor(defaultBgColor);
            g.fill(rc);
            g.setStroke(defaultOutline);
            g.setColor(defaultLineColor);
            g.draw(rc);
        }
        return rc;
    }

    /**
     * Calculates the position of the DrawHandle's box based on a starting point and the suggested direction to move away from that point.
     * The returned Point is in device space, but parameters are in model space, so the view transform is are required to convert them.
     *
     * @param location  a Point2D specifying the starting point (e.g. a nucleotide location) in model space.
     * @param direction a Point2D that represents a vector that indicates the direction to move away from the location.
     *                  This should be given in model space. It will be normalized (to a unit vector) and then multiplied
     *                  by an appropriate margin.
     * @param margin    A margin given in model-space that will be converted to device space and augmented by an
     *                  additional device-space margin.
     */
    protected Point2D calcBoxPos(Point2D location, Point2D direction, double margin, final DrawSettings settings, final AffineTransform view) {
        Point2D loc = view.transform(location, null);   // get device coords.
        Point2D dir = view.deltaTransform(direction, null); // get device direction (which would only differ if there is a view rotation)
        if (margin != 0)
            margin *= view.getScaleX(); // scaleY would work also. we assume they are equal. modelMargin is now in device coords
        margin += boxSize / Math.sqrt(2) + defaultOutline.getLineWidth() + DEFAULT_MARGIN;
        translate(loc, scale(normalize(dir), margin)); // move away from the nucleotide by a distance equal to the nucleotide radius plus the size of the handle.
        translate(loc, -boxSize / 2, -boxSize / 2); // move from center to top/left
        return loc;
    }

//    protected AffineTransform getReverseDeltaTr(final Point2D deviceDelta, AffineTransform worldToDevice) {
//        try {
//            _tmpTr.setTransform(worldToDevice);
//            _tmpTr.invert();
//            _tmpTr.deltaTransform(deviceDelta, _tmpPt);
//        } catch (NoninvertibleTransformException ex) {
//            ex.printStackTrace();
//        }
//        _tmpTr.setToTranslation(_tmpPt.getX(), _tmpPt.getY());
//        return _tmpTr;
//    }

    public void setSelection(final Nuc focus, final Collection<Nuc> selection) {
        targetNuc = focus;
        allNucs = selection;
        selectionChanged();
//        int foc = targetNuc == null ? 0 : 1 + targetNuc.indexInScene();
//        int hash = getSelectionHashCode(selection);
//        if (foc != prevFocusIndex || hash != prevSelectionHash) {
//            prevFocusIndex = foc;
//            prevSelectionHash = hash;
//            selectionChanged();
//        }
    }

    private static int getSelectionHashCode(Iterable<Nuc> selection) {
        int hash = 0;
        for (Nuc n : selection)
            hash = (hash ^ n.indexInScene()) * 16777619;
        return hash;
    }

    protected AffineTransform getModelSpaceDeltaTransform(Point start, Point end) {
        _tmpPt.setLocation(end.x - start.x, end.y - start.y);
        view.rayToModel(_tmpPt);
        _tmpTr.setToTranslation(_tmpPt.x, _tmpPt.y);
        return _tmpTr;
    }

    protected void moveLocation(Point start, Point end) {
        location.setLocation(startLocation.x + end.x - start.x, startLocation.y + end.y - start.y);
    }

    /** Returns true when this DragHandle is performing updates in response to being dragged.
     * Specifically, it is true during a call to {@link #drag(Point, Point, Point, SceneController.DragOpts)}
     */
    public boolean isUpdating() { return isPerformingDrag; }

    /** Called when this handle has been dragged by the user. */
    public final void drag(Point start, Point prev, Point current, final SceneController.DragOpts options) {
        isPerformingDrag = true;
        try {
            performDrag(start, prev, current, options);
        } finally {
            isPerformingDrag = false;
        }
    }
    /** Called when any nucleotides in the scene changed position. */
    public final void layoutUpdated() {
        if (!isPerformingDrag) nucPositionsChanged();
    }
    /** Called when this handle has been dragged by the user. */
    protected abstract void performDrag(Point start, Point prev, Point current, final SceneController.DragOpts options);
    /** Called when any nucleotides in the scene changed position,
     * EXCEPT when THIS DrawHandle is in the middle of a performDrag operation
     * (because any nucleotide changes then are likely due to its own drag operation).
     */
    protected abstract void nucPositionsChanged();
    protected abstract void prepareDraw(final boolean isDragging);
    protected abstract void selectionChanged();
    public abstract boolean isValid();
    protected abstract SceneUpdateInfo getCompletionEvent();

}
