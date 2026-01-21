package ur_rna.StructureEditor.services.drawing;

import ur_rna.StructureEditor.Program;
import ur_rna.StructureEditor.models.Bond;
import ur_rna.StructureEditor.models.Motif;
import ur_rna.StructureEditor.models.Nuc;
import ur_rna.StructureEditor.models.SceneUpdateInfo;
import ur_rna.StructureEditor.services.SceneController;
import ur_rna.Utilities.geom.Rectangles;

import java.awt.*;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.Collection;

/**
 * A Drawhandle that allows the user to rotate a group of selected Nucleotides.
 */
public class SelectionRotater extends DrawHandle {
    private CenterPoint center;
    public SelectionRotater(final SceneController controller) {
        super(controller);
        center = new CenterPoint(controller, this);
        setIcon(Program.getImage("rotate-handle"));
        // customCursor = Toolkit.getDefaultToolkit().createCustomCursor(Program.getImage("rotate-cursor"), new Point(16, 16), "rotate");
    }

    public static class CenterPoint extends DrawHandle {
        private final SelectionRotater owner;
        private Point2D.Float modelLocation = new Point2D.Float();  // both in WORLD coordinates.
        private Point2D.Float previousCenter = new Point2D.Float();  // both in WORLD coordinates.
        private Nuc customNuc;
        private Motif.Helix customHelix;
        private RotationCenter center = RotationCenter.CenterOfMass;
        private int selectionHash;

        public CenterPoint(final SceneController controller, SelectionRotater owner) {
            super(controller);
            this.owner = owner;
            setIcon(Program.getImage("rotate-center"));
            customCursor = Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR);
        }
        @Override
        protected void performDrag(final Point start, final Point prev, final Point current, final SceneController.DragOpts options) {
            getModelSpaceDeltaTransform(prev, current).transform(modelLocation, modelLocation);
            controller.controlsUpdated();
        }

        @Override
        public SceneUpdateInfo getCompletionEvent() {
            return null;
        }

        @Override
        protected void prepareDraw(final boolean isDragging) {
            view.toScreen(modelLocation, location);
        }

        @Override
        protected void nucPositionsChanged() {
            if (!owner.isUpdating()) {
                // calculate the change in center
                Point2D.Float newCenter = new Point2D.Float();
                calcLocation(newCenter);
                modelLocation.setLocation(modelLocation.x + newCenter.x - previousCenter.x, modelLocation.y + newCenter.y - previousCenter.y);
                previousCenter.setLocation(newCenter);
            }
        }

        @Override
        protected void selectionChanged() {
            int prevHash = selectionHash;

            selectionHash = 987;
            for(Nuc n : allNucs)
                selectionHash = 31 * selectionHash + n.indexInScene();

            if (prevHash == selectionHash) return;

            reset();
            calcLocation(modelLocation);
            previousCenter.setLocation(modelLocation);

        }

        private void calcLocation(Point2D.Float p) {
            calcLocation(p, this.center);
        }

        private void calcLocation(Point2D.Float p, RotationCenter rc) {
            float x=0, y=0;
            switch (rc) {
                case Default: {
                    customHelix = controller.findBaseHelix(allNucs);
                    if (customHelix == null)
                        calcLocation(p, RotationCenter.CenterOfMass);
                    else
                        calcLocation(p, RotationCenter.BaseOfHelix);
                    return;
                }
                case CenterOfMass: {
                    for(Nuc n : allNucs) {
                        x += n.location.x;
                        y += n.location.y;
                    }
                    x /= allNucs.size();
                    y /= allNucs.size();
                } break;
                case CustomNuc:
                    x = customNuc.location.x; y = customNuc.location.y;
                    break;
                case CenterOfHelix: {
                    Collection<Nuc> bases = customHelix.getBases();
                    for (Nuc n : bases) {
                        x += n.location.x;
                        y += n.location.y;
                    }
                    x /= bases.size();
                    y /= bases.size();
                } break;
                case BaseOfHelix: {
                    Bond b = customHelix.getPair(0);
                    x = (b.getNuc5().location.x + b.getNuc3().location.x) / 2;
                    y = (b.getNuc5().location.y + b.getNuc3().location.y) / 2;
                    //System.out.println("customHelix: " + x + "," + y);
                } break;
            }
            p.setLocation(x, y);
        }

        private void reset() {
            center = RotationCenter.Default;
            customNuc = null;
            customHelix = null;
        }

        @Override
        public boolean isValid() {
            return allNucs != null && allNucs.size() > 1;
        }
        public void saveRelativeLocation() {
            calcLocation(previousCenter);
        }
    }

    public void setSelection(final Nuc focus, final Collection<Nuc> selection) {
        super.setSelection(focus, selection);
        center.setSelection(focus, selection);
    }

    @Override
    protected void nucPositionsChanged() {
        // no need here. handled in prepareDraw
    }

    @Override
    public void dragComplete(final boolean success) {
        super.dragComplete(success);
        center.saveRelativeLocation();
    }
    @Override
    public boolean isValid() {
        return allNucs != null && allNucs.size() > 1;
    }

    @Override
    protected void prepareDraw(final boolean isDragging) {
        if (!isDragging)
            resetLocation();
    }

    @Override
    protected void selectionChanged() {
      resetLocation();
    }

    void resetLocation() {
        // Attempt to find a "nice" location for the drag handle.
        // it should be top-right of the entire selection, but not outside of the view.
        if (!isValid()) return;
        Rectangle2D.Float bounds = new Rectangle2D.Float();
        float box = getBoxSize();
        if (view.screenBounds.isEmpty())
            bounds.setFrame(0, 0, box, box);
        else {
            boolean init = false;
            for (Nuc n : allNucs)
                if (init)
                    bounds.add(n.location);
                else {
                    bounds.setFrame(n.location.x, n.location.y, 0, 0);
                    init = true;
                }
            Rectangles.grow(bounds, controller.getSettings().nucleotideRadius(null));
            view.toScreen(bounds);
            Rectangles.grow(bounds, box);
            Rectangle screen = new Rectangle(view.screenBounds);
            screen.grow(-view.screenMargin.width, - view.screenMargin.height);
            Rectangle2D.intersect(bounds, screen, bounds); // get intersection of screen and model points.
        }
        location.setLocation(bounds.x + bounds.width - box / 2, bounds.y + box / 2);
    }

    /** Called when this handle has been dragged by the user. */
    @Override
    protected void performDrag(Point start, Point prev, Point current, final SceneController.DragOpts options) {
        moveLocation(start, current);
        controller.rotateNucs(allNucs, center.modelLocation, start, current, options);
    }

    @Override
    public SceneUpdateInfo getCompletionEvent() { return SceneUpdateInfo.Rotate; }

    public CenterPoint getCenterPoint() {
        return center;
    }

    public enum RotationCenter {
        /** Selects BaseOfHelix or CenterOfMass if there is no helix. */
        Default,
        CenterOfMass,
        CustomXY,
        CustomNuc,
        BaseOfHelix,
        CenterOfHelix;
    }
}
