package ur_rna.StructureEditor.services;

import ur_rna.StructureEditor.Program;
import ur_rna.StructureEditor.Settings;
import ur_rna.StructureEditor.models.*;
import ur_rna.StructureEditor.services.drawing.DrawHint;
import ur_rna.StructureEditor.services.drawing.ICanvas;
import ur_rna.Utilities.geom.Vec2D;

import java.awt.*;
import java.awt.geom.Point2D;
import java.util.Collection;
import java.util.function.Consumer;

/**
 * Abstract some things away from RnaDrawController
 */
public interface SceneController {
    /**  flags passed to DrawHandle and dragNuc methods to indicate that the user is modifying the default behavior. */
    class DragOpts {
        public boolean expandMotif, disableHints, special;
        public DragOpts() {}
        public DragOpts(final boolean expand, final boolean special, final boolean hints) {
            expandMotif = expand;
            this.special = special;
            this.disableHints = !hints;
        }
    }

    void addUndo(SceneUpdateInfo event);
    boolean canUndo();
    boolean canRedo();
    void redo();
    void undo();
    void clearUndo();

    // Events that cause redraws but do not modify the scene
    void selectionUpdated();
    void viewUpdated();
    void historyUpdated();
    void controlsUpdated();

    // Events that indicate the scene has been modified
    void layoutUpdated(SceneUpdateInfo details);
    void layoutUpdated(SceneUpdateInfo details, boolean addHistoryEntry);
    void styleUpdated(SceneUpdateInfo details);
    void structureUpdated(SceneUpdateInfo details);
    void notifyUpdate(SceneUpdateInfo update);
    void notifyUpdate(SceneUpdateInfo update, boolean addHistoryEntry);

    boolean addSceneUpdateListener(Consumer<SceneUpdateEvent> listener);
    boolean removeSceneUpdateListener(Consumer<SceneUpdateEvent> listener);
    boolean addHistoryListener(Consumer<HistoryUpdateEvent> listener);
    boolean removeHistoryListener(Consumer<HistoryUpdateEvent> listener);

    ICanvas getCanvas();
    DrawSettings getSettings();
    void rotateNucs(Collection<Nuc> nucs, Point2D modelCenter, Point start, Point current, DragOpts options);
    void rotateNucs(Collection<Nuc> nucs, Point2D modelCenter, Point start, Point current, DragOpts options, Motif.Helix baseHelix, boolean showArc);
    void rotateNucs(Collection<Nuc> list, Point2D modelCenter, final double theta);
    Motif.Helix findBaseHelix(Collection<Nuc> range);

    DrawHint addHint(Shape shape);
    DrawHint addHint(String text, double centerX, double centerY);

    Vec2D calcNormalToNuc(Nuc cur);

    default Settings programSettings() { return Program.getInstance().settings(); }
}
