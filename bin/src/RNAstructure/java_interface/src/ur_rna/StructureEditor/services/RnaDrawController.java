package ur_rna.StructureEditor.services;

import ur_rna.StructureEditor.Settings;
import ur_rna.StructureEditor.models.*;
import ur_rna.StructureEditor.services.drawing.*;
import ur_rna.StructureEditor.services.drawing.export.Graphics2DRecorder;
import ur_rna.Utilities.EventSource;
import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.annotation.NotNull;
import ur_rna.Utilities.geom.Ellipses;
import ur_rna.Utilities.geom.GraphicsUtil;
import ur_rna.Utilities.geom.Vec2D;
import ur_rna.Utilities.swing.AcceleratorKey;
import ur_rna.Utilities.swing.InputAdapter;

import javax.swing.*;
import java.awt.*;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.geom.*;
import java.util.*;
import java.util.List;
import java.util.function.Consumer;

import static ur_rna.Utilities.geom.PointMath.diff;
import static ur_rna.Utilities.geom.Rectangles.fromPoints;

/**
 * The controller that interprets input from the screens, modifies the RNAscene,
 * keeps a history of changes, and reports updates to subscribed listeners.
 */
public class RnaDrawController implements SceneRenderer, SceneController {
    public void programSettingsChanged(final Settings settings) {
        handleLoopResize.setEnabled(!settings.DisableLoopResize);
        handleBranchSlide.setEnabled(!settings.DisableBranchSlide);
        handleRotate.setEnabled(!settings.DisableRotation);
        ((SelectionRotater)handleRotate).getCenterPoint().setEnabled(!settings.DisableRotation);
    }
    /** BehaviorMode defines two separate "behavior modes"
     * that combine selection behavior along with draw-handle and dragging behaviors.
     * ***********NOTE: NOT YET IMPLEMENTED***********
     *  ### Branch Slide Mode
     *      - Click – Selects Branch
     *      - Drag – Slides Branch (or Rotates helix if no slide is possible)
     *      - Shift-Drag – Moves Branch (with slide-hints)
     *      - Alt-Drag – Rotates Branch (at base)
     *      - Loop-Handle allow Loop resizing.
     *      - Press ‘R’ or select with band to show Rotation handles.
     *  ### Nuc Mode
     *      - Click Selects Nucleotide (Press spacebar to select Helix/Loop. Again to select branch)
     *      - Handles allow Loop resizing and arbitrary rotation.
     *      - Drag – Moves Selection
     *      - Shift-Drag – Rotates Selection (at base if helix is selected)
     */
    public enum BehaviorMode { NUC, BRANCH }

    private Component canvasComponent;
    private ICanvas canvas; // same object as canvasComponent, but via ICanvas interface.
    private RnaScene scene;
    private Set<Nuc> selection = new HashSet<>();
    private Collection<DrawHint> drawHints = new HashSet<>();
    private final DrawHandle handleRotate, handleBranchSlide, handleLoopResize;
    private final DrawHandle[] drawHandles;

    private Nuc focused;

    private DrawSettings settings = new DrawSettings();
    private MouseHandler mouseHandler = new MouseHandler();

    private List<IScreenObject> screenObjects = null;

    public EventSource.OneArg<SceneUpdateEvent> sceneUpdated = new EventSource.OneArg<>();
    public EventSource.OneArg<HistoryUpdateEvent> historyEvent = new EventSource.OneArg<>();

    public static int SELECTION_ADD_MASK = InputEvent.SHIFT_MASK,
            EDIT_CLICK_MASK = InputEvent.CTRL_MASK,
            DRAG_SPECIAL_MASK = AcceleratorKey.SystemMenuKeyModifier, // user can hold down a modifier key to affect the action of dragging a draw handle or nucleotide. Typically this simply disables sticky mode.
            DRAG_EXPANDED_MOTIF_MASK = InputEvent.SHIFT_MASK; // user can hold down a modifier key to affect the action of dragging a draw handle or nucleotide. Typically this expands the selection to include a larger structural motif (e.g. the full branch)


    private final Rectangle selectionBand = new Rectangle();
    private final Line2D.Float editBond = new Line2D.Float();
    private SelectionType selType = SelectionType.Individual;

    private BehaviorMode behaviorMode = BehaviorMode.NUC; // one of NUC_BEHAVIOR_MODE or BRANCH_BEHAVIOR_MODE

    public DrawSettings getSettings() { return settings; }
    public void setSettings(final DrawSettings settings) { this.settings = settings; }

    public RnaDrawController() {
        setScene(createDefaultScene());
        handleBranchSlide = new BranchSlider(this);
        handleRotate = new SelectionRotater(this);
        handleLoopResize = new LoopResizer(this);
        drawHandles = new DrawHandle[] { handleRotate, ((SelectionRotater)handleRotate).getCenterPoint(), handleLoopResize, handleBranchSlide };
        sceneUpdated.add(this::onSceneUpdate);
    }

//    char[] _selStrBuilder = new char[16];
//    final static int BITS_PER_CHAR = 16; // 0xFFFF => [0000 0000 0000 0000]
//    private String getSelStr() {
//        char[] arr = _selStrBuilder;
//        int used = 16;
//        while (used-- > 0)
//            arr[used]=0;
//        for (Nuc n : selection) {
//            int pos = n.index >> BITS_PER_CHAR;
//            int val = n.index % BITS_PER_CHAR;
//            if (pos > used) used = pos;
//            arr[pos] |= val;
//        }
//        return new String(arr, 0, used);
//    }

    private int[] EMPTY_SELECTION = new int[0];
    private int[] getSelectionState() {
        if (selection.size() == 0) return EMPTY_SELECTION;
        int i = 0;
        int[] sel = new int[selection.size() + 1];
        sel[0] = focused == null ? -1 : focused.indexInScene();
        for(Nuc n : selection)
            sel[++i] = n.indexInScene();
        return sel;
    }
    private void setSelectionState(int[] sel) {
        selection.clear();
        focused = null;
        if (sel != null && sel.length != 0) {
            if (sel[0] != -1) focused = scene.nucAtSceneIndex(sel[0]);
            for (int i = 1; i < sel.length; i++)
                selection.add(scene.nucAtSceneIndex(sel[i]));
        }
        selectionUpdated();
    }

    private RnaScene createDefaultScene() {
        RnaScene scene  = new RnaScene();
        Strand s = scene.strands.first();
        String seq = "AGUCUCCGCAGAUCAG";
        float r = settings.nucleotideRadius(null);
        float spacing = r * 2 * 1.4f;
        float margin = r * 1.5f;
        for(int i = 0; i < seq.length(); i++)
            s.add(new Nuc(""+seq.charAt(i), margin + spacing*i, margin));
        return scene;
    }


    public ICanvas getCanvas() {
        return canvas;
    }
    public void setCanvas(final ICanvas newCanvas) {
        canvas = newCanvas;
        setCanvasComponent(canvas.getComponent());
    }
    private void setCanvasComponent(final Component canvasComponent) {
        if (canvasComponent != null)
            mouseHandler.remove(canvasComponent);
        if (canvasComponent == null)
            throw new IllegalArgumentException("Canvas cannot be null.");
//        if (!(canvasComponent instanceof ICanvas))
//            throw new IllegalArgumentException("Canvas must implement the ICanvas interface.");
        this.canvasComponent = canvasComponent;
        mouseHandler.listen(canvasComponent);
    }
    public RnaScene getScene() {
        return scene;
    }
    public void setScene(final RnaScene scene) {
//        if (scene != null)
//            scene.updated.remove(this::onSceneUpdate);
        if (scene == null)
            throw new IllegalArgumentException("Scene cannot be null.");
        this.scene = scene;
        if (scene.history.size() == 0)
            addUndo(SceneUpdateInfo.Initial);
        //scene.updated.add(this::onSceneUpdate);
        resetStructure();
    }

    public boolean isUserActionInProgress() {
        return mouseHandler.state == MouseHandler.DRAGGING;
    }

    //@Override
    protected void onSceneUpdate(final SceneUpdateEvent eventInfo) {
        switch (eventInfo.type) {
            case Selection:
                handleSelectionChange();
                break;
            case Layout:
            case Style:
                break;
            case Structure:
                resetStructure();
                break;
            case Controls:
            case View:
                break;
            case History:
                return; //no need to redraw. Scene was not updated. Just history.
        }
        // Redraw for all events, except history
        redraw();
    }
    private void resetStructure() {
        // TODO: verify selection instead of clearing it.
        clearSelection();
        nucShapes = new Shape[scene.getNucCount()];
        nucScreenShapes = new ScreenNuc[nucShapes.length];
        scene.identifyPseudoknots();
    }
    private void clearSelection() {
        selection.clear();
        //disableDrawHandles();
        drawHints.clear();
        focused = null;
        selectionUpdated();
    }

    private void redraw() {
        if (_suspendRedraw==0 && canvas != null)
            canvas.repaintRequired();
    }

    Graphics2DRecorder _dummyGraphics = new Graphics2DRecorder();
    public Rectangle calcBounds(final Graphics2D graphics, View2D view, boolean includeInteractive) {
        _dummyGraphics.suspendRecording(); // turns off actual recording, so only the graphics context is updated.
        try {
            _dummyGraphics.getContext().copyFrom(graphics);
            List<IScreenObject> objs = draw(_dummyGraphics, view, this.settings, includeInteractive);
            int size = objs.size();
            if (size == 0) return null;
            Rectangle bounds = objs.get(0).getBounds();
            for (IScreenObject o : objs)
                bounds = bounds.union(o.getBounds());
            return bounds;
        } finally {
            _dummyGraphics.resumeRecording();
        }
    }

    private boolean useRecorder  = false;
    public Collection<IScreenObject> render(final Graphics2D graphics, View2D view, boolean includeInteractive) {
        if (useRecorder) {
            Graphics2DRecorder gr = new Graphics2DRecorder();
            gr.getContext().copyFrom(graphics);
            //gr.setClip(view.screenBounds);
            screenObjects = draw(gr, view, settings, includeInteractive);
            gr.replay(graphics);
        } else
            screenObjects = draw(graphics, view, settings, includeInteractive);

        return screenObjects;
    }

    private final AffineTransform _tmpTr = new AffineTransform();
    private final Point2D _tmpPt = new Point2D.Double();
    /**
     * Translates all selected nucleotides by the x and y coordinates of the delta point,
     * which is assumed to be in screen coordinates.
     */
    public void translateScreenNucs(Collection<Nuc> list, final Point start, final Point end) { translateScreenNucs(list, diff(end, start)); }
    /**
     * Translates all selected nucleotides by the x and y coordinates of the delta point,
     * which is assumed to be in screen coordinates.
     */
    public void translateScreenNucs(Collection<Nuc> list, final Point2D screenDelta) {
        canvas.getView().rayToModel(screenDelta, _tmpPt);
        _tmpTr.setToTranslation(_tmpPt.getX(), _tmpPt.getY());
        transformNucs(list, _tmpTr);
        //translateObjects(list, p.x, p.y);
    }

    private Point2D prevStart = new Point(); double prevTheta; private boolean prevBias;
    private final static double minBiasRange = 2 * Math.PI / 3;
    private final static double[] stickyRotationAngles = new double[] {
        0, Math.PI, -Math.PI,
            Math.PI/2, -Math.PI/2,
            Math.PI/4, -Math.PI/4,
            3 * Math.PI/4, -3 * Math.PI/4,
            Math.PI /3, 2 * Math.PI /3,
            -Math.PI /3, -2 * Math.PI /3
    }; private final static double stickyRotationEpsilon = Math.PI / 36; // 5 degrees
    @Override
    public void rotateNucs(Collection<Nuc> list, Point2D modelCenter, final Point start, final Point end, DragOpts options) { rotateNucs(list, modelCenter, start, end, options, null, true); }
    @Override
    public void rotateNucs(Collection<Nuc> list, Point2D modelCenter, final Point start, final Point end, DragOpts options, Motif.Helix baseHelix, boolean showArc) {
        if (list.size() == 0) return;
        if (baseHelix == null)
            baseHelix = findBaseHelix(list);

        Vec2D basePoint = null, orientation = null;
        if (baseHelix != null) {
            Bond baseBond = baseHelix.getPair(0);
            basePoint = Vec2D.getMidpoint(baseBond.n5.location, baseBond.n3.location);
            Bond endBond = baseHelix.getPair(baseHelix.size()-1);
            Vec2D endPoint = Vec2D.getMidpoint(endBond.n5.location, endBond.n3.location);
            orientation = new Vec2D(basePoint, endPoint);
        }

        if (modelCenter == null) {
            if (basePoint == null)
                modelCenter = calcNucGroupCenter(list);
            else
                modelCenter = basePoint;
        }

        final Point2D center = canvas.getView().toScreen(modelCenter, new Point2D.Float()); // this needs to be in screen coords for math with startPt and endPt (which are also in screen coords)
        Vec2D rayStart = new Vec2D(center, start), rayEnd = new Vec2D(center, end);

        boolean sticky = false;
        if (orientation != null && (!showArc || modelCenter.distanceSq(basePoint) < 6)) {
            double aDiff = -orientation.angleTo(rayStart);
            rayStart.rotate(aDiff);
            rayEnd.rotate(aDiff);
        }
        if (programSettings().SnapGuides ^ options.special /* "special" disables stickies*/) {
            double a = rayEnd.getAngle();
            for (double d : stickyRotationAngles)
                if (Math.abs(d - a) < stickyRotationEpsilon) {
                    rayEnd.rotate(d - a);
                    sticky = true;
                }
        }


        double thetaStart = rayStart.getAngle(), thetaEnd = rayEnd.getAngle();
        double theta = thetaEnd - thetaStart;

        // Assume that thetaStart is 90 degrees thetaEnd is 135 degrees.
        // The user could have gotten there by rotating -225 degrees (clockwise) or 45 degrees (counter-clockwise).
        // We would like to keep track of the user's actual path so that we can draw the appropriate arc.
        // In order to track the path, we store the previous theta and detect when the sign changes.
        // Note: -PI <= thetaStart <= PI  and -PI <= thetaEnd <= PI
        // So for theta (i.e. thetaEnd - thetaStart):  -2PI <= theta <= 2PI
        // If the user's mouse crosses over rayStart at the "top" of the circle (near the start point -- or more specifically on the same side of the center-point as the start point)
        //   then abs(theta) will be small (because thetaEnd will be slightly greater or less than theta.
        // If the user's mouse crosses over rayStart at the "bottom" of the circle (on the opposite side of the center-point than the start point)
        //   then abs(theta) will be large. In this case, the crossover should be ignored (because the direction of travel has not changed, even if the sign of theta has)
        // e.g. user crosses from 85 to 95 -- small change means user crossed over rayStart (which is at 90) and the direction of travel has changed -- it was clock and changed to anti-.
        // e.g. user crosses from 275 to 265 -- small change means user crossed over rayStart (which is at 90) and the direction of travel has changed -- it was clock and changed to anti-.


        boolean setBias;
        if (prevStart.equals(start)) {
            setBias = Math.abs(prevTheta) < minBiasRange &&
                    Math.abs(theta) < minBiasRange &&
                    (theta < 0 != prevTheta < 0);
        } else {
            // if prevStart is different, it means we are working with an entirely new rotation. so discard any previous bias
            setBias = true;
            prevStart = start;
        }
        prevTheta = theta;

        if (setBias) {
            prevBias = theta > 0 && theta <= Math.PI;
            //System.out.println("SETTING BIAS: " + prevBias);
        }


        if (prevBias && theta < 0)
            theta += 2 * Math.PI;
        else if (!prevBias && theta > 0)
            theta -= 2 * Math.PI;
//        if (theta - theta1 > Math.PI) theta -= 2 * Math.PI;
//        if (theta - theta1 < Math.PI) theta += 2 * Math.PI;

        // System.out.println(fmt("Theta: %s=%s\t%s=%s", theta1,theta1 * -180 / Math.PI, theta,theta * -180 / Math.PI));

//        System.out.println("Rotate theta: " + theta + " = " + (180 * theta / Math.PI));
        int centerRad = 5;
        int arcRadius = (int)(Math.max(rayStart.length(), rayEnd.length()) + 15) + 2 * centerRad;

        if (!options.disableHints) {
            addHint(Ellipses.fromCenter(center, centerRad)).fill(DrawHint.HintColorReference).noLine(); // the dot at the center of the rotation.
            //drawHints.add(new DrawHint(new Line2D.Double(center.getX(), center.getY(), start.getX(), start.getY()), DrawHint.cold, true));
            //drawHints.add(new DrawHint(new Line2D.Double(center.getX(), center.getY(), end.getX(), end.getY()), DrawHint.cold, true));
            if (showArc)
                addHint(new Arc2D.Double(center.getX() - arcRadius, center.getY() - arcRadius, 2 * arcRadius, 2 * arcRadius, thetaStart * -180 / Math.PI, theta * -180 / Math.PI, Arc2D.PIE));
            else
                addHint(rayEnd.toLine(center, arcRadius));
            if (sticky)
                addHint(rayEnd.toLine(center, arcRadius)).bold();
            //center.getX() + Math.cos(thetaStart + theta / 2) * (arcRadius + 15), center.getY() +Math.sin(thetaStart + theta / 2) * (arcRadius + 15)
            Point angleDisplay = rayEnd.getPointAlong(center, arcRadius + 20).toPoint();
            addHint(Math.round(thetaEnd * -180 / Math.PI) + " \u00B0", angleDisplay.x, angleDisplay.y);
        }

        if (theta != 0)
            rotateNucs(list, modelCenter, theta);
    }

    public DrawHint addHint(Shape shape) { DrawHint d = DrawHint.shape(shape); drawHints.add(d); return d; }
    public DrawHint addHint(String text, double centerX, double centerY) { DrawHint d = DrawHint.text(text, new Point2D.Double(centerX, centerY)); drawHints.add(d); return d; }

    public void rotateNucs(Collection<Nuc> list, Point2D modelCenter, final double theta) {
        _tmpTr.setToRotation(theta, modelCenter.getX(), modelCenter.getY());
        layoutUpdated(SceneUpdateInfo.Rotate, false);
        transformNucs(list, _tmpTr);
    }

    protected void transformNucs(Collection<Nuc> list, AffineTransform tr) {
        if (list.size() == 0) return;
        for (Nuc n : list)
            tr.transform(n.location, n.location);
        layoutUpdated(SceneUpdateInfo.Layout, false);
    }

    public Motif.Helix findBaseHelix(final Collection<Nuc> list) {
        for(Motif.Helix h : Motif.Helix.findAll(scene))
            if (list.contains(h.getPair(0).n5) && list.contains(h.getPair(0).n5))
                return h;
        return null;
    }

    private Point2D.Float calcNucGroupCenter(final Collection<Nuc> list) {
        Iterator<Nuc> i = list.iterator();
        if (!i.hasNext())
            return null;
        Nuc n = i.next();
        int count = 1;
        Point2D.Float p = new Point2D.Float(n.location.x, n.location.y);
        while(i.hasNext()) {
            n = i.next();
            count++;
            p.x += n.location.x;
            p.y += n.location.y;
        }
        p.x = p.x / count;
        p.y = p.y / count;
        return p;
    }
    public void selectSpecial(final SelectionType type) {
        if (focused == null) return;
        selectMotif(focused, true, type);
        selectionUpdated();
    }
    public void setSelectedBondType(final BondType selectedBondType) {
        for (Nuc n : selection) {
            if (n.isPaired())
                n.getPairBond().type = selectedBondType;
        }
        structureUpdated(SceneUpdateInfo.BondType.subType(String.format("Set Basepair Type (%s)", selectedBondType.getDescription())));
    }

    public void splitStrandsAtSelection() {
        if (selection.isEmpty()) return;
        for (Nuc n : getSelected()) {
            Nuc next = n.getNext();
            if (next != null)
                scene.divideStrand(next);
        }
        structureUpdated(SceneUpdateInfo.StrandDivisions.subType("Divided Strands"));
    }
    public void joinSelectionStrands() {
        if (selection.isEmpty()) return;
        Nuc[] sel = getSelected();
        boolean joined = false;
        for (int i = 0; i < sel.length; i++) {
            int s1i = sel[i].getStrand().getIndex();
            for (int j = i+1; j < sel.length; j++) {
                if (Math.abs(s1i-sel[j].getStrand().getIndex())==1) {
                    scene.joinStrands(sel[i].getStrand(), sel[j].getStrand());
                    joined = true;
                }
            }
        }
        if (joined)
            structureUpdated(SceneUpdateInfo.StrandDivisions.subType("Joined Strands"));
    }
    public void breakSelectedBonds() {
        if (selection.isEmpty()) return;
        for(Nuc n : getSelected())
            if (n.isPaired())
                scene.breakBond(n.getPairBond());
        structureUpdated(SceneUpdateInfo.Bonds.subType("Removed Basepairs"));
    }
    public Nuc[] getSelected() {
        return selection.toArray(new Nuc[selection.size()]);
    }
    public void deleteSelectedBases() {
        if (selection.isEmpty()) return;
        scene.removeNucs(getSelected());
        structureUpdated(SceneUpdateInfo.SequenceSize.subType("Removed Bases from Sequence"));
    }
    public void rotateNumber() {
        if (selection.isEmpty()) return;
        for(Nuc n : getSelected())
            n.style().numberOffset = (float)Vec2D.mod2PI(n.style().numberOffset + Math.PI / 8);
        structureUpdated(SceneUpdateInfo.Layout.subType("Rotated Base Number"));
    }
    private void rotateNumber(final Nuc nuc, final Point2D pos) {
        Vec2D norm = calcNormalToNuc(nuc);
        Vec2D dir = new Vec2D(nuc.location, pos);
        nuc.style().numberOffset = (float)norm.angleTo(dir);
    }

        /** Insert new bases before or after the specified nucleotide. If insertAt is null, the focused nucleotide will be used. **/
    public void insertBases(String sequence, Nuc ref, boolean insertAfter) {
        if (ref == null) ref = focused;
        if (ref == null) ref = ObjTools.first(selection, null);
        if (ref == null) throw new UnsupportedOperationException("Insert position was unspecified and no bases were selected.");
        Strand strand = ref.getStrand();
        Nuc adj = insertAfter ? ref.getNext() : ref.getPrev(); // get the adjacent base (for positioning later)

        int pos = ref.indexInStrand();
        if (insertAfter) pos++;
        int count = 0;
        for(int i = 0; i < sequence.length(); i++) {
            String symbol = sequence.substring(i, i+1);
            if (Character.isWhitespace(symbol.charAt(0)))
                continue;
            Nuc nuc = new Nuc(symbol);
            strand.add(pos+count++, nuc);
        }
        // now position the new nucleotides
        float len = settings.calcOptimalLoopDistance(count, 1); // length required to fit the new bases
        float radius = (float)(len / Vec2D.TwoPI); // required radius of circle to fit bases
        Vec2D center, start;
        if (adj == null) {
            // adj is null, which means ref was at the beginning (insertAfter=false) or end (insertAfter=true) of the
            // strand, so position nucs at start or end of strand.
            // First find the nuc on the OTHER side
            adj = insertAfter ? ref.getPrev() : ref.getNext();
            Vec2D dir;
            if (adj == null)
                dir = new Vec2D(1, 0);
            else
                dir = new Vec2D(adj.location, ref.location);
            start = dir.getPointAlong(ref.location, 2*settings.nucRadius+settings.loopSpacing);
            center = dir.getPointAlong(start, radius);
        } else {
            // position new nucs between two existing nucs
            Nuc n1 = insertAfter ? ref : adj, n2 = insertAfter ? adj : ref;
            Vec2D dir = new Vec2D(n1.location, n2.location);
            dir.rotate90(); // dir is the vector perpendicular to the ray from ref to adj.
            if (scene.drawFlipped)
                dir.negate();
            Vec2D mid = Vec2D.getMidpoint(n1.location, n2.location);
            start = dir.getPointAlong(n1.location, 2*settings.nucRadius+settings.loopSpacing);
            center = dir.getPointAlong(mid, 2*settings.nucRadius+settings.loopSpacing+radius);
        }
        Vec2D ray = new Vec2D(center, start);
        int dir = scene.drawFlipped ? 1 : -1;
        for(int i = 0; i < count; i++) {
            strand.get(pos + i).location.setLocation(ray.getPointAlong(center, radius));
            ray.rotate(dir * Vec2D.TwoPI / count);
        }
        structureUpdated(SceneUpdateInfo.SequenceSize.subType("Inserted Bases into Sequence"));
    }
    public void setBehaviorMode(final BehaviorMode behaviorMode) {
        this.behaviorMode = behaviorMode;
    }
    public BehaviorMode getBehaviorMode() {
        return behaviorMode;
    }
    public void removeColors(final boolean allNucs) {
        Collection<Nuc> nucs = allNucs ? scene.allNucs() : selection;
        for(Nuc n : nucs)
            if (n.style != null) {
                n.style.outlineColor = n.style.fillColor = n.style.textColor = null;
            }
        styleUpdated(SceneUpdateInfo.FormatBases.subType("Removed Color Annotations"));
    }
    // Reflect nuc positions throw a mirror line.
    public void flip(final boolean vertical, final boolean selectionOnly) {
        Collection<Nuc> nucs = selectionOnly ? selection : scene.allNucs();
        Point2D.Float center = calcNucGroupCenter(nucs);
        if (center == null) return; //return if nucs is empty
        for (Nuc n : nucs) {
            if (vertical)
                n.location.x = center.x - (n.location.x - center.x);
            else
                n.location.y = center.y - (n.location.y - center.y);
            if (n.style!=null && n.style.numberOffset != 0F) 
                n.style.numberOffset = -n.style.numberOffset;
        }
        scene.drawFlipped = !scene.drawFlipped;
        layoutUpdated(SceneUpdateInfo.Layout.subType(String.format("Flipped %s %s", selectionOnly ? "selection" : "scene", vertical ? "vertically" : "horizontally" )));
    }

    public void colorize(SceneColorizer sc, final boolean allNucs) {
        Collection<Nuc> nucs = allNucs ? scene.allNucs() : selection;
        if (sc.color(nucs))
            styleUpdated(SceneUpdateInfo.FormatBases.subType("Applied Color Annotations"));
    }

    public void redrawSelection() {
        NucLayout layout = new NucLayout(settings);
        int modified = 0;
        for (Motif.Helix h : scene.getHelices())
            for (Nuc n : h.getBases())
                if (selection.contains(n)) {
                    modified++;
                    layout.redrawHelix(h);
                    break; // back to outer helix loop
                }
        for (Motif.MultiLoop m : scene.getMultiLoops(1, false)) {
            //System.out.println(m.toString());
            for (Nuc n : m.getBases())
                if (selection.contains(n)) {
                    modified++;
                    layout.redrawMultiLoop(m, this);
                    break; // back to outer loop
                }
        }
        if (modified != 0)
            layoutUpdated(SceneUpdateInfo.Layout.subType("Redraw Selected Helices and Loops"));
    }
    public Nuc getFocused() {
        return focused;
    }

//    protected void translateObjects(Collection<? extends IDrawable> list, final float dx, final float dy) {
//        for (IDrawable d : list)
//            d.translate(dx, dy);
//        if (list.size() != 0)
//            scene.layoutUpdated();
//    }

    private class MouseHandler extends InputAdapter.Mouse {
        final static int SELECTING = 1, DRAGGING = 2, EDITING = 3, FORMATTING = 4;

        boolean dragIsValid;
        //final static int minDragDistSq = 0;//8; // i.e. at least 2px in both x and y OR 3 px in just x or just y
        int state; // SELECTING, DRAGGING, EDITING
        int prevModifiers;
        Point prevPos, curPos, startPos;
        Object dragObject = null;
        Nuc hoverNuc = null;
        Nuc deselectNuc;

        void storeMouse(MouseEvent e) {
            startPos = prevPos = curPos = e.getPoint();
            dragIsValid = false;
            deselectNuc = null;
            prevModifiers = e.getModifiers();
        }
        void updateMousePos(MouseEvent e) {
            prevPos = curPos;
            curPos = e.getPoint();
        }

        /**
         * Respond to mouse input.
         *
         * @param type The type of mouse input (e.g. click, move, etc)
         * @param e    A MouseEvent with additional information about the event.
         */
        @Override
        protected void onMouseInput(final InputAdapter.InputType type, final MouseEvent e) {
            //Component c = e.getComponent();
//                if (c != prevComponent) {
//                    screen = screens.get(c);
//                    clearMouse();
//                    prevComponent = c;
//                }
            int mod = e.getModifiers();

            // If the user presses or releases a modifier key (e.g. SHIFT) while they are
            // dragging a Nucleotide or DrawHandle, we interpret this as completing the previous operation and
            // starting a new one.
            if (state == DRAGGING && prevModifiers != mod && type == InputAdapter.InputType.Drag && !(dragObject instanceof DrawHandle)) {
                completeAction(prevModifiers);
                storeMouse(e);
                storeRnaCoords();
                //prevModifiers = mod;
//                MouseEvent up = new MouseEvent((Component)e.getSource(), e.getID(), e.getWhen(), prevModifiers, prevPos.x, prevPos.y, 0, false, e.getButton());
//                MouseEvent dn = new MouseEvent((Component)e.getSource(), e.getID(), e.getWhen(), mod, e.getX(), e.getY(), 0, false, e.getButton());
//                onMouseInput(InputAdapter.InputType.MouseUp, up);
//                onMouseInput(InputAdapter.InputType.MouseDown, dn);
            }

            switch (type) {
                case Move: {
                    //****************
                    // Mouse-Move -- Mouse is moving over the canvas but no mouse buttons are pressed.
                    // I.e. we are not dragging anything. So this is only for mouse-over or mouse-hover effects.
                    //****************
                    updateMousePos(e);
                    IScreenObject hit = hitTest(screenObjects, curPos, true);
                    if (hit == null) {
                        hoverNuc = null;
                        ((JComponent) canvasComponent).setToolTipText(null);
                    } else {
                        Nuc n = (Nuc) hit.getModel();
                        if (n == hoverNuc) break;
                        hoverNuc = n;
                        StringBuilder sb = new StringBuilder();
                        sb.append(n.symbol).append(' ').append(n.indexInStrand() + 1);
                        if (!n.isPaired())
                            sb.append(" Loop: ").append(n.getLoop().getType().name());
                        ((JComponent) canvasComponent).setToolTipText(sb.toString());
                    }
                } break;
                case MouseDown:
                    //****************
                    // Mouse-Down -- Mouse button has just been pressed down.
                    // Hit-test and respond if user has clicked a Nuc, DrawHandle or other
                    // interactive element.
                    //****************
                    storeMouse(e); // store the location where the mouse button was first pressed down.
                    if (e.getButton() == MouseEvent.BUTTON1) {
                        IScreenObject hit = hitTest(screenObjects, curPos, false);
                        dragObject = hit==null?null:hit.getModel();
                        if (dragObject == null) {
                            state = SELECTING;
                            selectionBand.setFrameFromDiagonal(startPos, startPos);
                        } else if ((mod & EDIT_CLICK_MASK) == EDIT_CLICK_MASK) {
                            state = EDITING;
                            if (dragObject instanceof Nuc) {
                                focused = (Nuc)dragObject;
                                Point2D p = canvas.getView().toScreen(focused.location, new Point2D.Float());
                                editBond.setLine(p, curPos);
                            }
                        } else if (dragObject instanceof Number) {
                            state = FORMATTING;
                        } else {
                            state = DRAGGING;
                            storeRnaCoords();
                            if (dragObject instanceof Nuc) {
                                /* Note: SELECT_ONLY is used here because we select on mouse-down (so the user
                                   can drag what they've selected), but we only deselect on mouse-up, because
                                   the user may have held down the SELECTION_ADD_MASK (e.g. SHIFT) in order to
                                   alter the method in which the drag alters the drawing
                                   (i.e. SELECTION_ADD_MASK might be the same as DRAG_SPECIAL_MASK)
                                   If we DEselect the target/focused Nuc here, then the rest of the selection will move,
                                   but not the focused Nuc (which is probably not what the user wants).
                                   So: Select on mouse-down; De-select on mouse-up (but only if the use hasn't dragged the focused Nuc)
                                 */
                                boolean append = 0!=(mod & SELECTION_ADD_MASK);
                                deselectNuc = updateSelection(append, (Nuc) dragObject, SELECT_ONLY);
                            } else if (dragObject instanceof DrawHandle) {
                                ((DrawHandle) dragObject).dragStarted(startPos);
                            }
                        }
                    }
                    break;
                case Drag:
                    //****************
                    // Mouse-Drag -- Mouse is dragging over the canvas.
                    // We have already handled the mouse-button-down event and are now in one of a few specific states.
                    //****************
                    updateMousePos(e); // updates prevPos and curPos

                    // There is a minimum number of pixels that the mouse must be dragged in x- or y- directions
                    // before a mouse-drag is considered valid. This is to avoid causing changes if the user is only
                    // clicking on a base etc. and didn't mean to drag it.
                    if (!dragIsValid) {
                        int dx = Math.abs(curPos.x-startPos.x);
                        int dy = Math.abs(curPos.y-startPos.y);
                        dragIsValid = dx>2 || dy>2 || dx+dy>3;
                    }
                    switch(state) {
                        case DRAGGING: // DRAGGING Specifically refers to dragging a Nucleotide or a DrawHandle, which modifies one or more nucleotide locations.
                            if (dragIsValid) {
                                resetRnaCoords();
                                DragOpts options = new DragOpts(0 != (mod & DRAG_EXPANDED_MOTIF_MASK), 0 != (mod & DRAG_SPECIAL_MASK), programSettings().ShowHints);
                                if (dragObject instanceof Nuc)
                                    dragNucs((Nuc) dragObject, startPos, prevPos, curPos, options);
                                else if (dragObject instanceof DrawHandle)
                                    ((DrawHandle) dragObject).drag(startPos, prevPos, curPos, options);
                            }
                            break;
                        case SELECTING: // SELECTING means drawing a selection-box
                            selectionBand.setFrameFromDiagonal(startPos, curPos);
                            redraw(); break;
                        case EDITING: // EDITING means creating or deleting base-pair bonds.
                            editBond.setLine(editBond.x1, editBond.y1, curPos.x, curPos.y);
                            redraw(); break;
                        case FORMATTING: // FORMATTING currently means changing the angle of a base's number label.
                            // In the future FORMATTING could also refer to other style changes, but importantly it
                            // differs from DRAGGING because nucleotide locations are NOT modified, so storing and
                            // resetting coordinates is not required.
                            if (dragIsValid) {
                                if (dragObject instanceof Number)
                                    rotateNumber(scene.getNuc( ((Integer) dragObject) -1), canvas.getView().toModel(curPos));
                                redraw();
                            }
                            break;
                    }
                    break;
                case MouseUp:
                    //****************
                    // Mouse-Up -- Mouse button was released.
                    // This usually signals the end of some drag-event.
                    // Perform completion actions and reset the state of the mouse handler.
                    //****************
                    updateMousePos(e);
                    completeAction(mod);
                    dragObject = null;
                    state = 0;
                    redraw();
                    break;
            }
        }

        private void completeAction(int mod) {
            if (state == SELECTING) {
                boolean add = 0!=(mod & SELECTION_ADD_MASK);
                updateSelection(add, fromPoints(startPos, curPos));
            } else if (state == EDITING) {
                IScreenObject hit = hitTest(screenObjects, curPos, true);
                if (hit != null && hit.getModel() != dragObject)
                    addBond((Nuc)dragObject, (Nuc)hit.getModel());
            } else if (state == FORMATTING) {
                addUndo(SceneUpdateInfo.FormatBases.subType("Rotated Number Label"));
            } else if (state == DRAGGING) {
                // Notify end of modification for undo
                if (dragIsValid && !startPos.equals(curPos)) {
                    if (dragObject instanceof Nuc) {
                        boolean altDrag = 0!=(mod & DRAG_SPECIAL_MASK);
                        // updateLayout has already been called, but in volotile mode (i.e. addHistoryEntry=false)
                        // so now we need to add the history entry
                        addUndo(altDrag?SceneUpdateInfo.Rotate:SceneUpdateInfo.Layout);
                    } else if (dragObject instanceof DrawHandle)
                        ((DrawHandle) dragObject).dragComplete(true);
                } else {
                    if (!dragIsValid && deselectNuc != null && 0!=(mod & SELECTION_ADD_MASK))
                        // i.e. mouse didn't move
                        // Note: Regarding the use of DESELECT_ONLY, please see the note about SELECT_ONLY in the MouseDown case above.
                        updateSelection(true, deselectNuc, DESELECT_ONLY);

                    if (dragObject instanceof DrawHandle) // tell the drag object it is no longer dragging.
                        ((DrawHandle) dragObject).dragComplete(false);
                }
            }
        }
    }

    private void dragNucs(Nuc dragNuc, Point startPos, Point prevPos, Point curPos, final DragOpts options) {
        Collection<Nuc> nucs = selection;
        if (options.expandMotif) { // include the branch
            INucGroup m = dragNuc.getBranch();
            if (m != null)
                nucs = m.getBases();
        }
        if (options.special && nucs.size() > 1)
            rotateNucs(nucs, null, startPos, curPos, new DragOpts(false, false, options.disableHints) ); //curPos.x - prevPos.x, curPos.y - prevPos.y);
        else
            translateScreenNucs(nucs, startPos, curPos); //curPos.x - prevPos.x, curPos.y - prevPos.y);
    }

    private void addBond(final Nuc n1, final Nuc n2) {
        SceneUpdateInfo operation;
        if (n1.getPaired() == n2) {
            n1.getPairBond().unpair();
            operation = SceneUpdateInfo.RemBonds;
        } else {
            Bond b1 = n1.getPairBond(), b2 = n2.getPairBond();
            operation = (b1==null&&b2==null)?SceneUpdateInfo.AddBonds:SceneUpdateInfo.Bonds;
            if (b1 != null) b1.unpair();
            if (b2 != null) b2.unpair();
            scene.addBond(n1, n2);
        }
        structureUpdated(operation);
    }

    private float[] coords = new float[128];
    private void resetRnaCoords() {
        int i = 0;
        for (Strand s : scene.strands)
            for(Nuc n : s) {
                n.location.x = coords[i++];
                n.location.y = coords[i++];
            }
    }
    private void storeRnaCoords() {
        int count = 0;
        for (Strand s : scene.strands)
            count += 2 * s.size();
        if (coords.length < count)
            coords = new float[count + count >> 1];
        int i = 0;
        for (Strand s : scene.strands)
            for(Nuc n : s) {
                coords[i++] = n.location.x;
                coords[i++] = n.location.y;
            }
    }


    private static class NucGroup implements INucGroup {
        public final Collection<Nuc> list;
        public NucGroup(Collection<Nuc> list) { this.list = list; }
        public NucGroup() { list = new ArrayList<>(); }
        @Override public Collection<Nuc> getBases() { return list; }
    }

    private void handleSelectionChange() {
        // Determine which drawHandles should be enabled
        //disableDrawHandles();
        drawHints.clear();

        handleLoopResize.setSelection(focused, selection);
        handleRotate.setSelection(focused, selection);
        handleBranchSlide.setSelection(focused, selection);

//        if (focused == null) {
//            // see if the user has selected a Loop, Helix, Domain, or Branch
//        } else {
//            if (!focused.isPaired())
//
//            if (selection.size() > 1)
//                handleRotate.setSelection(focused, selection);
//            handleBranchSlide.setSelection(focused, getMotif(focused, SelectionType.Branch, true));
//        }
//        int[] sel = getSelectionState();
//        System.out.print("Selection: (F:" + (focused == null ? "-" : focused.indexInScene()) + ") ");
//        for (int i = 0; i < sel.length; i++) {
//            System.out.print(sel[i] + ", ");
//        }
//        System.out.println();
        redraw();
    }
    private void selectMotif(Nuc nuc, boolean sel, SelectionType selectionType) {
        INucGroup ng = getMotif(nuc, selectionType, false);
        selectNucs(ng, sel);
    }
    private void selectNucs(final Iterable<Nuc> nucs, final boolean sel) {
        if (nucs == null) return;
        for (Nuc n : nucs) {
            if (sel)
                selection.add(n);
            else
                selection.remove(n);
        }
    }
    private INucGroup getMotif(Nuc nuc, SelectionType selectionType, boolean allowNull) {
        INucGroup m = null;
        switch (selectionType) {
            case Individual:
                break;
            case Helix:
                m = nuc.getHelix();
                break;
            case Loop:
                m = nuc.getLoop();
                break;
            case HelixOrLoop:
                if (nuc.isPaired())
                    m = nuc.getHelix();
                else
                    m = nuc.getLoop();
                break;
            case Domain:
                m = nuc.getDomain();
                break;
            case Branch:
                m = nuc.getBranch();
                break;
            case DomainOrBranch:
                if (nuc.isPaired())
                    m = nuc.getDomain();
                else
                    m = nuc.getBranch();
                break;
        }
        if (!allowNull && m == null)
            m = nuc;
        return m;
    }
    private static int SELECT_ONLY = 1;
    private static int DESELECT_ONLY = 2;
    /**
     * Updates the set of selected Nucs based on the current selection and the user's intention.
     * If the user has indicated they want to modify (aka append) the existing selection,
     * then the target nuc's selection status is flipped (i.e. if already selected, it should be deselected and vice-versa)
     * Otherwise, the target Nuc should be selected.
     *
     * If the target Nuc is to be DE-selected, but the allowMode is set to SELECT_ONLY, then the
     * target Nuc (which should be deselected at some point in the future, but not immediately)
     * will be returned by the function.
     *
     * @param append Whether to modify the selection (true) or clear the existing selection (false).
     * @param nuc The target Nuc.
     * @param allowMode Whether we are in SELECT_ONLY mode or DESELECT_ONLY.
     * @return The target Nuc if it is to be deselected, but we are in SELECT_ONLY mode. Otherwise Null.
     */
    private Nuc updateSelection(final boolean append, @NotNull final Nuc nuc, int allowMode) {
        boolean sel = selection.contains(nuc);
        //-debug: System.out.println("Sel: " + sel + " Append: " + append);
        if (sel && append) {
            //-debug: System.out.println("1");
            // Already selected, but should be deselected.
            if (allowMode == SELECT_ONLY) return nuc;
            //-debug: System.out.println("1");
            selectMotif(nuc, false, selType);
            selection.remove(nuc);
            if (focused == nuc) focused = null;
        } else if (!sel) {
            //-debug: System.out.println("2");
            if (allowMode == DESELECT_ONLY) return null;
            //-debug: System.out.println("2");
            if (!append)
                clearSelection();
            selectMotif(nuc, true, selType);
            focused = nuc;
        } else { // sel && !append
            //-debug: System.out.println("3");
            if (allowMode == DESELECT_ONLY) return null;
            //-debug: System.out.println("3");
            // Do not change selection if re-selecting without append.
            // I.e. assume user is just taking a break from pressing the mouse button or is re-adjusting focus.
            focused = nuc;
        }
        //-debug: System.out.println("selchanged");
        selectionUpdated();
        return null;
    }

    private void updateSelection(final boolean append, final Rectangle selRect) {
        Collection<IScreenObject> objects = hitTest(screenObjects, selRect, true);
        if (!append && (objects==null || objects.size()==0)) {
            clearSelection();
        } else  {
            List<Nuc> list = new ArrayList<>(objects.size());
            int selected = 0;
            boolean select;
            for (IScreenObject obj: objects)
                if (obj.getModel() instanceof Nuc)
                    list.add((Nuc)obj.getModel());

            if (!append) { // not appending to selection, so just clear selection and add those inside rect.
                selection.clear();
                select = true;
            } else {
                // we are appending, so select or deselect based on the state of the majority of nucs inside rect.
                // determine majority state
                for (Nuc n : list)
                    if (selection.contains(n))
                        selected++;
                select = selected < list.size() / 2;
            }
            selectNucs(list, select);
            if (select  && list.size() != 0 && (focused == null || !selection.contains(focused)))
                focused = list.get(0);
        }
        selectionUpdated();
    }

    private IScreenObject hitTest(Iterable<IScreenObject> list, final Point2D pt, boolean nucsOnly) {
        IScreenObject found = null;
        int order = IScreenObject.ZORDER_MIN;
        if (list == null) return null;
        for (IScreenObject obj : list)
            if ((!nucsOnly || obj.getModel() instanceof Nuc) && obj.getZOrder() > order && obj.hitTest(pt)) {
                if (obj.getZOrder() == IScreenObject.ZORDER_MAX) return obj;
                found = obj;
                order = found.getZOrder();
            }
        return found;
    }
    private List<IScreenObject> hitTest(Iterable<IScreenObject> list, final Rectangle rc, boolean nucsOnly) {
        List<IScreenObject> found = new ArrayList<>();
        if (list == null) return found;
        for (IScreenObject obj : list)
            if ((!nucsOnly || obj.getModel() instanceof Nuc) && obj.hitTest(rc))
                found.add(obj);
        return found;
    }
        // _pt.x = x; _pt.y = y;
//        try {
//            canvasComponent.getViewTransform().inverseTransform(_pt, _pt);
//        } catch (NoninvertibleTransformException ex) {
//            ex.printStackTrace();;
//        }
        //IDrawable found = hitTestHandles(_pt);
        //if (found != null)
        //    return found;

        //((ICanvas)canvasComponent).pointToModel(_pt);
        //return hitTestNuc(_pt);
//    }
//    private IDrawable hitTestHandles(final Point2D.Float pt) {
//        return null;
//    }
//    private IDrawable hitTestNuc(Point2D.Float pt) {
//        for (Strand s : scene.strands)
//        for (int i = s.size()-1; i > -1; i--) {
//            if (s.get(i).hitTest(pt, settings))
//                return s.get(i);
//        }
//        return null;
//    }
//    private IDrawable hitTestNuc(Rectangle2D.Float rc) {
//        for (Strand s : scene.strands)
//            for (int i = s.size()-1; i > -1; i--) {
//                if (s.get(i).hitTest(pt, settings))
//                    return s.get(i);
//            }
//        return null;
//    }
    public List<IScreenObject> draw(Graphics2D g, View2D view, DrawSettings settings, boolean interactive) {
        AffineTransform prevTransform = g.getTransform(); // store so it can be reset later.

        int objCount = scene.getNucCount() * 2; // include backbone segments
        //objCount += bonds.size(); // include bonds.
        List<IScreenObject> objects = new ArrayList<>(objCount);

        g.setBackground(settings.sceneBgColor);
        //Rectangle bounds = view.screenBounds;
        Rectangle clip = g.getClipBounds();
//        if (clip == null || clip.isEmpty())
//            clip = bounds;
        if (clip != null && !clip.isEmpty())
            g.clearRect(clip.x, clip.y, clip.width, clip.height);

        g.transform(view.trToScreen);

        // Debug transform
//        Point p = new Point(100,100);
//        Point2D pT = prevTransform.transform(p, null);
//        Point2D pV = view.transform(p, null);
//        Point2D pG = g.getTransform().transform(p, null);
//        System.out.println(fmt("Transform #%s:  T=(%s,%s)  V=(%s,%s)   G=(%s,%s)", ++drawCount, pT.getX(), pT.getY(), pV.getX(), pV.getY(), pG.getX(), pG.getY()));

        g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);

        // Draw the Basepairs (aka Bonds)
        drawBonds(g, view, objects);

        // Draw the backbone segments
        if (settings.backboneColor != null && (settings.drawBackboneInHelix | settings.drawBackboneInLoops)) {
            g.setColor(settings.backboneColor);
            g.setStroke(settings.backboneStroke);
            float minBackboneDist = 2 * (settings.nucRadius + settings.nucOutlineWidth);
            boolean[] backboneDone = new boolean[scene.getNucCount()];
            if (settings.drawBackboneCurved && settings.drawBackboneInLoops)
                drawCurvedBackbones(g, backboneDone, minBackboneDist);
            drawSimpleBackbones(g, backboneDone, minBackboneDist);
        }

        // Draw the nucleotides
        for (Strand strand : scene.strands) {
            int order = 0;
            NucStyle style = new NucStyle();
            for (Nuc n : strand) {
                settings.fillStyle(n, style);

//                if (n.isPaired()) {
//                    Nuc n2 = n.getPaired();
//                    g.setColor(n.compareTo(n3)<0?Color.CYAN:Color.YELLOW);
//                    g.setStroke(n.compareTo(n3)<0?new BasicStroke(3):new BasicStroke(1));
//                    g.drawLine((int) n.location.x, (int) n.location.y, (int) n3.location.x, (int) n3.location.y);
//                }

                Shape s = getNucShape(n);
                objects.add(getScreenObject(n, view, ++order));

                if (style.fillColor != null) {
                    g.setColor(style.fillColor);
                    g.fill(s);
                }

                if (style.textColor != null) {
                    g.setFont(style.font);
                    g.setColor(style.textColor);
                    GraphicsUtil.drawTextCentered(g, n.symbol, n.location); // + (n.indexInStrand()+1)
                }

                if (style.outlineColor != null) {
                    g.setColor(style.outlineColor);
                    g.setStroke(settings.nucOutlineStroke);
                    g.draw(s);
                }
            }
        }

        if (interactive) {
            Nuc bonding = null;
            if (mouseHandler.state == MouseHandler.EDITING) {
                IScreenObject so = hitTest(objects, editBond.getP2(), true);
                if (so != null) bonding = (Nuc)so.getModel();
            }

            if (mouseHandler.state == MouseHandler.SELECTING) {
                for (IScreenObject o : hitTest(objects, selectionBand, true))
                    if (o.getModel() != focused && o.getModel() != bonding)
                        drawNucOutline(g, (Nuc)o.getModel(), settings.nucSelColor);
            }
            for (Nuc n : selection)
                if (n != focused && n != bonding)
                    drawNucOutline(g, n, settings.nucSelColor);
            if (bonding != null && bonding != focused)
                drawNucOutline(g, bonding, settings.nucHitColor);
            if (focused != null)
                drawNucOutline(g, focused, settings.nucFocusColor);

            g.setTransform(prevTransform);
            for (DrawHint dh : drawHints)
                dh.draw(g, view, settings);

            drawHints.clear();

            g.setTransform(prevTransform); // i.e. draw in screen coordinates instead of model coords.

            //DrawHandle dragHandle = mouseHandler.dragObject instanceof DrawHandle ? (DrawHandle)mouseHandler.dragObject : null;

            for (DrawHandle obj : drawHandles)
                //if (obj != null && obj.isEnabled() && obj.isValid()) {
                if (obj != null && (mouseHandler.dragObject == null || mouseHandler.dragObject == obj) && obj.isEnabled() && obj.isValid()) {
                    objects.add(new ScreenShape(obj, obj.draw(g, view, settings, mouseHandler.dragObject == obj), IScreenObject.ZORDER_MAX));
//                    g.setColor(Color.RED);
//                    g.draw(objects.get(objects.size() - 1).shape);
                }

            if (mouseHandler.state == MouseHandler.SELECTING) {
                g.setStroke(settings.selBandOutline);
                g.setColor(settings.selBandLineColor);
                g.draw(selectionBand);
                g.setColor(settings.selBandFillColor);
                g.fill(selectionBand);
            } else if (mouseHandler.state == MouseHandler.EDITING) {
                g.setStroke(settings.editBondStroke);
                g.setColor(settings.bondColor);
                g.draw(editBond);
            }
        } // if (interactive)


        g.setTransform(prevTransform);
        if (settings.drawNumbers) {
            // draw nucleotide numbers
            g.transform(view.trToScreen);
            g.setColor(settings.numberLineColor);
            g.setStroke(settings.numberLineStroke);
            g.setFont(settings.numberFont);

            final int spacing = settings.drawNumbersInterval;
            for (int i = 0; i < scene.getNucCount(); i++) {
                if (i == 0 || ((i + 1) % spacing) == 0)
                    if (shouldDrawNumber(scene.getNuc(i)))
                        objects.add(drawNucNumber(g, view, scene.getNuc(i), i + 1));
            }
            int last = scene.getNucCount() - 1;
            if (((last + 1) % spacing > (spacing - 1) / 2))
                if (shouldDrawNumber(scene.getNuc(last)))
                    objects.add(drawNucNumber(g, view, scene.getNuc(last), last + 1));
        }

        return objects;
    }
    private boolean shouldDrawNumber(Nuc n) {
        return n.style==null || !Float.isInfinite(n.style.numberOffset);
    }

    private void drawBonds(Graphics2D g, View2D view, List<IScreenObject> objects) {
        List<Bond> bonds = scene.getBonds();
        Rectangle2D bondBounds = null; // track the bounds of the bonds. Linear draw mode draws bonds outside the bounds of the nucs, so we need to add a IScreenObject so that the total bounds can be calculated.
        Line2D.Float bondLine = new Line2D.Float();
        for (Bond b : bonds) {
            switch (b.getType()) {
                case Pseudo:
                    g.setColor(settings.bondColorPseudo);
                    g.setStroke(settings.bondStrokePseudo);
                    break;
                case Prohibited:
                case Forced:
                case Standard:
                case Special:
                default:
                    g.setColor(getBondColor(b));
                    g.setStroke(settings.bondStroke);
                    break;
            }
            Shape bondShape;
            switch (scene.drawMode) {
                case Circular:
                    bondShape = NucLayout.getCircularBond(b);
                    //addHint(Ellipses.fromCenter(((QuadCurve2D.Float)bondShape).getCtrlPt(), settings.nucRadius / 2)).fromModel().fill(Colors.hot);
                    break;
//                case Flattened:
//                    bondShape = null;
//                    break;
                case Linear:
                    bondShape = NucLayout.getLinearBond(b);
                    break;
                case Standard:
                    bondLine.setLine(b.n5.location, b.n3.location);
                    bondShape = bondLine;
                    break;
                default:
                    bondShape = null;
                    break;
            }
            if (bondShape == null) continue;
            if (bondBounds == null)
                bondBounds = bondShape.getBounds2D();
            else
                Rectangle2D.union(bondBounds, bondShape.getBounds2D(), bondBounds);
            g.draw(bondShape);
            //g.drawLine((int) n1.location.x, (int) n1.location.y, (int) n2.location.x, (int) n2.location.y);
        }
        if (bondBounds != null) {
            //g.setColor(Color.RED);
            //g.draw(bondBounds);
            objects.add(new ScreenShape(bonds, view.trToScreen.createTransformedShape(bondBounds), IScreenObject.ZORDER_MIN));
        }
    }

    private Vec2D _backboneDir = new Vec2D(), _backboneStart = new Vec2D(), _backboneEnd = new Vec2D();
    Line2D.Float _backboneLine = new Line2D.Float();
    private void drawSimpleBackbones(Graphics2D g, final boolean[] done, final float minDist) {
        for (Strand strand : scene.strands) {
            if (strand.size() > 1) {
                for (int i = 1; i < strand.size(); i++) {
                    Nuc prev = strand.get(i - 1);
                    Nuc n = strand.get(i);
                    if (done[n.indexInScene()]) continue;
                    if (n.isPaired(false) && prev.isPaired(false)) {
                        // In helix
                        if (!settings.drawBackboneInHelix) continue;
                    } else {
                        // In loop OR connecting loop to helix.
                        if (!settings.drawBackboneInLoops) continue;
                    }
                    _backboneDir.setTo(n.location).subtr(prev.location);
                    float length = (float) _backboneDir.length();
                    if (length < minDist)
                        continue;
                    _backboneDir.setLength(1);
                    _backboneStart.setTo(prev.location).addScaled(_backboneDir, settings.nucRadius + settings.nucOutlineWidth);
                    _backboneEnd.setTo(prev.location).addScaled(_backboneDir, length - (settings.nucRadius + settings.nucOutlineWidth));
                    _backboneLine.setLine(_backboneStart, _backboneEnd);
                    //g.drawLine((int) backboneStart.x, (int) backboneStart.y, (int) backboneEnd.x, (int) backboneEnd.y);
                    g.draw(_backboneLine);
                }
            }
        }
    }
    // TODO: Implement curved bond drawing.
    private void drawCurvedBackbones(Graphics2D g, final boolean[] done, final float minDist) {
//        NucLayout layout = new NucLayout();
//        for(Motif.MultiLoop m : scene.getMultiLoops(0, true)) {
//            Ellipses.Circle fullCircle = layout.calcMultiLoopBestFitCircle(m);
//            List<Bond> bonds = m.getBonds();
//            if (bonds.size() == 1) {
//                // hairpin
//                Bond b = bonds.get(0);
//                if (allNucsFitCircle(b.n5, b.n3, fullCircle))
//                    drawCircleBond(g, b.n5, b.n3, fullCircle);
//            } else {
//                for (int i = 0; i < bonds.size(); i++) {
//                    Bond b = bonds.get(i);
//                    Bond next = i + 1 < bonds.size() ? bonds.get(i + 1) : bonds.get(0);
//                }
//            }
//        }
    }
    private boolean allNucsFitCircle(Nuc start, Nuc end, Ellipses.Circle c) {
        for(int i = start.indexInScene(); i <= end.indexInScene(); i++)
            if (!nucFitsCircle(i, c)) return false;
        return true;
    }
    private boolean nucFitsCircle(int sceneIndex, Ellipses.Circle c) {
        return Math.abs(scene.getNuc(sceneIndex).location.distance(c.center()) - c.radius()) > settings.nucRadius / 2;
    }


    private Color getBondColor(Bond b) {
        if (b.n5.style != null && b.n5.style.bondColor != null)
            return b.n5.style.bondColor;
        if (b.n3.style != null && b.n3.style.bondColor != null)
            return b.n3.style.bondColor;
        return settings.bondColor;
    }

    public Vec2D calcNormalToNuc(Nuc cur) {
        Vec2D dir;
        if (cur.isPaired() && scene.drawMode == SceneDrawMode.Standard) {
            dir = new Vec2D(cur.getPaired().location, cur.location);
        } else {
            Nuc prev = cur.getPrev(), next = cur.getNext();
            if (next == null && prev == null)
                dir = new Vec2D(1, 0);
            else if (next == null || prev == null) {
                // This is a terminal nuc, so return a vector that is perpendicular to the ray from this nuc to the previous/next one.
                // It's arbitrary whether to rotate by 90 or -90, so take a guess, depending on draw direction and flip status.
                boolean isEnd = next==null;
                prev = isEnd ? prev : next;
                dir = new Vec2D(prev.location, cur.location).rotate90();
                if (isEnd == scene.drawFlipped)
                    dir.negate();
            } else {
                Vec2D vp = new Vec2D(prev.location, cur.location).normalize(1,0),
                        vn = new Vec2D(next.location, cur.location).normalize(1,0);
                dir = Vec2D.getMidpoint(vp, vn);
                //double a1 = Math.atan2(prev.location.y-cur.location.y, prev.location.x-cur.location.x);
                //double a2 = Math.atan2(next.location.y-cur.location.y, next.location.x-cur.location.x);
                //double avg = (a1+a2)/2;
                //double diff = Math.abs((a1<0?a1+Math.PI*2:a1) - (a2<0?a2+Math.PI*2:a2));

//                System.out.printf("Aprev:%.0f  Anext:%.0f  sum:%.0f  avg:%.0f  diff:%.0f angle:%.0f\n",
//                        a1*180/Math.PI,
//                        a2*180/Math.PI,
//                        (a1+a2)*180/Math.PI,
//                        avg*180/Math.PI,
//                        diff*180/Math.PI,
//                        (Math.PI+avg)*180/Math.PI
//                        );
                if (dir.length() < 0.0000001)
                    dir = new Vec2D(prev.location, cur.location).rotate90();

//                if (Math.abs(diff-Math.PI) < 0.001)
//                    dir = new Vec2D(prev.location, cur.location).rotate90();
//                else
//                    dir = Vec2D.fromAngle(Math.PI + avg);

                //Vec2D mid = Vec2D.getMidpoint(prev.location, next.location);
                //dir = new Vec2D(mid, cur.location);
            }
        }
        if (dir.length() == 0) dir = new Vec2D(1, 0);
        return dir.normalize();
    }

    private final int ZORDER_NUMBER = -1;
    private IScreenObject drawNucNumber(final Graphics2D g, View2D view, final Nuc cur, final int number) {
        float r = settings.nucRadius;
        Vec2D dir = calcNormalToNuc(cur);

        if (cur.style!=null && cur.style.numberOffset!=0)
                dir.rotate(cur.style.numberOffset);

        float x=cur.location.x, y = cur.location.y, dx=(float)dir.x, dy=(float)dir.y;
        float x1=x+dx*(r+settings.numberLineMargin), y1=y+dy*(r+settings.numberLineMargin), x2, y2; //x2=x+dx*(numberDistance-numberLineMargin), y2=y+dy*(numberDistance-numberLineMargin);

        Rectangle2D.Float bounds = GraphicsUtil.drawTextCentered(g, Integer.toString(number), x+dx*(r+settings.numberDistance), y+dy*(r+settings.numberDistance));

        float MARGIN = (float)Math.sqrt(r+settings.numberLineMargin);
        x2 = dx > 0 ? bounds.x - MARGIN : bounds.x + bounds.width + MARGIN;
        y2 = dy > 0 ? bounds.y - MARGIN : bounds.y + bounds.height + MARGIN;

        // Calculate the maxim length of the line that keeps it outside of the text box
        float L;
        if (Math.abs(dx) < 0.00001)
            L = (y2-y)/dy;
        else if (Math.abs(dy) < 0.00001)
            L = (x2-x)/dx;
        else
            L = Math.max((x2-x)/dx, (y2-y)/dy);

        // recalculate x2 and y2 from L
        x2=x+dx*L; y2=y+dy*L;

        g.drawLine((int)x1, (int)y1, (int)x2,(int) y2);

        // Create a shape that roughly bounds the number
        // Shape bounds = Ellipses.fromCenter(text, nucRadius());
        // Return a IScreenObject (of the bounds transformed for the screen)
        return new ScreenShape(number, view.trToScreen.createTransformedShape(bounds), ZORDER_NUMBER);
    }
    //private void drawLine(Graphics2D g, float x1, float y1, float x2, float y2) { g.drawLine((int)x1, (int)y1, (int)x2, (int)y2); }

    Ellipse2D.Float tmpNucShape = new Ellipse2D.Float();
    private void drawNucOutline(Graphics2D g, Nuc n, Color c) {
        Ellipse2D.Float shp = (Ellipse2D.Float)nucShapes[n.indexInScene()];
        if (shp == null)
            return;
        float margin = 0;
        tmpNucShape.setFrame(shp.x - margin, shp.y - margin, shp.width + 2*margin, shp.height + 2*margin);
        g.setStroke(settings.nucOutlineSelStroke);
        g.setColor(c);
        g.draw(tmpNucShape);
    }

    public abstract static class Colors {
        public static Color hot = new Color(0xFF004E);
        public static Color warm = new Color(0xFF9900);
        public static Color cool = new Color(0x23B3FF);
        public static Color cold = new Color(0x28519D);
        public static Color earth = new Color(0x663713);
        public static Color sand = new Color(0xAAA621);
        public static Color grass = new Color(0x084208);
    }

    private Shape[] nucShapes;
    private ScreenNuc[] nucScreenShapes;
    private Shape getNucShape(Nuc n) {
        Ellipse2D.Float s = (Ellipse2D.Float)nucShapes[n.indexInScene()];
        float r = settings.nucleotideRadius(n);
        if (s == null)
            nucShapes[n.indexInScene()] = s = new Ellipse2D.Float();
        s.x = n.location.x-r;
        s.y = n.location.y-r;
        s.width = s.height = r*2;
        return s;
    }


    private static class ScreenNuc  implements IScreenObject {
        public final int selectionBorder = 3; // user can select from a few px away. In screen pixels. (Not affected by view scale)
        final float OUTLINE_WIDTH = 0.5f; // width of outline stroke (which gets magnified by the view's transform scale)
        public final Nuc nuc;
        public int zorder;
        public DrawSettings settings;
        public final Point p = new Point();
        float screenRadius;
        int hitDistance, hitDistanceSq;
        Rectangle bounds = new Rectangle();
        public ScreenNuc(final Nuc nuc, final DrawSettings settings) {
            this.nuc = nuc;
            this.settings = settings;
        }
        public void setZorder(final int newZorder) { this.zorder = newZorder; }
        @Override
        public void updateView(final View2D view) {
            view.trToScreen.transform(nuc.location, p);
            screenRadius = (float) ((settings.nucRadius+((BasicStroke)settings.nucOutlineSelStroke).getLineWidth()/2) * view.trToScreen.getScaleX());
            hitDistanceSq = (int)Math.ceil((screenRadius+selectionBorder)*(screenRadius+selectionBorder));
            hitDistance = (int)Math.ceil(screenRadius+selectionBorder);
            bounds.setBounds(Math.round(p.x - screenRadius), Math.round(p.y - screenRadius), Math.round(2*screenRadius), Math.round(2*screenRadius));
        }
        @Override
        public boolean hitTest(final Point2D point) {
            return Vec2D.distanceSq(p.x, p.y, point.getX(), point.getY()) < hitDistanceSq;
        }
        @Override
        public boolean hitTest(final Rectangle rc) {
            int dx = Math.abs(p.x - (rc.x+rc.width/2));
            int dy = Math.abs(p.y - (rc.y+rc.height/2));

//            if (nuc.indexInScene()==10)
//                System.out.printf("p:%s  rc:%s  dx:%s dy%s rcw/2:%s, hitDist:%s\n", p, rc, dx, dy, rc.width/2, hitDistance);

            // If the distance between the center of the rectangle and the circle is larger than
            // the radius (in x or y direction) then they cannot possibly be intersecting.
            if (dx > (rc.width/2 + hitDistance)) return false;
            if (dy > (rc.height/2 + hitDistance)) return false;

            // If the distance between the center of the rectangle and the circle is smaller than
            // half of the rectangle size (in x or y direction) then they MUST be intersecting
            // (because the circle center is INSIDE the rectangle)
            if (dx <= (rc.width/2) || dy <= (rc.height/2)) return true;

            // The circle and rectangle are close enough to be possibly touching,
            // but only at the corner of the rectangle.
            // Determine if the distance from the center of the circle to the corner is less than the radius.
            double cornerDistSq = (dx - rc.width/2)*(dx - rc.width/2) + (dy - rc.height/2)*(dy - rc.height/2);
            return cornerDistSq <= hitDistanceSq;
        }
        @Override
        public Nuc getModel() {
            return nuc;
        }
        @Override
        public Rectangle getBounds() {
            return bounds;
        }
        @Override
        public int getZOrder() {
            return zorder;
        }
    }

    private IScreenObject getScreenObject(Nuc n, View2D view, int zorder) {
        ScreenNuc sn = nucScreenShapes[n.indexInScene()];
        if (sn == null)
            sn = nucScreenShapes[n.indexInScene()] = new ScreenNuc(n, settings);
        sn.updateView(view);
        sn.setZorder(zorder);
        return sn;
    }

    public boolean addSceneUpdateListener(Consumer<SceneUpdateEvent> listener) {        return sceneUpdated.add(listener);    }
    public boolean removeSceneUpdateListener(Consumer<SceneUpdateEvent> listener) {        return sceneUpdated.remove(listener);    }
    public boolean addHistoryListener(Consumer<HistoryUpdateEvent> listener) {        return historyEvent.add(listener);    }
    public boolean removeHistoryListener(Consumer<HistoryUpdateEvent> listener) {        return historyEvent.remove(listener);    }

    private int _suspendUpdatesCounter = 0;
    public void suspendUpdates() {
        _suspendUpdatesCounter++;
    }
    public void resumeUpdates() {
        if (_suspendUpdatesCounter==0)
            throw new IllegalStateException("resumeUpdates called too many times");
        _suspendUpdatesCounter--;
    }
    private int _suspendRedraw = 0;
    public void suspendRedraw() {
        _suspendRedraw++;
    }
    public void resumeRedraw(boolean forceRedraw) {
        if (_suspendRedraw==0)
            throw new IllegalStateException("resumeRedraw called too many times");
        _suspendRedraw--;
        if (forceRedraw && _suspendRedraw==0) redraw();
    }


    @Override public void selectionUpdated() {  notifyUpdate(SceneUpdateInfo.Selection, false); }
    @Override public void viewUpdated() { notifyUpdate(SceneUpdateInfo.View, false); }
    @Override public void controlsUpdated() { notifyUpdate(SceneUpdateInfo.Controls, false); }
    @Override public void historyUpdated() { notifyUpdate(SceneUpdateInfo.EditHistory, false); }
    //@Override public void styleUpdated() { notifyUpdate(SceneUpdateInfo.FormatBases, true); }

    @Override public void structureUpdated(SceneUpdateInfo details) { notifyUpdate(details); }
    @Override public void styleUpdated(SceneUpdateInfo details) { notifyUpdate(details); } // note: addHistoryEntry is default:true
    @Override public void layoutUpdated(SceneUpdateInfo details) { layoutUpdated(details, true); }
    @Override public void layoutUpdated(SceneUpdateInfo details, boolean addHistoryEntry) {
        suspendUpdates();
        try {
            for (DrawHandle h : drawHandles)
                if (h != null && h.isEnabled() && h.isValid())
                    h.layoutUpdated();
        } finally {
            resumeUpdates();
        }
        notifyUpdate(details, addHistoryEntry);
    }
    //@Override public void notifyUpdate(SceneUpdateCategory update) { notifyUpdate(update, null); }
    @Override public void notifyUpdate(SceneUpdateInfo update, boolean addHistoryEntry) {
//        if (update.category == SceneUpdateCategory.Structure)
//            scene.identifyPseudoknots();
        if (_suspendUpdatesCounter==0) {
            if (addHistoryEntry) addUndo(update);
            SceneUpdateEvent e = new SceneUpdateEvent(this, update, scene);
            sceneUpdated.invoke(e);
        }
    }
    @Override public void notifyUpdate(SceneUpdateInfo details) {
        notifyUpdate(details, true);
    }

    public void clearUndo() {
        scene.history.clear();
        addUndo(SceneUpdateInfo.Initial);
    }
    public void addUndo(final SceneUpdateInfo type) {
        HistoryState state = new HistoryState(scene.getState(), type);
        state.put("selection", getSelectionState());
        //state.put("focused", focused == null ? -1 : focused.indexInScene());
        scene.history.store(state);
        historyEvent.invoke(new HistoryUpdateEvent(this, true, state));
        historyUpdated();
    }
    public boolean canUndo() { return scene.history.canUndo(); }
    public boolean canRedo() { return scene.history.canRedo(); }
    public void undo() {
        if (scene.history.canUndo()) {
            applyHistoryState(scene.history.undo());
        }
    }
    public void redo() {
        if (scene.history.canRedo()) {
            applyHistoryState(scene.history.redo());
            //HistoryState u = scene.history.redo();
            //scene.setState(u.state);
        }
    }

    private void applyHistoryState(HistoryState state) {
        suspendRedraw();
        try {
            scene.loadState(state.state);
            historyEvent.invoke(new HistoryUpdateEvent(this, false, state));
            notifyUpdate(SceneUpdateInfo.GotoHistory, false);
            historyUpdated();
            setSelectionState(state.get("selection", EMPTY_SELECTION));
        } finally {
            resumeRedraw(true);
        }
//        int i = state.get("focused", -1);
//        focused = i == -1 ? null : scene.nucAtSceneIndex(i);
    }

    public void setSelType(SelectionType sel) {
        this.selType = sel;
    }

    /**
     * If a nucleotide is selected, this expands to the helix or loop.
     * If a helix or loop is selected, this expands to the containing branch.
     */
    public void expandSelection() {
        // For each nucleotide in the selection, determine if its enclosing loop or branch is selected
        boolean marked[] = new boolean[scene.getNucCount()]; // cache the results of each nucleotide in the branch, so we don't add the same branch over and over.
        boolean useBranch = isFullMotifSelected();
        // select each branch (if useBranch is true) or loop/helix
        for (Nuc n : selection.toArray(new Nuc[selection.size()]))
            if (!marked[n.indexInScene()]) {
                IMotif m = useBranch ? n.getBranch() : n.isPaired() ? n.getHelix() : n.getLoop();
                if (m == null) continue;
                Collection<Nuc> bases = m.getBases();
                for (Nuc mn : bases)
                    marked[mn.indexInScene()] = true;
                selection.addAll(bases);
            }
        selectionUpdated();
    }

    // determines whether the entire loop or helix is selected for each nucleotide in the selection.
    private boolean isFullMotifSelected() {
        boolean marked[] = new boolean[scene.getNucCount()]; // cache the results of each nucleotide in the motifs, so we don't test the same loop/helix multiple times.
        for (Nuc n : selection)
            if (!marked[n.indexInScene()]) {
                IMotif m = n.isPaired() ? n.getHelix() : n.getLoop();
                if (m != null)
                    for (Nuc mn : m.getBases()) {
                        if (selection.contains(mn))
                            marked[mn.indexInScene()] = true;
                        else
                            return false;
                    }
            }
        return true;
    }

    public void selectAll() {
        selection.clear();
        for(Strand s : scene.strands)
            selection.addAll(s);
        selectionUpdated();
    }
    public enum SelectionType {
        Individual,
        Helix,
        Loop,
        HelixOrLoop,
        Domain,
        Branch,
        DomainOrBranch
    }
}
