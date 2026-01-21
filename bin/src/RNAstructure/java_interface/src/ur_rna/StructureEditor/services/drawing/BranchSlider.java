package ur_rna.StructureEditor.services.drawing;

import ur_rna.StructureEditor.Program;
import ur_rna.StructureEditor.models.Bond;
import ur_rna.StructureEditor.models.Motif;
import ur_rna.StructureEditor.models.Nuc;
import ur_rna.StructureEditor.models.SceneUpdateInfo;
import ur_rna.StructureEditor.services.SceneController;
import ur_rna.Utilities.geom.Ellipses.Circle;
import ur_rna.Utilities.geom.Vec2D;

import java.awt.*;
import java.awt.geom.Point2D;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import static ur_rna.Utilities.geom.PointMath.magnitude;
import static ur_rna.Utilities.geom.PointMath.scale;

/**
 * A DrawHandle that allows a user to "slide" an RNA branch (i.e. a helix and everything on one side of it.) along the smooth arc of the circular loop at the helix's base.
 */
public class BranchSlider extends DrawHandle {
    private Motif.Helix helix;
    private Bond base, b5, b3;
    private Motif.Loop l5, l3;
    private Set<Nuc> branchNucs = new HashSet<>();

    private Point handleOffset = new Point();
    public BranchSlider(final SceneController controller) {
        super(controller);
        boxSize = 20;
        setIcon(Program.getImage("slide-branch-handle"));
    }

    @Override
    public boolean isValid() {
        return helix != null;
    }
    @Override
    public void setSelection(final Nuc focus, final Collection<Nuc> selection) {
        super.setSelection(focus, selection);
    }

    @Override
    protected void nucPositionsChanged() {
        // do nothing. handled by prepareDraw
    }
    private Point2D.Float _ptModel = new Point2D.Float();

    /** Called when this handle has been dragged by the user. */
    @Override
    protected void performDrag(Point start, Point prev, Point current, final SceneController.DragOpts options) {
        moveLocation(start, current);
        getModelSpaceDeltaTransform(start, current).transform(targetNuc.location, _ptModel);
        positionNucs(_ptModel, options);
        //targetNuc.transform(getModelSpaceDeltaTransform(start, current));
    }
    private void positionNucs(final Point2D.Float newTargetLoc, final SceneController.DragOpts options) {
        if (base == null) return;
        Vec2D pb, p5, p3;
        pb = Vec2D.getMidpoint(base.n5.location, base.n3.location);
        p5 = b5 == null ? null : Vec2D.getMidpoint(b5.n5.location, b5.n3.location);
        p3 = b3 == null ? null : Vec2D.getMidpoint(b3.n5.location, b3.n3.location);

        if (p5 == null && l5 != null && l5.getType() == Motif.Loop.LoopType.Terminal)
            p5 = new Vec2D(l5.getNuc5().location);
        if (p3 == null && l3 != null && l3.getType() == Motif.Loop.LoopType.Terminal)
            p3 = new Vec2D(l3.getNuc3().location);

        Circle circ = null;
        if (p5 != null && p3 != null) {
            if (p5.equals(p3) || p5.distanceSq(p3) < 5) {
                // use a loop nuc if possible
                if (l3 != null)
                    circ = Circle.calcFromPoints(p3, pb, l3.getNucAt(l3.size()/2).location);
                else if (l5 != null)
                    circ = Circle.calcFromPoints(p3, pb, l5.getNucAt(l5.size()/2).location);
                else // Note that if l3 is null, b3 cannot be. same with l5 and b5 (because otherwise p5 or p3 would be null)
                    circ = Circle.calcFromPoints(pb, b3.n3.location, b3.n5.location);
            } else {
                circ = Circle.calcFromPoints(p3, p5, pb);
            }
        }
        if (circ == null) {
            positionNucsRotate(newTargetLoc, options);
            return;
        }
        if (!options.disableHints)
            controller.addHint(circ).fromModel();

        float rad = controller.getSettings().nucleotideRadius(targetNuc);

        double distBase = Math.min(Vec2D.distanceSq(circ.center(), base.n5.location), Vec2D.distanceSq(circ.center(), base.n3.location)); // distance between center and nearest base nuc
        double distTarget = Vec2D.distanceSq(circ.center(), targetNuc.location);  // Distance of target nucleotide before user's operation

        // determine whether the helix is "flipped" inside the loop
        boolean flippedInside = distTarget < distBase-rad*rad;

        // determine whether the drawing is "flipped" (drawn counter-clockwise) the helix
        Vec2D centerToBase = new Vec2D(circ.center(), pb);
        Vec2D normalAtBase = new Vec2D(base.n5.location, base.n3.location).rotate90();

        boolean flipped5to3 = centerToBase.dot(normalAtBase) < 0;
        if (!options.disableHints)
            controller.addHint(normalAtBase.toLine(circ.center(), centerToBase.length() * 1.25)).color(flipped5to3 ? Color.RED : Color.ORANGE).fromModel();

        if (flippedInside != flipped5to3)
            normalAtBase.negate();
        double angle = normalAtBase.angleTo(centerToBase);
        // first orient the branch to face the circle.
        if (Math.abs(angle) > Math.PI / 90)
            controller.rotateNucs(branchNucs, pb, angle);

        // Now rotate around the multi-branch or internal loop circle.
        controller.rotateNucs(branchNucs, circ.center(), startLocation, location, options, helix, false);

        controller.addHint(new Circle(base.n5.location, 4)).fromModel().color(Color.RED);
        controller.addHint(new Circle(base.n3.location, 4)).fromModel().color(Color.ORANGE);

//        if (b5!=null) {
//            controller.addHint(new Circle(b3.midpoint(), 8)).fromModel().color(Color.GREEN);
//            controller.addHint(new Circle(b5.left(base.n3).location, 12)).fromModel().color(Color.BLUE);
//            controller.addHint(new Circle(b5.right(base.n3).location, 12)).fromModel().color(Colors.Purple);
//            controller.addHint(new Circle(b3.left(base.n3).location, 12)).fromModel().color(Color.BLUE);
//            controller.addHint(new Circle(b3.right(base.n3).location, 12)).fromModel().color(Colors.Purple);
//        }

//        // position the single-stranded nucs on the 5' side of the helix
//        // i.e. between the 3'nuc of the 5'-side helix and the 5' nuc of THIS helix
        if (l5 != null)
            positionLoopNucs(new Motif.Segment(b5 == null ? l5.start : b5.right(base.n5), base.n5), circ, flipped5to3);
//        // position the single-stranded nucs on the 3' side of the helix
//        // i.e. between the 3' nuc of THIS helix and the 5'nuc of the 3'-side helix
        if (l3 != null && l3 != l5)
            positionLoopNucs(new Motif.Segment(base.n3, b3==null? l3.end : b3.left(base.n3)), circ, flipped5to3);

        controller.controlsUpdated();
    }

    private void positionLoopNucs(Motif.Segment segment, Circle circ, boolean flipped) {
        Point2D center = circ.center();
        // segment includes fixed start- and end-point nucs.
        Vec2D a5 = new Vec2D(center, segment.start.location);
        Vec2D a3 = new Vec2D(center, segment.end.location);
        double angExtent = a5.arcAngleTo(a3, !flipped);
        double angStart = a5.getAngle();
        //Vec2D mid = a5.getRotated(angExtent/2);
        int count = segment.size(); // count is the number of total nucs. Only (count-2) are moved (2 ar fixed) and there are (count-1) gaps in-between them.
        //double dist = 0;
        //DrawSettings settings = controller.getSettings();
        for (int i = 1; i < count - 1; i++) {
            double theta = angStart + angExtent * i / (count-1);
            segment.getNucAt(i).location.setLocation(circ.getCircleX(theta), circ.getCircleY(theta));
            // dist += settings.nucleotideRadius(segment.getNucAt(i)) * 2;
        }

//        // Determine if the nucleotides all fit. If not, push them successively out further away from the center of the circle so they form more of a parabola-like or elipse-like shape.
//        double arclen = circ.radius() * Math.abs(angExtent);
//        dist += settings.nucleotideRadius(segment.start) + settings.nucleotideRadius(segment.end);
//        if (dist > arclen) {
//            Vec2D ray = new Vec2D();
//            for (int i = 1; i < (count+1) / 2 ; i++) {
//                Nuc np = segment.getNucAt(i-1);
//                Nuc n = segment.getNucAt(i);
//                double d = Vec2D.distance(np.location, n.location);
//                double r = settings.nucleotideRadius(np)+settings.nucleotideRadius(n);
//
//            }
//        }
    }

    private void positionNucsRotate(final Point2D.Float newTargetLoc, final SceneController.DragOpts options) {
        controller.rotateNucs(new HashSet<Nuc>(targetNuc.getBranch().getBases()), null, startLocation, location, options);
    }

    @Override
    protected void prepareDraw(final boolean isDragging) {
        if (!isValid()) return;
        if (!isDragging) {
            double scale = controller.getSettings().nucleotideRadius(null) * view.trToScreen.getScaleX() + getBoxSize() * Math.sqrt(2) / 2;
            scale(handleOffset, scale / magnitude(handleOffset));
            view.toScreen(targetNuc.location, location);
            location.translate(Math.round(handleOffset.x), Math.round(handleOffset.y));
        }
    }

    @Override
    protected void selectionChanged() {
        base = b5 = b3 = null;
        l5 = l3 = null;
        helix = targetNuc == null ? null : targetNuc.getHelix();
        branchNucs.clear();
        if (helix != null) {
            branchNucs.addAll(targetNuc.getBranch().getBases());
            base = helix.getPair(0);
            Nuc n = base.n5.getPrev(); // get the base on the 5' side of the helix
            if (n != null) {
                l5 = n.getLoop(); // the base can be part of a helix or part of a loop.
                b5 = n.getPairBond();
                if (l5 != null && b5 == null) {  // the base was part of a loop, so try to get the next pair at the 5' end of the loop
                    n = l5.getNuc5().getPrev();
                    if (n != null)
                        b5 = n.getPairBond();
                }
            }
            n = base.n3.getNext(); // get the base on the 3' side of the helix
            if (n != null) {
                l3 = n.getLoop(); // the base can be part of a helix or part of a loop.
                b3 = n.getPairBond();
                if (l3 != null && b3 == null) {  // the base was part of a loop, so try to get the next pair at the 3' end of the loop
                    n = l3.getNuc3().getNext();
                    if (n != null)
                        b3 = n.getPairBond();
                }
            }
        }
        calcHandleOffset();
        view.rayToScreen(handleOffset);
    }

    private void calcHandleOffset() {
        if (targetNuc == null || targetNuc.getPaired() == null) return;
        Vec2D dir = new Vec2D(targetNuc.getPaired().location, targetNuc.location);
        handleOffset.setLocation(dir);
        //handleOffset.setLocation(dir.getPointAlongScaled(targetNuc.location, 1));
    }

    private static final SceneUpdateInfo SlideBranch = SceneUpdateInfo.Layout.subType("Rotated/Adjusted Branch");
    @Override
    public SceneUpdateInfo getCompletionEvent() { return SlideBranch; }
}
