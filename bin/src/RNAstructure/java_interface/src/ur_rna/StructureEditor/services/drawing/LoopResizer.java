package ur_rna.StructureEditor.services.drawing;

import ur_rna.StructureEditor.Program;
import ur_rna.StructureEditor.models.*;
import ur_rna.StructureEditor.services.RnaDrawController.Colors;
import ur_rna.StructureEditor.services.SceneController;
import ur_rna.Utilities.annotation.Nullable;
import ur_rna.Utilities.geom.Ellipses;
import ur_rna.Utilities.geom.Ellipses.Circle;
import ur_rna.Utilities.geom.PointMath;
import ur_rna.Utilities.geom.Vec2D;
import ur_rna.Utilities.swing.ImageUtil;

import javax.imageio.ImageIO;
import java.util.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Arc2D;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.List;

import static ur_rna.Utilities.geom.PointMath.dist;
import static ur_rna.Utilities.geom.PointMath.midpoint;

/**
 * A DrawHandle that allows a user to resize an RNA Loop (a single stranded region) in a uniform way.
 */
public class LoopResizer extends DrawHandle {
    private Motif.Loop loop;
    private Motif.Loop.LoopType type;
    private Vec2D handleOffsetDir = new Vec2D(); // direction of the normal vector for the target nucleotide
    private static double radToDeg = 180/Math.PI;
    private BufferedImage[] directionalIcons = new BufferedImage[8]; // N,NE,E,SE,S,SW,W,NW  (clockwise because graphics coordinates are left-handed due to reversed Y-axis)
    private Motif.MultiLoop multiLoop;
    private BufferedImage resizeSegmentIcon, resizeMultiLoopIcon;

    private BufferedImage getDirectionalIcon(double screenAngle) {
        // graphics coordinates are left-handed due to reversed Y-axis, so NORTH is -90 degrees instead of 90 and rotation goes clockwise.
        final double arc = 2*Math.PI/directionalIcons.length; // the arc-length (in radians) of each direction segment.
        // -90 degrees ==> NORTH ==> [0]  (actually goes from (-90 - arc) to (-90 + arc)
        // -45 degrees ==> N/E   ==> [1]
        //   0 degrees ==> EAST  ==> [2]
        // ...
        // Adding PI/2 means position 0 will correspond to NORTH (-90 degrees) instead of EAST (0 degrees)
        // Adding arc/2 means each section goes from -arc/2 to arc/2 instead of from 0 to arc  (e.g. EAST goes from  -22.5 degrees to 22.5 degrees instead of from 0 to 45 degrees).
        int pos = (int)Math.floor((screenAngle+Math.PI/2+arc/2)/arc);
        pos %= 8; if (pos < 0) pos += 8;
        // String[] dirNames = "N,NE,E,SE,S,SW,W,NW".split(",");
        // System.out.println("Direction: " + ((screenAngle*radToDeg)%360) + "==>" + pos + "==>" + dirNames[pos]);
        return directionalIcons[pos];
    }
    private void createDirectionalIcons(BufferedImage northIcon) {
        directionalIcons[0] = northIcon;
        double angle = 2 * Math.PI / directionalIcons.length;
        int w = northIcon.getWidth(), h = northIcon.getHeight();

        AffineTransform tr = new AffineTransform();
        for(int i=1; i < directionalIcons.length; i++) {
            tr.setToRotation(i*angle, w/2.0, h/2.0);
            BufferedImage ri = ImageUtil.createCompatibleImage(w, h, true);
            Graphics2D g = ri.createGraphics();
            g.drawImage(northIcon, tr, null);
            g.dispose();
            directionalIcons[i] = ri;
        }
//        for(int i = 0; i < directionalIcons.length; i++)
//            try {
//                ImageIO.write(directionalIcons[i], "PNG", new File("LoopResizer_" + i + ".png"));
//            } catch (IOException ex) {
//                ex.printStackTrace();
//            }
    }
    public LoopResizer(final SceneController controller) {
        super(controller);
        resizeSegmentIcon = Program.getImage("loop-resize-handle");
        resizeMultiLoopIcon = Program.getImage("multiloop-resize-handle");
        createDirectionalIcons(resizeSegmentIcon);
        super.setIcon(directionalIcons[0]);
    }

    @Override
    public boolean isValid() {
        return loop != null;
    }

    @Override
    protected void nucPositionsChanged() {
       calculateOffset();
    }
    private Point2D.Float _ptModel = new Point2D.Float();
    /** Called when this handle has been dragged by the user. */
    @Override
    protected void performDrag(Point start, Point prev, Point current, final SceneController.DragOpts options) {
        moveLocation(start, current);
        getModelSpaceDeltaTransform(start, current).transform(targetNuc.location, _ptModel);
        positionNucs(_ptModel, options);
    }
    private void positionNucs(final Point2D.Float newPos, final SceneController.DragOpts options) {
        Nuc n1 = loop.getNuc5();
        Nuc n2 = loop.getNuc3();
        Nuc p1 = n1.getPrev();
        Nuc p2 = n2.getNext();
        Nuc nStart = p1 == null ? n1 : p1;
        Nuc nEnd = p2 == null ? n2 : p2;

        boolean enableSnap = controller.programSettings().SnapGuides ^ options.special;

        final Nuc n = targetNuc;
        switch (type) {
            case Terminal:
                if (n == nStart)
                    // terminal loop with user moving terminal nucleotide.
                    positionTerminalEndPoint(newPos, nStart, nEnd, 1, enableSnap);
                else if (n == nEnd)
                    positionTerminalEndPoint(newPos, nEnd, nStart, -1, enableSnap);
                else
                    positionTerminalLoop(newPos, nStart, nEnd, enableSnap);
                break;
//            case Hairpin:
//                positionLoopBetweenHelices(newPos, nStart, nStart.getPaired(), enableSnap);
//                break;
            default:
                if ((options.expandMotif ^ controller.programSettings().ExpandMotif) && multiLoop != null)
                    // Hairpins are specifically excluded by selectionChanged()
                    expandMultiLoop(controller, multiLoop, targetNuc.location, newPos, null, options);
                else
                    positionLoopBetweenHelices(newPos, nStart, nEnd, enableSnap);
                break;
        }
        controller.layoutUpdated(ResizeLoop, false);
    }
    public static void expandMultiLoop(SceneController ctrl, Motif.MultiLoop mloop, Point2D.Float ptFrom, Point2D.Float ptDest, @Nullable Circle loopCircle, SceneController.DragOpts options) {
        final float screenScale = (float)ctrl.getCanvas().getView().trToScreen.getScaleX(); // used for DrawHints and "sticky" determination.
        final float textOffset = 20 / screenScale;

        Bond trunk = mloop.getTrunk();
        Point2D mid = trunk.midpoint(null);
        Vec2D norm = trunk.normal(null);
        if (trunk.getScene().drawFlipped) norm.negate();

        //  ri^2 = (xi-x0)^2 + (yi-y0)^2
        //  ri^2 = (xi-x0)^2 + (yi-y0)^2

        if (loopCircle==null) {
            // Find the best circle that describes the EXISTING multi-loop
            // Use the existing targetNuc along with each helix and calculate a best-fit.
            List<Bond> helices = mloop.getBonds();
            Point2D[] pts = new Point2D[helices.size()+1];
            for (int i = 0; i < pts.length-1; i++) {
                pts[i] = helices.get(i).midpoint(null);
            }
            pts[pts.length-1] = ptFrom;
            loopCircle = Circle.calcLeastSquaresFit(pts);
            if (loopCircle == null) {
                ctrl.addHint("Circle calculation failed.", ptDest.x, ptDest.y + textOffset).color(Color.RED).fromModel();
                return;
            }
        }
        // define the new destination circle based on the new point
//        Vec2D cnorm = new Vec2D(mid, loopCircle.center()); // define a normal vector pointing along
//        if (cnorm.dot(norm)<0) {
//            if (!options.disableHints)
//                ctrl.addHint(loopCircle).fromModel().color(Color.RED); // error--bad circle.
//        }

        Vec2D ptReflect = Vec2D.reflect(ptDest, norm, mid);
        Circle circDest;
        if (ptReflect.distanceSq(ptDest)<1||mid.distanceSq(ptDest)<1)
            circDest = new Circle(Vec2D.getMidpoint(mid, ptFrom), mid.distance(ptDest)/2);
        else
            circDest = Circle.calcFromPoints(mid, ptDest, ptReflect);

        if (circDest==null) {
            ctrl.addHint(loopCircle).fromModel().color(Color.RED); // error--bad circle.
            ctrl.addHint("Circle calculation failed.", ptDest.x, ptDest.y + textOffset).color(Color.RED).fromModel();
            return;
        }

        boolean sticky = false;
        final float MIN_RADIUS = ctrl.getSettings().nucRadius, MAX_RADIUS = ctrl.getSettings().nucRadius * 1000,  MAX_ANGLE = (float)(88*Math.PI/180);
        Vec2D vdest = new Vec2D(mid, ptDest);
        if (circDest.radius() < MIN_RADIUS)
            circDest.setCircle(norm.getPointAlong(mid, MIN_RADIUS), MIN_RADIUS);
        else if (norm.angleBetween(vdest) > MAX_ANGLE) {
            circDest.setCircle(norm.getPointAlong(mid, MAX_RADIUS), MAX_RADIUS);
            ctrl.addHint("Loop size out of bounds.", ptDest.x, ptDest.y + textOffset).color(Color.RED).fromModel();
            ctrl.addHint(circDest).fromModel().color(Color.RED).bold(); // error--bad circle.
            return;
        } else if (options.special ^ ctrl.programSettings().SnapGuides) {
            float optimal = new NucLayout(ctrl.getSettings()).calcOptimalMultiLoopRadius(mloop);
            if (Math.abs(optimal - circDest.radius())<5/screenScale) {
                circDest.setCircle(norm.getPointAlong(mid, optimal), optimal);
                sticky = true;
            }
        }

        Point2D.Float newMid = new Point2D.Float();
        AffineTransform tr = new AffineTransform();
        if (!options.disableHints) {
            ctrl.addHint(loopCircle).fromModel().color(Colors.warm).dashed();
            if (sticky)
                ctrl.addHint(circDest).fromModel().bold();
            else
                ctrl.addHint(circDest).fromModel();
            ctrl.addHint(new Circle(circDest.center(), 2)).fill(Colors.cold).fromModel();
        }

        double dR = circDest.radius() - loopCircle.radius();
        Vec2D tmp = new Vec2D(); Point2D.Float c = loopCircle.center(), cDest = circDest.center();
        for(Motif.Loop loop : mloop.getLoops())
            for(Nuc n : loop.getBases())
                expandPointAboutCenter(tmp, n.location, c, cDest, dR);
        for(Bond b : mloop.getBonds()) {
            if (b == trunk) continue;
            b.midpoint(mid);
            newMid.setLocation(mid);
            expandPointAboutCenter(tmp, newMid, c, cDest, dR);
            newMid.x -= mid.getX();
            newMid.y -= mid.getY();
            tr.setToTranslation(newMid.x, newMid.y);
            b.n5.getBranch().transformAll(tr);
        }
    }
    private static void expandPointAboutCenter(Vec2D ray /* reusable dummy vector */, Point2D pt, Point2D oldCenter, Point2D newCenter, double dR) {
        //controller.addHint(new Line2D.Float(pt, oldCenter)).fromModel().color(Colors.hot    );
        ray.setTo(pt).subtr(oldCenter);
        double dist = ray.length() + dR;
        //controller.addHint(ray.toLine(newCenter)).fromModel().color(Colors.grass);
        pt.setLocation(ray.getPointAlong(newCenter, dist < 0 ? 0 : dist));
    }

    private void positionLoopBetweenHelices(final Point2D.Float ptNew, final Nuc nBase1, final Nuc nBase2, final boolean enableSnap) {
        Point2D p1, p2, m1, m2;
        p1 = nBase1.location;
        p2 = nBase2.location;

        boolean isHairpin = nBase1.getPaired() == nBase2;
        // If this is a hairpin then require that the layout arc includes the positions of both nucleotides in the closing base pair:
//        if (isHairpin) {
//
//        } else {
            // The loop is terminated on the right and left by paired nucleotides.
            //     nBase1 is on the "left" (i.e. its index is lower than any in the loop)
            //     and nBase2 is on the "right" (i.e. its index is higher than any in the loop)
            // The locations of these two nucleotides ARE NOT directly used to calculate the layout arc. Instead the mid-points between them and their paired nucleotides are used for this --i.e. the midpoints of the bonds.
            m1 = midpoint(p1, nBase1.getPaired().location);
            m2 = midpoint(p2, nBase2.getPaired().location);

//        }

        Point2D mid = midpoint(p1, p2); // midpoint between p1 and p2
        Vec2D vNorm = new Vec2D(p1, p2).rotate90().normalize(); //points along axis of helix
        Vec2D vUser = new Vec2D(mid, ptNew); // points either parallel or anti-parallel to vNorm, depending on whether user's point was above or below the base-pair.
        Vec2D vUserNorm = vUser.projectedOn(vNorm).normalize(); // points either parallel or anti-parallel to vNorm, depending on whether user's point was above or below the base-pair.

        if (isHairpin) {
            // If there is only one base-pair, the loop could be "above" or "below" it (i.e. in the direction of vNorm or opposite it).
            //     N  N       -----B--P-----
            //   N      N        N      N
            //  N        N      N        N            ^  vNorm
            //   N      N        N      N             |
            //-----B--P-----       N  N         ------|-----
            // However, if an adjacent base-pair exists, the loop must be on the opposite side of it (which could still be above or below the bond)
            Nuc nAdj = nBase1.getPrev();
            if (nAdj != null && nAdj.isPaired() && nAdj.getPaired() == nBase2.getNext()) {
                Vec2D vAdj = new Vec2D(mid, midpoint(nAdj.location, nAdj.getPaired().location)); // (could have been done with the vector normal to the Adj bond, but this would add ambiguity in the case the the adj bond is "twisted" with respect to this one. So using the midpoint is better.
                if (vAdj.dot(vUserNorm) > 0)
                    vUserNorm.negate(); // vNorm should point AWAY from the adjacent base-pair.
                // TODO: possibly still allow flipping the loop if the user is in free-draw mode (ie. when snap is disabled)
            }
        }

        int count = loop.size();
        boolean snap = false;
        float zoom = (float) view.trToScreen.getScaleX();
        DrawSettings settings = controller.getSettings();
        float nucRad = settings.nucRadius;

        // Calculate a circle that contains both base nucs as well as the new location (defined by the user's mouse drag)
        Circle circUser = Circle.calcFromPoints(p1, p2, ptNew);
        if (circUser == null || (!isHairpin && enableSnap && Math.abs(vUser.dot(vNorm)) < SNAP_DISTANCE / zoom)) {
            // points are co-linear, so arrange on a line
            positionLinear(nBase1, nBase2, true);
            return;
        }
        //System.out.println(String.format("circUser: %s %s [%s]" , circUser.center(), circUser.radius(), circUser.getFrame()));

        // controller.addHint(circUser).fromModel().color(Color.RED);
        // testGetRadius();

        double chordLength = PointMath.dist(p1, p2);

        if (enableSnap) {
            Circle circIdeal = null;
            // Determine the ideal circle.
            switch (type) {
                case Hairpin:
                    // Calculate the partial circumference, which includes the loop nucleotides, but not the helix base-pair
                    double idealArcLength =  nucRad * 2 * (count + 1)  // add one to count because we have a half-nuc from each of the base-pair nucleotides.
                            + settings.loopSpacing * (count + 1); // there is a gap BEFORE the first loop nuc, BETWEEN each loop nuc (i.e. count-1) and after the last NUC = 1 + (count -1) + 1
                    circIdeal = new Circle(Circle.getRadiusFromArcAndChord(idealArcLength, chordLength));
                    circIdeal.center(vUserNorm.getPointAlong(mid, getChordHeight(idealArcLength, chordLength, circIdeal.radius())));
                    break;
                case Multibranch: {
                    // Get the best circle that could pass through p1 and p2, whilst having a center close to the
                    // average center from each other base-pair in the multi-branch loop
                    Point2D.Float avgCenter = new Point2D.Float();
                    int branchCount = 0;
                    Collection<Bond> branches = Motif.Loop.getMultiBranchBasePairs(loop.getNuc3(), false);
                    if (branches != null)
                    for (Bond branch : branches) {
                        Point2D m3 = midpoint(branch.getNuc5().location, branch.getNuc3().location);
                        if (branch == nBase1.getPairBond() || branch == nBase2.getPairBond()) {
                            controller.addHint(Ellipses.fromCenter(m3, nucRad / 2)).color(null, Colors.cold).fromModel();
                            continue;
                        } else
                            controller.addHint(Ellipses.fromCenter(m3, nucRad / 2)).color(null, Colors.cool).fromModel();
                        Circle c = Circle.calcFromPoints(p1, p2, m3);
                        if (c == null) continue;
                        avgCenter.x += c.getCenterX();
                        avgCenter.y += c.getCenterY();
                        branchCount++;
                    }
                    if (branchCount > 0) {
                        avgCenter.x = avgCenter.x / branchCount;
                        avgCenter.y = avgCenter.y / branchCount;
                        circIdeal = new Circle(avgCenter, dist(avgCenter, p1));
                    }
                }
                    break;
            }
            if (circIdeal == null)
                circIdeal = new Circle(midpoint(m1, m2), dist(m1, m2) / 2);

             //controller.addHint(Ellipses.fromCenter(circIdeal.center(), 2)).fill(Color.PINK).noLine().fromModel();
             //controller.addHint(circIdeal).color(Color.PINK).fromModel();

            snap = PointMath.dist(circUser.center(), circIdeal.center()) <= Math.max(2, SNAP_DISTANCE / Math.sqrt(zoom));
            if (snap)
                circUser = circIdeal;
        }

        // enforce a minimum size
        if (isHairpin) {
            // Calculate the minimum arc length required to accommodate the nucleotides in the loop (used for hairpins).
            // Includes a half-nucleotide width for the base-pair and each loop nuc
            double minArcLength = nucRad * (count + 0.5f);
            Circle circMin = minArcLength > chordLength ? new Circle(Circle.getRadiusFromArcAndChord(minArcLength, chordLength)) : null;
            if (circMin != null) {
                circMin.center(vUserNorm.getPointAlong(mid, getChordHeight(minArcLength, chordLength, circMin.radius())));
                // controller.addHint(circMin).color(new Color(0x009900)).fromModel();
                if (vUserNorm.dot(circUser.center(), mid) < vUserNorm.dot(circMin.center(), mid) || vUserNorm.dot(ptNew, mid) < nucRad)
                    circUser = circMin;
            } else if (vUserNorm.dot(ptNew, mid) < nucRad && enableSnap) {
                positionLinear(nBase1, nBase2, true);
                return;
            }
        }

        Point2D center = circUser.center();
        Vec2D vStartAngle = new Vec2D(center, p1);
        Vec2D vEndAngle = new Vec2D(center, p2);

//        boolean centerAboveMid = vNorm.dot(center, midpoint(nBase1.location, nBase2.location)) > 0; // Determine whether the arc should go from v1 to v2 or reversed.
//        boolean reverseNorm = vNorm.dot(vUserNorm) < 0;
        boolean reverseArc = vNorm.dot(vUserNorm) > 0; // centerAboveMid != reverseNorm;
        double angExtent = vStartAngle.arcAngleTo(vEndAngle, reverseArc);
        double angStart = vStartAngle.getAngle();
        //System.out.println(" centerAboveMid: " + centerAboveMid + " reverseNorm: " + reverseNorm + " reverseArc: " + reverseArc);
        //System.out.println(" v1-angle: " + vStartAngle.getAngle()*radToDeg + " v2=angle" + vEndAngle.getAngle()*radToDeg + " trueAngle: " + vStartAngle.angleTo(vEndAngle) * radToDeg + "  arc: " + angExtent * radToDeg);

        //controller.addHint(vStartAngle.toLine(center)).color(Color.BLUE).fromModel();
        //controller.addHint(vEndAngle.toLine(center)).color(Color.CYAN).fromModel();
        //controller.addHint(vNorm.clone().scale(20).toLine(mid)).color(Color.DARK_GRAY).fromModel();
        //controller.addHint(vUserNorm.clone().scale(20).toLine(mid)).color(Color.LIGHT_GRAY).fromModel();
        //controller.addHint(circUser).color(reverseArc ? Color.GREEN : new Color(0x008800)).fromModel();

        controller.addHint(circUser).dashed().fromModel();
        if (snap)
            controller.addHint(circUser.arc(angStart, angExtent)).fromModel().bold();


        for (int i = 1; i < count + 1; i++) {
            double theta = angStart + angExtent * i / (count+1);
            nBase1.getNext(i).location.setLocation(circUser.getCircleX(theta), circUser.getCircleY(theta));
        }
    }

//    private void positionHairpin(final Point2D.Float ptNew, final Nuc nBase, final boolean enableSnap) {
//        // nBase and nPair are the "base nucs" -- the nucleotides in the base-pair directly adjacent to the hairpin loop.
//        // The index of nBase is always lower than nPair.
//        Nuc nPair = nBase.getPaired();
//
//        Point2D.Float mid = (Point2D.Float)midpoint(nBase.location, nPair.location, new Point2D.Float()); // midpoint between helix nucs (nBase and nPair)
//        //Vec2D base = new Vec2D(nBase.location, nPair.location).normalize();
//        Vec2D vNorm = new Vec2D(nBase.location, nPair.location).rotate90().normalize(); //points along axis of helix
//        Vec2D vUserNorm = new Vec2D(mid, ptNew).projectedOn(vNorm).normalize(); // points either parallel or anti-parallel to vNorm, depending on whether user's point was above or below the base-pair.
//
//        // If there is only one base-pair, the loop could be "above" or "below" it (i.e. in the direction of vNorm or opposite it).
//        //     N  N       -----B--P-----
//        //   N      N        N      N
//        //  N        N      N        N            ^  vNorm
//        //   N      N        N      N             |
//        //-----B--P-----       N  N         ------|-----
//
//        // If an adjacent base-pair exists, the loop must be on the opposite side of it (which could still be above or below the bond)
//        Nuc nAdj = nBase.getPrev();
//        if (nAdj != null && nAdj.isPaired() && nAdj.getPaired() == nPair.getNext()) {
//            Vec2D vAdj = new Vec2D(mid, midpoint(nAdj.location, nAdj.getPaired().location)); // (could have been done with the vector normal to the Adj bond, but this would add ambiguity in the case the the adj bond is "twisted" with respect to this one. So using the midpoint is better.
//            if (vAdj.dot(vUserNorm) > 0)
//                vUserNorm.negate(); // vNorm should point AWAY from the adjacent base-pair.
//            // TODO: possibly still allow flipping the loop if the user is in free-draw mode (ie. when snap is disabled)
//        }
//
//        // Calculate a circle that contains both base nucs as well as the new location (defined by the user's mouse drag)
//        Circle circUser = calcCircle(nBase.location, nPair.location, ptNew);
//        if (circUser == null) return;  // points are co-linear
//        //System.out.println(String.format("circUser: %s %s [%s]" , circUser.center(), circUser.radius(), circUser.getFrame()));
//
//        // controller.addHint(circUser).fromModel().color(Color.RED);
//        // testGetRadius();
//
//        int count = loop.getLength();
//        boolean snap = false;
//        float zoom = (float) view.trToScreen.getScaleX();
//
//        IDrawSettings settings = controller.getSettings();
//        double minArcLength = settings.nucleotideRadius(null) * (count + 0.5f);  // includes a half-nucleotide width for the base-pair and each loop nuc
//        double chordLength = PointMath.dist(nBase.location, nPair.location);
//        Circle circMin = new Circle(getRadiusFromArcAndChord(minArcLength, chordLength));
//        circMin.center(vUserNorm.getPointAlong(mid, getChordHeight(minArcLength, chordLength, circMin.radius())));
//        // controller.addHint(circMin).color(new Color(0x009900)).fromModel();
//
//        if (enableSnap) {
//            // Determine the ideal circle.
//            // Calculate the partial circumference, which includes the loop nucleotides, but not the helix base-pair
//            double idealArcLength = settings.nucleotideRadius(null) * 2 * (count + 1)  // add one to count because we have a half-nuc from each of the base-pair nucleotides.
//                                 + settings.nucleotideSpacing() * (count + 1);        // there is a space BEFORE the first loop nuc, BETWEEN each loop nuc (i.e. count-1) and after the last NUC = 1 + (count -1) + 1
//            // Calculate the chord length, which is the distance between the paired bases.
//
//            Circle circIdeal = new Circle(getRadiusFromArcAndChord(idealArcLength, chordLength));
//            Vec2D pt = vUserNorm.getPointAlong(mid, getChordHeight(idealArcLength, chordLength, circIdeal.radius()));
//            // controller.addHint(Ellipses.fromCenter(pt, 8)).color(Color.BLACK).fromModel();
//            circIdeal.center(pt);
//
//            //  controller.addHint(circIdeal).color(new Color(0xFF66FF)).fromModel();
//            snap = PointMath.dist(circUser.center(), circIdeal.center()) <= Math.max(2, SNAP_DISTANCE / Math.sqrt(zoom));
//            if (snap)
//                circUser = circIdeal;
//        }
//
//        if (vUserNorm.dot(circUser.center(), mid) < vUserNorm.dot(circMin.center(), mid) || vUserNorm.dot(ptNew, mid) < settings.nucleotideRadius(null))
//            circUser = circMin;
//
//        //DrawHint dh = controller.addHint(Ellipses.fromCenter(center, r)).fromModel().line(2 / zoom);
////        if (snap)
////            dh.line(2 / zoom).color(Colors.cold);
//
//        Vec2D vStartAngle = new Vec2D(circUser.center(), nBase.location);
//        Vec2D vCircCenter = new Vec2D(mid, circUser.center());
//
//        boolean flip = vCircCenter.dot(vNorm) > 0;
//        //controller.addHint(vStartAngle.toLine(circUser.center())).color(Color.CYAN).fromModel();
//        controller.addHint(vNorm.clone().scale(20).toLine(ptNew)).color(Color.DARK_GRAY).fromModel();
//        controller.addHint(vUserNorm.clone().scale(20).toLine(ptNew)).color(Color.LIGHT_GRAY).fromModel();
//        //controller.addHint(vCircCenter.toLine(mid)).color(flip ? Color.GREEN : Color.RED).fromModel();
//
//        double angStart = vStartAngle.getAngle();
//        double angExtent = vStartAngle.arcAngleTo(new Vec2D(circUser.center(), nPair.location), flip);
//        System.out.println("angExtent: " + angExtent * 180 / Math.PI);
//
////        DrawHint arc = controller.addHint(circUser.arc(angStart, angExtent)).fromModel().line(2 / zoom);
////        if (snap)
////            arc.line(3 / zoom).color(Colors.cold);
//
//        if (snap)
//            controller.addHint(circUser.arc(angStart, angExtent)).fromModel().line(3 / zoom);
//
//        for (int i = 1; i < count + 1; i++) {
//            double theta = angStart + angExtent * i / (count+1);
//            nBase.getNext(i).location.setLocation(circUser.getCircleX(theta), circUser.getCircleY(theta));
//        }
////        controller.addHint(Ellipses.fromCenter(mid, 4)).fromModel();
////        controller.addHint(base.clone().scale(x).toLine(mid)).color(Color.ORANGE).fromModel();
////        controller.addHint(norm.clone().scale(y).toLine(mid)).color(Color.PINK).fromModel();
////        controller.addHint(user.toLine(mid)).color(Color.GREEN).fromModel();
//    }

    private double getChordHeight(double arcLength, double chordLength, double radius) {
        double height = Math.sqrt(radius * radius - chordLength * chordLength / 4); // formula for height of isosceles triangle with vertex at center of circle.
        return  arcLength < Math.PI * radius ? -height : height; // if arcLength < half the circumference, the triangle is "uspside down" and we need to shift the circle down instead of up to position the chord relative to the center.
    }

//    private static boolean tested = false;
//    private void testGetRadius() {
//        if (tested) return;
//        System.out.println("********* testGetRadius *********");
//        //System.out.println("angle\tradius\tarclength\tfound-r\terr");
//
//
//        double circ = 10000;
//        double r = circ / (2 * Math.PI);
//        //double chord  = 1;
//        System.out.println("angle\tr\tchord\tarclength\tcirc\tfound-r\terr");
//
//        for (int i = 1; i < 180 * 4; i++) {
//            double angle = i * Math.PI / (180 * 4);
//            double chord = 2 * r * Math.sin(angle);
////            double r = 1 / (2 * Math.sin(angle));
//            double arc = 2 * r * (Math.PI - angle);
//            double foundRadius = getRadiusFromArcAndChord(arc, chord);
//            double err = Math.abs(r - foundRadius) / r;
//            System.out.println(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s", angle, r, chord, arc, circ, foundRadius, err));
//        }
//        tested = true;
//    }

    private void positionTerminalLoop(final Point2D.Float ptNew, final Nuc nStart, final Nuc nEnd, final boolean enableSnap) {
        final double radToDeg = -180 / Math.PI;
        final float zoom = (float)view.trToScreen.getScaleX();

        Vec2D vNorm = new Vec2D(nStart.location, nEnd.location).rotate90().normalize();
        Point2D.Double mid = (Point2D.Double)midpoint(nStart.location, nEnd.location, new Point2D.Double());
        Vec2D vUser = new Vec2D(mid, ptNew);

        Circle circ = Circle.calcFromPoints(nStart.location, nEnd.location, ptNew);
        if (circ == null || (enableSnap && Math.abs(vUser.dot(vNorm)) < SNAP_DISTANCE / Math.sqrt(zoom))) {
            positionLinear(nStart, nEnd, true);
            return;  // points are co-linear
        }
        Point2D.Float center = circ.center();
        double r = circ.radius();

        boolean snap = enableSnap && Point2D.distance(center.x, center.y, mid.x, mid.y) <= Math.max(2, 2 * SNAP_DISTANCE / zoom);
        if (snap) {
            center.setLocation(mid);
            r = PointMath.dist(center, nStart.location);
            // controller.addHint(Ellipses.fromCenter(center, r)).fromModel();}
        }

        Vec2D v1 = new Vec2D(center, nStart.location), v2 = new Vec2D(center, nEnd.location);
        boolean reverseArc = vUser.dot(vNorm) > 0;

        // controller.addHint(norm.toLine(midpoint(nStart.location, nEnd.location))).color(Color.CYAN).fromModel();
        double angStart = v1.getAngle(), angExtent = v1.arcAngleTo(v2, reverseArc);

        controller.addHint(Ellipses.fromCenter(center, 2)).fill(DrawHint.HintColorReference).noLine().fromModel();
        DrawHint arg = controller.addHint(new Arc2D.Double(center.x-r, center.y-r, r*2, r*2, angStart * radToDeg, angExtent * radToDeg, Arc2D.OPEN)).fromModel();
        if (snap)
            arg.color(DrawHint.HintColorReference).bold();

//            controller.addHint(v1.toLine(center)).color(Color.PINK).fromModel();
//            controller.addHint(v2.toLine(center)).color(Color.RED).fromModel();

//            Point2D p = v1.clone().scale(0.5).add(center).toPoint();
//            controller.addHint(""+Math.round(angStart* radToDeg), p.getX(), p.getY()).color(Color.RED).fromModel();
//            p = v2.clone().scale(0.5).add(center).toPoint();
//            controller.addHint(""+Math.round(v2.getAngle() * radToDeg  ), p.getX(), p.getY()).color(Color.RED).fromModel();
//            p = v1.clone().add(v2).scale(0.25).add(center).toPoint();
//            controller.addHint(""+Math.round(angExtent * radToDeg) + ", " + (flip ? "FLIP" : "--"), p.getX(), p.getY()).color(Color.ORANGE).fromModel();
//            controller.addHint(""+Math.round(angExtent * radToDeg) + ", " + (flip ? "FLIP" : "--") + ", v1->v2 " + (v1.angleTo(v2)*radToDeg) + ", v2->v1 " + (v2.angleTo(v1)*radToDeg), ptNew.getX()+ 10, ptNew.getY()+ 15).color(Color.BLACK).fromModel();

        // controller.addHint(new Arc2D.Double(xc-r, yc-r, r*2, r*2, a1, adiff, Arc2D.OPEN)).color(snap ? Colors.cold : Colors.cool).fromModel();

        int gaps = nEnd.indexInStrand() - nStart.indexInStrand();
        double theta;

//        boolean reversed = nStart.getPrev() == null;
//        if (!enableSnap) {
//            Vec2D vNew = new Vec2D(center, ptNew);
//            angExtent = v1.angleTo(vNew);
//            //controller.addHint(vNew.toLine(center)).fromModel();
//
//            if (flip && angExtent > 0)
//                angExtent -= 2 * Math.PI;
//            else if (!flip && angExtent <= 0)
//                angExtent += 2 * Math.PI;
//            if (reversed) {
//                angStart = v2.getAngle();
//                angExtent = -angExtent;
//            }
//        }
        for (int i = 1; i < gaps; i++) {
            //if (enableSnap)
                theta = angStart + angExtent * i / gaps;
//            else
//                theta = angStart + angExtent * i / (targetNuc.index - nStart.index);
                //controller.addHint(Ellipses.fromCenter(xc + r * Math.cos(thetaNuc), yc + r * Math.sin(thetaNuc), 2)).color(Colors.hot, 0x33).fromModel();
            nStart.getNext(i).location.setLocation(center.x + r * Math.cos(theta), center.y + r * Math.sin(theta));
        }
//        if (!enableSnap) {
//            (reversed ? nStart : nEnd).location.setLocation(center.x + r * Math.cos(angStart), center.y + r * Math.sin(angStart));
//        }


        //controller.addHint(Ellipses.fromCenter(new Point.Double(xc, yc), r)).fromModel();

//        Point2D.Float mid = new Point2D.Float();
//        midpoint(nStart.location, nEnd.location, mid);
//        Vec2D base = new Vec2D(nStart.location, nEnd.location);
//        float baseLength = (float)base.length(); // distance between the fixed terminal nucleotides.
//        int gaps = nEnd.index - nStart.index;  // number of gaps between nucleotides
//        float minSpacing = baseLength / gaps;  // minimum possible distance between nucleotides (e.g. if arranged on a line)
//        // Calculate the position where the targetNuc would be if arranged on the base line.
//        Point2D.Float ptRef = base.clone().setLength(minSpacing * (targetNuc.index - nStart.index)).add(nStart.location).toPointF();
//        controller.addHint(Ellipses.fromCenter(ptRef, 8)).color(Color.BLUE, 0x99).fromModel();
//
//        Vec2D normal = base.getRotated90();
//        Vec2D vecUser = new Vec2D(ptRef, ptNew);
//        controller.addHint(base.toLine(nStart.location)).color(Color.GREEN).fromModel();
//        controller.addHint(normal.toLine(mid)).color(Color.ORANGE).fromModel();
//        controller.addHint(vecUser.toLine(ptRef)).color(Color.PINK).fromModel();
//
//        Vec2D vVert = vecUser.projectedOn(normal);
//        Vec2D vHorz = vecUser.projectedOn(base);
//        controller.addHint(vVert.toLine(ptRef)).color(Color.YELLOW).fromModel();
//        controller.addHint(vHorz.toLine(vVert.clone().add(ptRef).toPoint())).color(Color.CYAN).fromModel();
//
//        double v = vVert.length() * 1.5;
////        double xp = vHorz.length();
////        double yp = vecUser.dot(normal) / normal.length();
////        double offset = yp - Math.sqrt((1 - xp * xp / baseLength * 2) * v);
//
//        double b = 1/(v*v);
//        double  x1 = baseLength / 2, y1 =0, r1 = 1 - b * y1 * y1,
//                x2 = vHorz.length(), y2 = vecUser.dot(normal) / normal.length(), r2 = 1 - b * y2 * y2,
//                x3 = 0, y3 = v, r3 = 1 - b * y3 * y3;
//        Matrix3D m = new Matrix3D(
//                x1 * x1, y1, 1,
//                x2 * x2, y2, 1,
//
//                );
//        m.invert();
//        double[] r = new double[ ]{ r1, r2, r3 };
//        m.transform(r, 0, 0);
//
//        double h = Math.sqrt(1/r[0]);
//        double j = r[1] / (-2 * b);
//        double jj = Math.sqrt(r[2] / b);
//
//        controller.addHint(String.format("h %s  j %s  jj %s", h, j, jj), 400, 20).color(Color.BLACK);
//
//        double a = r[0], c = r[1], d = r[2];
//        //b * j * j
//
//       // assert d == b * j * j : " bj2 = " + (b * j * j) + " d=" + d;
//        assert roundDbl(a * x1 * x1 +  b * y1 * y1 + c * y1 + d) == 1 : "result = " + (a * x1 * x1 +  b * y1 * y1 + c * y1 + d);
//        assert roundDbl(a * x2 * x2 +  b * y2 * y2 + c * y2 + d) == 1 : "result = " + (a * x2 * x2 +  b * y2 * y2 + c * y2 + d);
//        assert roundDbl(a * x3 * x3 +  b * y3 * y3 + c * y3 + d) == 1 : "result = " + (a * x3 * x3 +  b * y3 * y3 + c * y3 + d);
//
//        Point2D center = vVert.getNormalized().scale(j).add(mid).toPoint();
//        Ellipse2D e = new Ellipse2D.Double(center.getX() - h / 2, center.getY() - v, h, 2 * v);
//        Shape s = AffineTransform.getRotateInstance(base.getAngle(), center.getX(), center.getY()).createTransformedShape(e);
//
//        controller.addHint(s).fromModel();
    }

    private void positionTerminalEndPoint(final Point2D.Float ptNew, final Nuc nMoved, final Nuc nFixed, final int dir, final boolean enableSnap) {
        Point2D.Float ptFixed = nFixed.location; //new Point2D.Float();
//        if (end.isPaired()) {
//            midpoint(end.location, end.getPaired().location, pEnd); // set pEnd to be the midpoint between the two paired bases.
//        } else
        // ptEnd.setLocation(nEnd.location);

        int spaces = Math.abs(nMoved.indexInStrand() - nFixed.indexInStrand()); // start and end remain fixed. this is the number of spaces between
        if (enableSnap) {
            // if this is NOT a modified drag, then use "snap" positions:
            //  1. Snap to same X coordinate as start point
            //  2. Snap to same Y coordinate as start point
            //  3. Snap to a distance where nucleotides have "optimal" spacing
            DrawSettings settings = controller.getSettings();
            float addHintLength = settings.nucleotideRadius(nMoved) * 2;
            float snapDist = (float)(SNAP_DISTANCE / view.trToScreen.getScaleX());

            // Snap to perfectly horizontal or vertical with the fixed nuc. Also snap to the lines parallel and perpendicular to the helix
            Vec2D[] snapLines = new Vec2D[]{new Vec2D(0, 1), new Vec2D(1, 0), null, null };
            if (nFixed.isPaired()) {
                Vec2D vecHelix = new Vec2D(nFixed.location, nFixed.getPaired().location).normalize();
                if (Math.abs(vecHelix.x) > 0.1 && Math.abs(vecHelix.y) > 0.1) {
                    snapLines[2] = vecHelix;
                    snapLines[3] = vecHelix.getRotated90();
                }
            }
            Vec2D vecUser = new Vec2D(ptFixed, ptNew);
            double distUser = vecUser.length();
            for (Vec2D vecSnap : snapLines) {
                if (vecSnap == null) continue;
                double projection = vecSnap.dot(vecUser);
                if (Math.sqrt(distUser * distUser - projection * projection) <= snapDist) {
                    // The projection gives the distance along the snap line that one must travel
                    // to get to the point closest to the user's point.
                    ptNew.setLocation(vecSnap.x * projection + ptFixed.x, vecSnap.y * projection + ptFixed.y);
                    controller.addHint(vecSnap.toLine(ptFixed,  (projection < 0 ? -1 : 1) * (distUser + addHintLength))).fromModel().bold();
                }
            }

            // Snap to the optimal distance
            float optimalDistance = settings.calcOptimalNuc2NucDistance(nFixed, nMoved);
            vecUser = new Vec2D(ptFixed, ptNew);
            if (Math.abs(vecUser.length() - optimalDistance) <= snapDist) {
                ptNew.setLocation(vecUser.setLength(optimalDistance).add(ptFixed));
                controller.addHint(Ellipses.fromCenter(ptFixed, optimalDistance)).fromModel().bold();
            }
        }

        float x = (ptFixed.x - ptNew.x) / spaces, y = (ptFixed.y - ptNew.y) / spaces;
        for (int i = 0; i < spaces; i++)
            nMoved.getNext(dir * i).location.setLocation(ptNew.x + x * i, ptNew.y + y * i);
    }

    void positionLinear(Nuc nStart, Nuc nEnd, boolean boldHint) {
        positionLinear(nStart, nEnd, nStart.location, nEnd.location, false, false, boldHint);
    }
    void positionLinear(Nuc nStart, Nuc nEnd, Point2D p1, Point2D p2, boolean gapAtStart, boolean gapAtEnd, boolean boldHint) {
        int gaps = nEnd.indexInStrand() - nStart.indexInStrand();
        int dir = gaps < 0 ? -1 : 1;
        if (dir == -1) gaps = -gaps;
        if (gapAtStart) gaps++;
        if (gapAtEnd) gaps++;
        double dx = (p2.getX() - p1.getX()) / gaps, dy = (p2.getY() - p1.getY()) / gaps;
        double startX = p1.getX(), startY = p1.getY();
        int gapIndex = gapAtStart ? 1 : 0;
        Nuc n = nStart;
        while(true) {
            n.location.setLocation(startX + gapIndex * dx, startY + gapIndex * dy);
            gapIndex++;
            if (n == nEnd) break;
            n = n.getNext(dir);
        }
        DrawHint h = controller.addHint(new Vec2D(p1, p2).toLine(p1)).fromModel();
        if (boldHint) h.bold();
    }

    @Override
    protected void prepareDraw(final boolean isDragging) {
        if (!isValid()) return;
        if (!isDragging) {
            handleOffsetDir.setLength(controller.getSettings().nucleotideRadius(null) * view.trToScreen.getScaleX() + getBoxSize() * Math.sqrt(2) / 2);
            view.toScreen(targetNuc.location, location);
            location.translate((int)handleOffsetDir.x, (int)handleOffsetDir.y);
            if (multiLoop != null && controller.programSettings().ExpandMotif)
                super.setIcon(resizeMultiLoopIcon);
            else
                super.setIcon(getDirectionalIcon(handleOffsetDir.getAngle()));
        }
    }
    @Override
    protected void selectionChanged() {
        loop = targetNuc == null ? null : targetNuc.getLoop();
        multiLoop = loop == null ? null : targetNuc.getMultiLoop(2, false);
        if (loop == null) return;
        type = loop.getType();
        calculateOffset();
    }

    private void calculateOffset() {
        if (targetNuc!=null) view.rayToScreen(controller.calcNormalToNuc(targetNuc), handleOffsetDir);
    }

    @Override
    public void dragComplete(final boolean success) {
        super.dragComplete(success);
        calculateOffset();
    }

//    public Shape draw(final Graphics2D g, final View2D view, final IDrawSettings settings) {
//        nucRad = settings.nucleotideRadius(null);
//
//        Motif.Loop l = targetNuc.getLoop();
//        Nuc s = l.getNuc5(), e = l.getNuc3(); // get start end end of loop
//        Nuc sh = s.getPrev(), eh = e.getNext(); // get paired bases at either end
//
//        java.util.List<Point2D.Float> points = new ArrayList<>(l.getLength() + 2);
//        for (Nuc n : l.getBases())
//            points.add(n.location);
//
//        if (sh != null) points.add((Point2D.Float) midpoint(sh.location, sh.getPaired().location));
//        if (eh != null) points.add((Point2D.Float) midpoint(eh.location, eh.getPaired().location));
//
//        if (points.size() < 3)
//            throw new RuntimeException("Cannot determine circle");
//
//        //http://www.had2know.com/academics/best-fit-circle-least-squares.html
//        double sx2 = 0, sy2 = 0, sxy = 0, sx = 0, sy = 0, sr1 = 0, sr2 = 0;
//        for (Point2D.Float p : points) {
//            double x2 = p.x * p.x;
//            double y2 = p.y * p.y;
//            sx2 += x2;
//            sy2 += y2;
//            sxy += p.x * p.y;
//            sx += p.x;
//            sy += p.y;
//            sr1 += p.x * (x2 + y2);
//            sr2 += p.y * (x2 + y2);
//        }
//        Matrix3D mx = new Matrix3D(
//                sx2, sxy, sx,
//                sxy, sy2, sy,
//                sx, sy, points.size());
//        mx.invert();
//        double[] vec = new double[]{sr1, sr2, sx2 + sy2};
//        mx.transform(vec);
//
//        double cx = vec[0] / 2;
//        double cy = vec[1] / 2;
//        double r = Math.sqrt(4 * vec[2] + vec[0] * vec[0] + vec[1] * vec[1]) / 2;
//
//        Shape circle = new Ellipse2D.Double(cx - r, cy - r, 2 * r, 2 * r);
//        circle = view.createTransformedShape(circle);
//        g.setColor(RnaDrawController.Colors.cool);
//        g.draw(circle);
//
//        Point2D.Double center = new Point2D.Double(cx, cy);
//        g.setColor(RnaDrawController.Colors.cold);
//        g.fill(dot(center, Math.max(nucRad * 0.6, 3 / view.getScaleX()), view));
//
//        g.setColor(RnaDrawController.Colors.setAlpha(RnaDrawController.Colors.hot, 0x99));
//        if (sh != null) g.fill(nucShape(sh, view));
//        if (eh != null) g.fill(nucShape(eh, view));
//
//        g.setStroke(DrawHint.HintLineDefault);
//        g.setColor(RnaDrawController.Colors.cool);
//        g.draw(line(center, targetNuc.location, view));
//
//        setBoxPos(calcBoxPos(targetNuc.location, diff(targetNuc.location, center), nucRad, settings, view));
//
//        return drawBox(g);
//    }

    private static final SceneUpdateInfo ResizeLoop = SceneUpdateInfo.Loop; // SceneUpdateInfo.Layout.subType("Resized Loop");
    @Override
    public SceneUpdateInfo getCompletionEvent() {
        return ResizeLoop;
    }
}
