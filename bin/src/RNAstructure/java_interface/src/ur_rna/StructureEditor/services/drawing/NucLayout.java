package ur_rna.StructureEditor.services.drawing;

import ur_rna.RNAstructure.RnaBackendException;
import ur_rna.RNAstructure.backend.RNA;
import ur_rna.StructureEditor.Program;
import ur_rna.StructureEditor.models.*;
import ur_rna.StructureEditor.services.SceneController;
import ur_rna.StructureEditor.services.SceneDrawMode;
import ur_rna.StructureEditor.services.fileIO.RnaFileIO;
import ur_rna.Utilities.geom.Ellipses;
import ur_rna.Utilities.geom.Ellipses.Circle;
import ur_rna.Utilities.geom.Vec2D;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Arc2D;
import java.awt.geom.Point2D;
import java.awt.geom.QuadCurve2D;

import static java.lang.Math.PI;
import static ur_rna.Utilities.geom.Vec2D.TwoPI;

/**
 * Performs nucleotide layouts of various types.
 */
public class NucLayout {
    private DrawSettings drawSettings;

    public NucLayout() { this(new DrawSettings());  }
    public NucLayout(final DrawSettings drawSettings) {
        this.drawSettings = drawSettings;
    }

    public void redrawRadial(RnaScene scene) throws RnaBackendException {
        RNA rna = RnaFileIO.createRNAfromScene(scene, true, true, true);
        backendRedrawRadial(scene, rna, 1);
//            double[] oldPoints , newPoints;
//            oldPoints = s.getLocations((double[])null, 0);
//            newPoints = s.getLocations((double[])null, 0);
//            double[] params = new double[3];
//            if (getBestOrientation(oldPoints, newPoints, params, nucleotideRadius)) {
//                AffineTransform tr = new AffineTransform();
//                tr.translate(params[0], params[1]);
//                tr.rotate(params[2], params[3], params[4]);
//                tr.transform(newPoints, 0, newPoints, 0, 0);
//            }
    }
    private void backendRedrawRadial(RnaScene scene, RNA rna, int structureNumber) {
        final int RNA_IMPORT_MARGIN = (int) drawSettings.nucRadius;

        boolean flip = Program.getInstance().settings().DrawClockwise;

        int numBases = rna.GetSequenceLength();
        if (numBases == 0) return;
        if (!scene.hasBonds()) {
            redrawCircular(scene);
            return;
        }
        scene.drawMode = SceneDrawMode.Standard;
        scene.drawFlipped = flip;

        //DEBUG_TIMING: StopWatch timer = new StopWatch(true);

        // always remove pseudoknots before drawing.
        // TODO: Update RNAstructure to do this automatically.
        RnaFileIO.breakBackendPseudoknots(scene, rna, structureNumber);

        //DEBUG_TIMING: timer.println("BreakPseudoknot").restart();
        //rna.DetermineDrawingCoordinates(iceil(2*drawSettings.nucRadius), iceil(drawSettings.bondLength), iceil(drawSettings.nucSpacing), iceil(drawSettings.nucSpacing), structureNumber);
        final int diameter = iceil(2*drawSettings.nucRadius);
        rna.DetermineDrawingCoordinates(diameter, diameter, structureNumber);
        //DEBUG_TIMING: timer.println("DetermineDrawingCoordinates").restart();
        Rectangle bounds = new Rectangle();

        bounds.setLocation(rna.GetNucleotideXCoordinate(1), rna.GetNucleotideYCoordinate(1));

        int nucIndex = 0, strandIndex = 0;
        Strand strand = scene.strands.get(strandIndex);
        char[] sequence = rna.GetSequence().toCharArray();
        for (int i = 1; i <= numBases; i++) {
            if (RnaFileIO.isInterStrandLinker(sequence[i - 1])) { // skip over linkers
                if (nucIndex != 0) {
                    nucIndex = 0;
                    strand = scene.strands.get(++strandIndex);
                    while(strand.isEmpty()) // under some strange circumstances, a strand could be empty.
                        strand = scene.strands.get(++strandIndex);
                }
                continue;
            }
            int x = rna.GetNucleotideXCoordinate(i);
            int y = rna.GetNucleotideYCoordinate(i);
            if (flip) { x = -x; }
            strand.get(nucIndex).location.setLocation(x, y);
            bounds.add(x, y);
            nucIndex++;
        }

        //DEBUG_TIMING: timer.println("Read Coords").restart();

        // reset coordinate origin
        float dx = RNA_IMPORT_MARGIN - bounds.x, dy = RNA_IMPORT_MARGIN - bounds.y;
        for (Strand s : scene.strands)
            for (Nuc n : s)
                n.translate(dx, dy);

        //DEBUG_TIMING: timer.println("Translate").restart();
    }
    private static int iceil(double number) { return (int)Math.ceil(number); }

    public void redrawRadial(final RnaSceneGroup group) throws RnaBackendException {
        //DEBUG_TIMING: StopWatch timer = new StopWatch(true);
        RNA rna = RnaFileIO.writeRNAstructureGroup(group, true, true);
        //DEBUG_TIMING: timer.println("writeRNAstructureGroup").restart();
        for (int i = 0; i < group.size(); i++) {
            backendRedrawRadial(group.get(i), rna, 1+i);
            //DEBUG_TIMING: timer.println("layout "+i).restart();
        }
    }

    public void redrawCircular(final RnaScene scene) {
        int spots = scene.getNucCount() + scene.strands.size(); // leave a space in between the start and end of each strand.
        float circ = drawSettings.calcOptimalLoopDistance(spots, 1);  // total circumference   //// (drawSettings.nucDiameter()+drawSettings.nucSpacing) * spots;
        double radius = circ / TwoPI;  // radius of circle required for
        int flip = Program.getInstance().settings().DrawClockwise ? -1 : 1;
        final double startAngle = drawSettings.circularStartAngle;
        int pos = 0;
        for (Strand s : scene.strands) {
            for (Nuc n : s) {
                double theta = startAngle + flip * pos * TwoPI / spots;
                n.location.setLocation(radius * Math.cos(theta), radius * -Math.sin(theta));
                pos++;
            }
            pos++; // skip an extra spot between strands
        }
        scene.drawMode = SceneDrawMode.Circular;
        scene.drawFlipped = flip == -1;
    }

    public static QuadCurve2D getCircularBond(Bond b) {
        //RnaScene scene = b.getScene();
        //int spots = scene.getNucCount() + scene.strands.size();
        //int midIndex = (b.n5.indexInScene() + b.n3.indexInScene()) / 2;

        Vec2D v1 = new Vec2D(b.n5.location);
        Vec2D v2 = new Vec2D(b.n3.location);
        Vec2D mid = Vec2D.getMidpoint(v1, v2);
        double angle = v1.angleBetween(v2);

        //Vec2D pMid = (Vec2D)b.midpoint(null);
        //float part = ((b.n3.indexInScene()+b.n3.getStrandIndex()) - (b.n5.indexInScene()+b.n5.getStrandIndex())) / (float)spots;

        // scaling to 0 moves the control point to the center (0,0). scaling to 1 leaves it as-is (the mid-point)
        // angle is from 0 to PI, so 2 * angle / TwoPI goes from 0 to 1.
        // 2 * (angle*2) / (TwoPI+angle) goes from 0 to 1.
        mid.scale(Math.min(1, 1 - 4 * angle / (3 * PI)));

//        Vec2D vNorm = b.normal(null);
//        Vec2D ray = new Vec2D(scene.getNuc(midIndex).location, pMid);
//
//        boolean flipped = vNorm.dot(ray) < 0;
//        if (flipped)
//            vNorm.negate();
//
//        float circ = (settings.nucDiameter()+settings.nucSpacing) * spots; // total circumference
//        double radius = circ / TwoPI;  // radius of circle required
//
//        float part = ((b.n5.indexInScene()+b.n5.getStrandIndex()) - (b.n3.indexInScene()+b.n3.getStrandIndex())) / (float)spots;

        return new QuadCurve2D.Float(b.n5.location.x, b.n5.location.y, (float)mid.x, (float)mid.y, b.n3.location.x, b.n3.location.y);
    }

    public void redrawLinear(final RnaScene scene) {
        float step = drawSettings.nucDiameter()+drawSettings.nucSpacing;
        int pos = 0;
//        if (Program.getInstance().settings().DrawClockwise)
//            step = -step;
        for (Strand s : scene.strands) {
            for (Nuc n : s) {
                n.location.setLocation(step*pos, 0);
                pos++;
            }
            pos++; // skip an extra spot between strands
        }
        scene.drawMode = SceneDrawMode.Linear;
        scene.drawFlipped = false;
    }

    public static Arc2D getLinearBond(Bond b) {
        Vec2D mid = Vec2D.getMidpoint(b.n5.location, b.n3.location);
        Vec2D v5 = new Vec2D(mid, b.n5.location);
        Vec2D v3 = new Vec2D(mid, b.n3.location);
        double arc = v5.arcAngleTo(v3, false);
        return Ellipses.arc(mid, v5.length(), v5.getAngle(), arc, Arc2D.OPEN); //new Arc2D.Float()
    }

    public void setDrawSettings(final DrawSettings drawSettings) {
        this.drawSettings = drawSettings;
    }

    public void redrawHelix(final Motif.Helix h) {
        // get midpoints of start and end bonds in the helix
        Bond bFirst = h.first(), bLast = h.last();

        // get midpoints of start and end bonds in the helix as well as the midpoint between them
        Point2D pFirst = bFirst.midpoint(null);
        Point2D pLast = bLast.midpoint(null);
        Vec2D pMid = Vec2D.getMidpoint(pFirst, pLast);

        Vec2D dirV = new Vec2D(pFirst, pLast); // get vector that points along helix (from the midpoint of the first bond to the midpoint of the last bond)
        Vec2D dirNormal = bFirst.normal(null); // get vector normal to the first bond (in an ideal helix, this would point parallel or anti-parallel to dirV)

        if (pFirst.distanceSq(pLast) < 1) // handle the case of pFirst and pLast being at the same location.
            dirV = dirNormal;

        boolean flipped = dirNormal.dot(dirV) < 0; // determine if the helix is drawn "clockwise"/reflected (i.e. the bond-normal is anti-parallel to the helix-normal)

        int count = h.size();
        float stepLength = drawSettings.calcOptimalNuc2NucDistance(bFirst.n5, bLast.n5) / (count-1); // get ideal nuc-to-nuc separation.
        dirV.setLength(stepLength); // adding dirV to a point will advance it along the "vertical" length of the helix.
        float hDist = drawSettings.bondLength / 2 + drawSettings.nucRadius; // get horizontal distance between nucs (i.e. between paired nucs). hDist is HALF the distance across the bond
        Vec2D dirH = dirV.getRotated90().setLength(flipped?-hDist:hDist); // direction from midpoint of a bond to the 5' nuc  (and -dirV is the direction to the 3' nuc)

        pMid.add(dirV.x * -0.5*(count-1), dirV.y * -0.5*(count-1)); // move pMid to the ideal midpoint of the first bond
        Vec2D p = new Vec2D();
        for (int i = 0; i < count; i++) {
            Bond b = h.getPair(i);
            b.n5.location.setLocation(p.setTo(pMid).add(dirH));
            b.n3.location.setLocation(p.setTo(pMid).subtr(dirH));
            pMid.add(dirV); // advance along length of helix
        }
    }

    public float calcOptimalMultiLoopRadius(Motif.MultiLoop m) {
        int loopNucs = m.loopNucCount();
        int helixCount = m.helixCount();
        float chordLength = drawSettings.calcOptimalOuterDistance(2*helixCount, 0),
                arcLength = drawSettings.calcOptimalLoopDistance(loopNucs, 2);
        double optimalRadius = Circle.getRadiusFromArcAndChord(arcLength, chordLength);
        if (Double.isNaN(optimalRadius))
            optimalRadius = (arcLength + chordLength * 1.4) / Vec2D.TwoPI;
        return (float)optimalRadius;
    }

    public Circle calcLoopBestFitCircle(Motif.Loop m) {
        // get outer paired nucs, if present.
        Nuc b5 = m.getNuc5().getPrev(), b3 = m.getNuc3().getNext();
        Point2D[] pts = new Point2D[m.size()+(b5==null?0:1)+(b3==null?0:1)];
        int count = 0;
        if (b5!=null)
            pts[count++] = b5.getPairBond().midpoint();
        for(int i = 0; i < m.size(); i++)
            pts[count++] = m.getNucAt(i).location;
        if (b3!=null)
            pts[count] = b3.getPairBond().midpoint();
        return Circle.calcLeastSquaresFit(pts);
    }

    public Circle calcMultiLoopBestFitCircle(Motif.MultiLoop m) {
        Point2D[] pts = new Point2D[m.loopNucCount()+m.helixCount()];
        int ptCount = 0;
        for(Object o : m.elements)
            if (o instanceof Motif.Loop)
                for (Nuc n : ((Motif.Loop) o).getBases())
                    pts[ptCount++]=n.location;
            else if (o instanceof Bond)
                pts[ptCount++]=((Bond) o).midpoint(null);
        return Circle.calcLeastSquaresFit(pts);
    }

    public void redrawMultiLoop(final Motif.MultiLoop m, SceneController controller) {
        // count total nucs.
        // position nucs and branches.
        int loopNucs = m.loopNucCount();
        int helixCount = m.helixCount();
        float chordLength = drawSettings.calcOptimalOuterDistance(2*helixCount, 0),
            arcLength = drawSettings.calcOptimalLoopDistance(loopNucs, 2);
        //System.out.printf("Nucs:%s, Helix:%s\n", loopNucs, helixCount);
        double optimalRadius = Circle.getRadiusFromArcAndChord(arcLength, chordLength);
        if (Double.isNaN(optimalRadius))
            optimalRadius = (arcLength + chordLength * 1.4) / Vec2D.TwoPI;

        final float MIN_LOOP_RADIUS = Math.min(drawSettings.nucRadius,  drawSettings.calcOptimalOuterDistance(loopNucs/2+helixCount, 0));
        double rad;
        Circle bestFit = calcMultiLoopBestFitCircle(m);
        if (bestFit == null)
            rad = optimalRadius;
        else {
            rad = bestFit.radius();
            // If it's too far off, assume it is broken and reset to optimal.
            if (rad < MIN_LOOP_RADIUS || rad < optimalRadius / 10)
                rad = optimalRadius;
            else if (rad > optimalRadius * 100)
                rad = optimalRadius;
        }

        Bond trunk = m.getTrunk();
        Vec2D mid = trunk.midpoint();
        Vec2D normal = trunk.normal();
        final boolean flip = trunk.getScene().drawFlipped;
        if (flip) normal.negate();

        // controller.addHint(normal.toLine(mid, 100)).fromModel().color(Color.RED);

        Circle circ = new Circle(normal.getPointAlong(mid, rad), rad);

        // Calculate the angle reserved for each loop nucleotide
        double nucAngle = arcLength/optimalRadius/loopNucs;
        double helixAngle = (Vec2D.TwoPI * optimalRadius - arcLength)/optimalRadius/helixCount;
        // Calculate the angle reserved for each helix
        double angle = new Vec2D(circ.center(), mid).getAngle(); //  + helixAngle / 2
        //rdc.addHint(new Ellipses.Circle(circ.getPoint(angle) , 20)).fromModel().color(Color.BLUE);
        int dir = flip ? 1 : -1;
        //rdc.addHint(Vec2D.fromAngle(angle+dir*helixAngle/2).toLine(circ.center(), circ.radius())).fromModel().color(Color.ORANGE);
        //rdc.addHint(circ).fromModel().color(RnaDrawController.Colors.cold);

        angle+=dir*helixAngle/2; // move for Trunk branch
        Vec2D norm = new Vec2D(), translate = new Vec2D();
        AffineTransform tr = new AffineTransform();
        for(Object o : m.elements) {
            if (o instanceof Motif.Loop)
                // simply move each loop nucleotide into place
                for (Nuc n : ((Motif.Loop) o).getBases()) {
                    angle += dir*nucAngle / 2;
                    n.location.setLocation(circ.getCircleX(angle), circ.getCircleY(angle));
                    angle += dir*nucAngle / 2;
                }
            else if (o instanceof Bond) {
                // For the bonds, we need to:
                //   1) determine new x, y coordinates for the midpoint and determine translation
                //   2) determine rotation for helix
                //   3) apply rotate and translate to the entire branch
                if (o == trunk) continue;
                Bond b = (Bond) o;
                angle += dir*helixAngle / 2;
                // set translate to the new midpoint, then subtract the old midpoint to get the translate vector
                circ.getPoint(angle, translate);
                translate.subtr(b.midpoint(mid));
                Motif.Branch br = b.getScene().getBranch(b.n5);
                b.normal(norm); // these are all exit branches
                if (flip) norm.negate();
                tr.setToTranslation(translate.x, translate.y);
                tr.rotate(angle - norm.getAngle(), mid.x, mid.y);
                br.transformAll(tr);
                angle += dir*helixAngle / 2;
            }
        }
    }

    public void redrawLoop(final Motif.Loop m) {
        // get midpoints
//        Bond b = h.first();
//        Vec2D base = Vec2D.getMidpoint(b.n5.location, b.n3.location);
//        Vec2D dirH = new Vec2D(b.n5.location, b.n3.location).setLength(2*nucRadius());
//        Vec2D dirV = dirH.getRotated90().setLength(2*nucRadius()+settings.nucleotideSpacing());
//        boolean flipped = new Vec2D(base,  h.last().midpoint(null)).dot(dirV) < 0;
//        if (flipped) dirV.negate();
//
//        for (int i = 0; i < h.size(); i++) {
//            b = h.getPair(i);
//            Vec2D mid = dirV.getPointAlongScaled(base, i);
//            b.n5.location.setLocation(dirH.getPointAlongScaled(mid, -1));
//            b.n3.location.setLocation(dirH.getPointAlongScaled(mid, 1));
//        }
    }
}
