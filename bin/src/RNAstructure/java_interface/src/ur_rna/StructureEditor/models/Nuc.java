package ur_rna.StructureEditor.models;

import ur_rna.Utilities.annotation.NotNull;

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.util.Collections;
import java.util.Set;

/**
 * Represents a single nucleotide.
 */
public class Nuc implements Cloneable, INucGroup, Comparable<Nuc> {
//    public static Map<String,Object> EMPTY_MAP = Collections.emptyMap();
    public Nuc() { this(null, 0 ,0); }
    public Nuc(final String symbol) { this(symbol, 0, 0); }
    public Nuc(final String symbol, float x, final float y) {
        this.symbol = symbol;
        this.location = new Point2D.Float(x, y);
    }

    // These transient properties are set from outer scopes and are liable to change.
    // They define the "environment" of the Nuc rather than the properties of the nuc itself.
    // package-local properties -- Should be set only from StrandList or RnaScene
    transient Strand strand;
    /** Stores the index of the nucleotide within the full scene.  */
    transient int indexInScene;
    Bond pair;

    // These are the actual owned properties of the nucleotide.
    public String symbol;
    public Point2D.Float location;
    public int number = -1; // stores historic or alternate number (for annotation only. not used by this program)
    public NucStyle style;

    public RnaScene getScene() { return strand==null?null:strand.getScene(); }
    public Strand getStrand() { return strand; }
    public int getStrandIndex() { return strand==null?-1:strand.getIndex(); }
    public Bond getPairBond() { return pair; } //pairIndex == -1 ? null : scene.bonds.get(pairIndex); }
    public Nuc getPaired() { return pair == null ? null : getPairBond().getOther(this); }
    public Nuc getPaired(boolean includePseudo) { return pair != null && (includePseudo||pair.isTruePair()) ? getPairBond().getOther(this) : null; }
    public boolean isPaired() { return pair != null; }
    public boolean isPaired(boolean includePseudo) { return pair != null && (includePseudo||pair.isTruePair()); }
    /** Returns the style for this nucleotide, creating a new one if necessary. */
    public @NotNull NucStyle style() {
      if (style==null) style = new NucStyle();
        return style;
    }

    /**
     * @return The next nucleotide in this strand, or null if this is the last (3') nucleotide.
     */
    public Nuc getNext() { return strand.getNext(this); }

    /**
     * @return The nucleotide in this strand that is the specified number of nucleotides
     * (steps) away from this one.
     * A negative number of steps will return nucleotides before this one in the strand,
     * while a positive number will return nucleotides after this one.
     */
    public Nuc getNext(int steps) { return strand.getNext(this, steps); }

    /**
     * @return The previous nucleotide in this strand, or null if this is the first (5') nucleotide.
     */
    public Nuc getPrev() { return strand.getNext(this, -1); }

    /** Compares based on index in scene. Returns -1 if this Nuc comes before other, 0 if they are at the same index, or 1 if this Nuc comes after other. */
    public int compareTo(final Nuc other) {
//        if (this.strandIndex == other.strandIndex)
//            return this.index - other.index; //  e.g. 5 - 2 = 3 (positive->second comes after)
//        return this.strandIndex - other.strandIndex;
        return this.indexInScene - other.indexInScene;
    }

    /** Compares based on index in scene. Returns -1 if lhs comes before rhs, 0 if they are at the same index, or 1 if lhs comes after rhs. */
    public static int compare(final Nuc lhs, final Nuc rhs) { return lhs.indexInScene - rhs.indexInScene; } // return   lhs.compareTo(rhs); }

    public Nuc clone() {
        try {
            return (Nuc)super.clone();
        } catch (CloneNotSupportedException ex) {
            throw new InternalError(); // this shouldn't happen, since we are Cloneable
        }
    }

    //@Override
    public void translate(final float dx, final float dy) {
        location.x += dx;
        location.y += dy;
    }

    /**
     * Returns a string representation of the object. In general, the
     * {@code toString} method returns a string that
     * "textually represents" this object. The result should
     * be a concise but informative representation that is easy for a
     * person to read.
     * It is recommended that all subclasses override this method.
     * <p>
     * The {@code toString} method for class {@code Object}
     * returns a string consisting of the name of the class of which the
     * object is an instance, the at-sign character `{@code @}', and
     * the unsigned hexadecimal representation of the hash code of the
     * object. In other words, this method returns a string equal to the
     * value of:
     * <blockquote>
     * <pre>
     * getClass().getName() + '@' + Integer.toHexString(hashCode())
     * </pre></blockquote>
     *
     * @return a string representation of the object.
     */
    @Override
    public String toString() {
        return ("S"+ (strand==null?"?":""+(getStrandIndex()+1)) + ".") + (symbol == null ? "N" : symbol) + (indexInStrand()+1);
    }

    //    /**
//     * Determine if any part of the motif is under the specified point.
//     * (Note that the motif may allow a small margin around each physical element).
//     *
//     * @param point The location to test.
//     * @param settings The settings, which specify values such as nucleotide radius that
//     *                 vary depending on the graphics scale and other settings.
//     * @return True if any part of the motif is under the specified point, or false otherwise.
//     */
//    public boolean hitTest(Point2D point, final IDrawSettings settings) {
//        double r = settings.nucleotideRadius(this) + settings.getMargin(this, IDrawSettings.MARGIN_HITTEST, 0);
//        double x = point.getX()-location.getX();
//
//        if (Math.abs(x) > r) return false;
//
//        double y = point.getY()-location.getY();
//        if (Math.abs(y) > r) return false;
//
//        return x*x + y*y <= r*r;
//    }
//
//    /**
//     * Gets a {@link Rectangle2D} representing the bounding box of the motif. See {@link Shape#getBounds2D()}.
//     *
//     * @param settings The settings, which specify values such as nucleotide radius that
//     *                 vary depending on the graphics scale and other settings.
//     * @return A {@link Rectangle2D} object representing the bounding rectangle of the motif.
//     */
//    @Override
//    public Rectangle2D getBounds(final IDrawSettings settings) {
//        float margin = settings.getMargin(this, IDrawSettings.MARGIN_OUTLINE, 0);
//        float radius = settings.nucleotideRadius(this);
//        float r = radius+margin;
//        return new Rectangle2D.Float((float)location.getX()-r, (float)location.getY()-r,r*2,r*2);
//    }
//
//    @Override
//    public Shape getOutline(IDrawSettings settings) {
//        float margin = settings.getMargin(this, IDrawSettings.MARGIN_OUTLINE, 0);
//        float radius = settings.nucleotideRadius(this);
//        float r = radius+margin;
//        return new Ellipse2D.Float((float)location.getX()-r, (float)location.getY()-r,r*2,r*2);
//    }


    public Motif.Helix getHelix() { return Motif.Helix.getHelix(this); }
    public Motif.Loop getLoop() { return Motif.Loop.getLoop(this); }
    public Motif.MultiLoop getMultiLoop(final int minimumBranches, final boolean includeExterior) {
        return Motif.MultiLoop.getFor(this, minimumBranches, includeExterior);
    }
    public Motif.Domain getDomain() { return Motif.Domain.getDomain(this); }
    public Motif.Branch getBranch() { return Motif.Branch.getBranch(this); }
    /**
     * Gets the overall index of this nucleotide, given the starting indices for each strand as specified in strandIndexes
     * For example:
     *    strandIndexes = int[] { 1, 20, 100 }
     *    for a Nuc with index=0 in the first strand, indexInScene would return 1.
     *    for a Nuc with index=5 in the second strand, indexInScene would return 25.
     *    for a Nuc with index=2 in the third strand, indexInScene would return 102.
     */
    public int indexInScene(int[] strandIndexes) {
        return strandIndexes[strand.getIndex()] + indexInStrand();
    }

    /** Implementation of {@link INucGroup#getBases()} */
    @Override public Set<Nuc> getBases() {
        return Collections.singleton(this);
    }
    /** Implementation of {@link INucGroup#transformAll(AffineTransform)}  */
    @Override public void transformAll(final AffineTransform tr) { transform(tr); }

    public void transform(final AffineTransform tr) { tr.transform(location, location); }
    public int indexInStrand() {
        return strand==null?-1:strand.sceneToStrandIndex(indexInScene);
    }
    /**
     * Gets the overall index of this nucleotide, including all nucleotides in previous strands.
     */
    public int indexInScene() {
        //if (strandIndex == 0) return index;
        //return scene.getNucCount(strandIndex) + index;
        return indexInScene;
    }

    public String toString(final String format) {
        return format.replace("$s", symbol)
                .replace("$N", Integer.toString(indexInScene+1))
                .replace("$n", Integer.toString(indexInScene))
                .replace("$S", Integer.toString(getStrandIndex()+1))
                .replace("$s", Integer.toString(getStrandIndex()))
                .replace("$i", Integer.toString(indexInStrand()))
                .replace("$I", Integer.toString(indexInStrand()+1));
    }
    public Nuc nextInScene() { return nextInScene(1, false); }
    public Nuc nextInScene(final int steps, final boolean loopToBeginning) {
        return strand.getNextInScene(this, steps, loopToBeginning);
    }
}
