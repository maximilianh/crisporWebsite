package ur_rna.StructureEditor.models;

import ur_rna.Utilities.annotation.NotNull;
import ur_rna.Utilities.annotation.Nullable;
import ur_rna.Utilities.geom.Vec2D;

import java.awt.geom.Point2D;
import java.util.NoSuchElementException;

/**
 * Represents a bond between two nucleotides.
 */
public class Bond implements Comparable<Bond>, Cloneable  {
    public RnaScene getScene() { return n5.getScene(); }
    //transient int index;


    /**  The index of this bond within the collection of bonds.    */

    /**  S1 and S2 represent the strand indices, N1 and N2 represent the nucleotide indices within each respective strand.  */
    //public int S1, N1, S2, N2;
    public BondType type;
    //public int style;
    public Nuc n5, n3;

    //public int getIndex() { return index; }
    public Nuc getNuc5() {
        return n5; //scene.strands.get(S1).get(N1);
    }
    public Nuc getNuc3() {
        return n3; //scene.strands.get(S2).get(N2);
    }

//    public Nuc getNuc5() {
//        return isOrdered() ? getNuc5() : getNuc3();
//    }
//    public Nuc getNuc3() {
//        return isOrdered() ? getNuc3() : getNuc5();
//    }

    public Bond clone() {
        try {
            return (Bond) super.clone();
        } catch (CloneNotSupportedException ex) {
            throw new InternalError(); //should never happen
        }
    }
    /**
     * Given a nucleotide on one side of the bond, this returns the nucleotide on the other side.
     * If the given nucleotide is not on either side, this returns null.
     */
    public @Nullable Nuc getOther(Nuc n) {
        if (n == n5) return n3;
        if (n == n3) return n5;
        throw new NoSuchElementException("The specified nucleotide is not associated with this bond.");
//
//        if (n.strandIndex == S1 && n.index == N1)
//            return getNuc3();
//        if (n.strandIndex == S2 && n.index == N2)
//            return getNuc5();
//        return null;
    }

    public Motif.Helix getHelix() {
        return n5.getHelix();
    }

//    public Bond getNext() { return getNext(1); }
//    public Bond getPrev() { return getNext(-1); }
//    public Bond getNext(int step) {
//        step = index + step;
//        if (step < 0 || step >= scene.bonds.size())
//            return null;
//        return scene.bonds.get(step);
//    }

    /**
     * Whether this bond represents an inter-molecular bond (i.e. whether the bond connects nucleotides from different RNA strands)
     * @return False if both nucleotides belong to the same strand. True otherwise (i.e. the second nucleotide belongs to a different strand than the first).
     */
    public boolean isHybrid() { return n5.strand != n3.strand; /*S1 != S2;*/  }

 //   public Bond() {}
    public Bond(Nuc n1, Nuc n2) { this(n1, n2, BondType.Default); }
    public Bond(final Nuc n1, final Nuc n2, final BondType type) {
        //this.scene = n5.strand.getScene();
        this.n5 = n1;
        this.n3 = n2;
//        S1 = n5.strandIndex;
//        N1 = n5.index;
//        S2 = n3.strandIndex;
//        N2 = n3.index;
//        scene = n5.scene;
        this.type = type;
        makeOrdered();
    }
//    public Bond(final RnaScene scene, int n5, int n3, BondType type) {
//        this.scene = scene;
//        this.n5 = scene.getNuc(n5);
//        this.n3 = scene.getNuc(n3);
//        this.type = type;
//        makeOrdered();
//    }
//    public Bond(final RnaScene scene, BondInfo info) {
//        this(scene, info.n5, info.n3, info.type);
//    }

    @NotNull
    public BondType getType() {
        if (type==null) type = BondType.Default;
        return type;
    }

    public void unpair() { getScene().breakBond(this); }

    /**
     * Assume we are traveling around a multi-branch loop from 5' to 3'.
     * When we reach the "entrance" branch (aka "Trunk Helix") we encounter a bond that is
     * reversed compared to all other bonds in the loop (which all belong to "exit" branches
     * (aka "Limb Helices")
     * In order to keep traversing the multi-branch loop (in the same direction)
     * we need to go from the 3' side of the bond to the 5' side.
     *
     * In order to simplify the logic of traversing multi-branch loops,
     * we define the "right" nucleotide of a Bond to be
     *      - null: if Bond is internal (i.e. not at the start or end of its helix)
     *      - n3: if Bond is the first bond in its helix (i.e.: Bond == Bond.getHelix().getPair(0))
     *      - n5: if Bond is the last bond in its helix (i.e.: Bond == Bond.getHelix().getPair(0))
     * The "left" nucleotide is converse: n5 and n3 are switched.
     *
     * Using this definition, a multi-branch loop can be traversed by
     * calling {@link Nuc#getNext()} on unpaired nucleotides and
     * jumping across branches by calling {@link Bond#right(Nuc)} ()} on paired Nucs.
     * (or alternatively by calling {@link Nuc#getPrev()} on unpaired nucleotides and
     * {@link Bond#left(Nuc)} ()}.
     *
     * @param ref Any other nucleotide that is in the same MultiLoop as this Nuc.
     *            This is used for reference in cases where it is ambiguous which
     *            loop is being traversed. This ONLY occurs when this nucleotide is
     *            in a SINGLE-BP helix, therefore its basepair is both the FIRST and
     *            the LAST basepair.
     * */
    public Nuc right(Nuc ref) {
        switch (getMultiLoopRole(ref)) {
            case Entrance: return n5;
            case Exit: return n3;
            default: return null;
        }
    }
    public Nuc left(Nuc ref) {
        switch (getMultiLoopRole(ref)) {
            case Entrance: return n3;
            case Exit: return n5;
            default: return null;
        }
    }

    public enum MultiLoopRole {
        /** This bond is an internal one. It is not part of any Multi-Branch, Internal or Hairpin Loop */
        None,
        /** This bond is the "entrance" to a multi-branch loop. I.e. it is the last bond of the "Trunk" Helix. */
        Entrance,
        /** This bond is one of the "exits" out of a multi-branch loop. I.e. it is the first bond of a "Limb" Helix. */
        Exit,
        /** The role of this bond could not be determined with the provided information. */
        Unknown
    }

    /**
     * Determine whether this bond is the "entrance" ("Trunk") to a Multi-Branch Loop by passing in
     * another reference nucleotide that is in the same loop.
     * @param ref Any other nucleotide that is in the same MultiLoop as this Nuc.
     *            This is used for reference in cases where it is ambiguous which
     *            loop is being traversed. This ONLY occurs when this nucleotide is
     *            in a SINGLE-BP helix, therefore its basepair is both the FIRST and
     *            the LAST basepair.
     * @return True if this is an entrance/trunk base-pair or false if this is an exit/limb pair.
     */
    public MultiLoopRole getMultiLoopRole(@Nullable Nuc ref) {
        boolean first = isFirst();
        boolean last = isLast();
        if (!first && !last)
            return MultiLoopRole.None; // there is a helix bond before and after this one, so it not in any loop.
        if (!first)
            return MultiLoopRole.Entrance; // there is a helix bond before and after this one, so it not in any loop.
        if (!last) return MultiLoopRole.Exit; // there is a helix bond before and after this one, so it not in any loop.
        if (ref == null) return MultiLoopRole.Unknown;
        if (ref.indexInScene == n5.indexInScene || ref.indexInScene == n3.indexInScene)
            throw new IllegalArgumentException("The reference nucleotide cannot be one of the bases in this pair.");
        return (ref.indexInScene > n5.indexInScene && ref.indexInScene < n3.indexInScene) ?
                MultiLoopRole.Entrance : MultiLoopRole.Exit;
    }
    public Bond nextInHelix() { return nextInHelix(1); }
    public Bond prevInHelix() { return nextInHelix(-1); }
    public Bond nextInHelix(int direction) {
        // Verify that the next nuc after the 5' nuc is paired
        // and that its pair is the same as the previous neighbor of the 3' nuc.
        Nuc n; return (null!=(n= n5.nextInScene(direction, false))
            && (null != (n = n.getPaired(false)))
            && n == n3.nextInScene(-direction, false))
                ? n.getPairBond() : null;
    }
    public boolean isFirst() { return prevInHelix() == null; }
    public boolean isLast() { return nextInHelix() == null; }
    /**
     * Returns the midpoint between the two nucleotides in the bond.
     *  If  assignTo  is NOT null, its location is set to the midpoint
     *  and it is returned by the function.
     *  Otherwise, a new Vec2D is created and returned.
     */
    public Point2D midpoint(@Nullable Point2D assignTo) {
        if (assignTo == null) return midpoint();
        assignTo.setLocation((n5.location.x + n3.location.x)/2, (n5.location.y + n3.location.y)/2);
        return assignTo;
    }
    public Vec2D midpoint() { return Vec2D.getMidpoint(n5.location, n3.location); }
    /**
     * Returns a Vec2D that represents the direction normal to the bond -- i.e.
     * orthogonal to the vector from n5 to n3.
     * If the structure was drawn counter-clockwise (the default),
     * the normal vector of the first bond in a helix will point
     * along the direction of the helix. If the structure was drawn clockwise,
     * the normal vector will point in the opposite direction
     * (along the helix, but away from it).
     *
     *  If  assignTo  is NOT null, its coordinates are set to the normal,
     *  and it is returned by the function.
     *  Otherwise, a new Vec2D is created and returned.
     */
    public Vec2D normal(@Nullable Vec2D assignTo) {
        if (assignTo == null)
            assignTo = new Vec2D(n5.location, n3.location);
        else
            assignTo.setLocation(n3.location.x-n5.location.x, n3.location.y - n5.location.y);
        return assignTo.rotate90();
    }
    public Vec2D normal() { return normal(null); }

    /** Returns true for any bond that is not a PseudoKnot bond or
     * a Prohibited-Constraint Bond.
     */
    public boolean isTruePair() {
        switch (type) {
            case Prohibited:
            case Pseudo:
                return false;
            default: return true;
        }
    }

    /**
     * A class for storing the indices of bonded nucleotides instead of direct references.
    *  Sometimes used during initial construction or deserialization of RnaScenes
    */
    public static class BondInfo {
        public int n1, n2;
        public BondType type;
        public BondInfo(final int n1, final int n2, final BondType type) {
            this.n1 = n1;
            this.n2 = n2;
            this.type = type;
        }
    }

//    public Bond(final int s1, final int n5, final int s2, final int n3) {
//        S1 = s1;
//        N1 = n5;
//        S2 = s2;
//        N2 = n3;
//        type = BondType.Default;
//    }

    public boolean contains(Nuc test) {
        return n5 ==test|| n3 ==test;
//        return S1 == test.strandIndex && N1 == test.index ||
//                S2 ==test.strandIndex && N2 == test.index;
    }

    public void makeOrdered() {
        if (!isOrdered())
            swap();
    }

//    /** If the nucleotides are on the same strand, return the nucleotide closest to the 5' end.
//     * Otherwise, return the nucleotide with the lowest strand-index
//     * (i.e. strand 1 nucleotide before strand 2 nucleotide) regardless of the
//     * nucleotide positions within each respective strand.
//      */
//    public Nuc get5End() {
//        return isOrdered() ? N1 : N2;
//    }
//    public Nuc get3End() {
//        return isOrdered() ? N2 : N1;
//    }

    /**
     * Swap the order of the nucleotides. (e.g. to put then in canonical order.
     */
    private void swap() {
        Nuc tmp = n3;
        n3 = n5;
        n5 = tmp;
//        int tmp = S2;
//        S2 = S1; S1 = tmp;
//        tmp = N2;
//        N2 = N1; N1 = tmp;
    }
    /**
     * Returns true if the nucleotides are in canonical order.
     *
     * Canonical order is hereby defined as:
     * Nucleotides are .
     * All nucleotides from strand 1 are listed listed in increasing order by index
     * (this corresponds to listing them from 5' to 3'). If there are subsequent strands,
     * then the nucleotides from each strand are listed after those of previous strands.
     * (E.g. all strand 1 nucleotides come before strand 2 nucleotides.)
     *
     * The following function returns true if N1 comes before N2 in the above ordering.
     * I.e. {@code N1.strand < N2.strand} or {@code N1.strand == N2.strand && N1.index <= N2.index}.
     * The function also returns true if both nucleotides are null or if N2 is null (but N1 is not).
     */
    public boolean isOrdered() {
//        if (N1 == null)
//            return N2 == null;
        return getNuc5().compareTo(getNuc3()) <= 0;
    }

    public int compareTo(Bond other) {
        int result = this.getNuc5().compareTo(other.getNuc5());
        if (result == 0)
            return this.getNuc3().compareTo(other.getNuc3());
        return result;
    }
    public static int compare(Bond lhs, Bond rhs) {
        return lhs.compareTo(rhs);
    }

    @Override
    public String toString() {
        return "Bond " + n5 + ":" + n3 + (type==null||type==BondType.Default?"":" ("+type+")");
    }
}
