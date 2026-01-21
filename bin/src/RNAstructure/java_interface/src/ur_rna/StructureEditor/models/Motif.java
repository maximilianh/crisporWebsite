package ur_rna.StructureEditor.models;

import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.geom.Vec2D;

import java.util.*;

/**
 * Represents a structural element (motif) of an RNA structure -- e.g. a base, pair, bulge, loop, helix, hairpin, domain, etc.
 */
public abstract class Motif {
    //public static void transform(IMotif g) {

//    /**
//     * Determine if any part of the motif is under the specified point.
//     * (Note that the motif may allow a small margin around each physical element).
//     *
//     * @param point The location to test.
//     * @param settings The settings, which specify values such as nucleotide radius that
//     *                 vary depending on the graphics scale and other settings.
//     * @return True if any part of the motif is under the specified point, or false otherwise.
//     */
//    @Override
//    public boolean hitTest(final Point2D point, final IDrawSettings settings) {
//        // Perform a faster test to exclude points outside the bounds.
//        Collection<Nuc> bases = getBases();
//        if (!getBounds(bases, settings).contains(point)) return false;
//
//        for (Nuc n : bases)
//            if (n.hitTest(point, settings))
//                return true;
//        return false;
//    }
//    /**
//     * Gets a Shape representing the outline of the motif which can be used by the drawing
//     * framework to draw a selection indicator etc.
//     *
//     * @param settings The settings, which specify values such as nucleotide radius that
//     *                 vary depending on the graphics scale and other settings.
//     * @return A Shape object representing the outline of the motif.
//     */
//    @Override
//    public Shape getOutline(final IDrawSettings settings) {
//        return null;
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
//        return getBounds(getBases(), settings);
//    }
//
//    /**
//     * Returns the smallest bounding rectangle that contains all of the bounding rectangles of the
//     * nucleotides in this motif.
//     */
//    protected Rectangle2D getBounds(Collection<Nuc> bases, final IDrawSettings settings) {
//        Iterator<Nuc> i = bases.iterator();
//        Rectangle2D rc = i.next().getBounds(settings);
//        while(i.hasNext())
//            rc.add(i.next().getBounds(settings));
//        return rc;
//    }
//
//    @Override
//    public void translate(float dx, float dy) {
//        for (Nuc n : getBases())
//            n.translate(dx, dy);
//    }

    /**
     * A helix represents a series of adjacent nucleotides which are
     * all paired to a another series of adjacent nucleotides on the same
     * RNA strand or a different RNA strand.
     * <p>
     * From this definition, it follows that a helix may be composed of
     * either inter- or intra-molecular pairs but never both, and that a helix
     * can include bases from at most two RNA strands.
     * <p>
     * Due to this definition, the only information required to store a helix unambiguously
     * are the first and last bases of the series on either strand. Let these be
     * called n5 and n3 (for the 5' and 3' bases respectively).
     * It follows that n5 must belong to same strand as n3.
     * Similarly n5.pairedTo must belong to the same strand as n3.pairedTo,
     * but this is not necessarily the same strand that n5 and n3 belong to.
     */
    public static class Helix implements IMotif {
        /** The 5'-most nucleotide in the helix.  */
        Nuc nuc5;
        /** The number of basepairs in the helix. */
        int pairCount;

        public RnaScene getScene() { return nuc5.getScene(); }

        protected Helix(final Nuc nuc5, int pairCount) {
            this.nuc5 = nuc5;
            this.pairCount = pairCount;
        }
        //protected Helix(Bond start, Bond end) { this(start.getScene(), get Math.min(start.index, end.index), Math.abs(start.index - end.index) + 1); }

        /**
         * Get the number of base-pairs in the helix (which is half the number of nucleotides).
         *
         * @return The number of base-pairs in the helix.
         */
        public int size() {
            return pairCount;
        }
        public Bond getPair(int index) {
            if (index < 0 || index >= pairCount)
                throw new IndexOutOfBoundsException();
            return nuc5.getNext(index).getPairBond();
        }

        public Bond first() { return getPair(0); }
        public Bond last() { return getPair(pairCount-1); }

        public boolean isClockwise() {
            Vec2D helixNorm = new Vec2D(first().midpoint(), last().midpoint());
            return first().normal().dot(helixNorm) < 0;
        }

        /**
         * Whether or not the helix contains inter-molecular pairs
         * (i.e. a hybridization of two distinct strands).
         *
         * @return true if n5.strand is different than n3.strand and false otherwise.
         */
        public boolean isHybrid() {
            return getPair(0).isHybrid();
//            Bond b = getPair(0);
//            return b.isHybrid(); // .S1 != b.S2;
        }

        @Override
        public Collection<Nuc> getBases() {
            List<Nuc> list = new ArrayList<>(pairCount * 2);
            for (int i = 0; i < pairCount; i++) {
                Bond b = getPair(i);
                list.add(b.getNuc5());
                list.add(b.getNuc3());
            }
            return list;
        }

        /**
         * Find all distinct helices in specified list of strands.
         *
         * @return a list of helices.
         */
        public static List<Helix> findAll(RnaScene scene) {
            List<Helix> found = new ArrayList<>();
            Bond first = null, last = null;
            ObjTools.RefInt orientation = new ObjTools.RefInt(); // whether the second strand is reversed relative to the first.

            // assume bonds are sorted and nucleotides are ordered.
            /* Each bond connects nuc1 (on strand 1) to nuc2 (on strand 2)
               for a bond "B" to be in a helix with the previous bond "A",
               the following has to be true:
               1) B.nuc1 has to be ADJACENT to LAST.nuc1
                  ADJACENT means that B.nuc1 has to be on the same strand as A.nuc1 and,
                  the index of B.nuc1 must be immediately before or after A.nuc1.
                  Due to the way bonds are ordered, B.nuc1.index is always >= A.nuc1.index
                  (assuming they are on the same strand). Therefore B.nuc1.index must be equal to
                  A.nuc1.index + 1,
               2) B.nuc2 must be ADJACENT to LAST.nuc2
                  Again, this requires that they are on the same strand,
                  but strand 2 can be oriented parallel or anti-parallel (aka "reversed") relative
                  to strand 1. This means that B.nuc2 can be either before or after A.nuc2
                  (i.e. B.nuc2.index can be either A.nuc2.index+1 or A.nuc2.index-1).
               3) Once a second pair has been added to a helix, the orientation of a potential third
                  bond "C" must follow suit. E.g. if B.nuc2 is one nucleotide before A.nuc2, then C.nuc2 must
                  be one nucleotide before B.nuc2. Conversely if B.nuc2 is one nucleotide after A.nuc2,
                  then C.nuc2 must be one nucleotide after B.nuc2.
             */
            for (Bond current : scene.getBonds(false)) {
                if (first != null && bondIsNextInHelix(last, current, orientation))
                    last = current; // extend the existing helix to include the current pair. (the helix is defined by "first", "last", and all pairs in-between)
                else {
                    // add the existing helix to the list if it is composed of two or more pairs.
                    if (last != first)
                        found.add(new Helix(first.n5, last.n5.indexInScene - first.n5.indexInScene + 1));
                    // initiate a new helix
                    last = first = current;
                    orientation.value = 0;
                }
            }
            // add the final helix, if there was one.
            if (last != first)
                found.add(new Helix(first.n5, last.n5.indexInScene - first.n5.indexInScene + 1));
            return found;
        }

        protected static boolean isHelixBond(Bond b) {
            return b != null && b.isTruePair();
        }

        /**
         * Return the helix that contains the given nuc or null if nuc is not in a helix (ie. it is unpaired)
         */
        public static Helix getHelix(final Nuc nuc) {
            if (!isHelixBond(nuc.pair)) return null;
            ObjTools.RefInt orientation = new ObjTools.RefInt();
            Nuc n, start = nuc.pair.getNuc5(), end = start;

            // Find 5' end
            while ((n = start.getPrev()) != null && isHelixBond(n.pair) && bondIsNextInHelix(n.pair, start.pair, orientation))
                start = n;
            while ((n = end.getNext()) != null && isHelixBond(n.pair) && bondIsNextInHelix(end.pair, n.pair, orientation))
                end = n;
            return new Helix(start, end.indexInScene - start.indexInScene+1);
        }
        /**
         * Determines whether the "current" bond is immediately after the "prev" bond in a helix.
         * Specifically:
         * 1. current.type must be a valid helix bond (base pair)
         * 2. The "from" nucleotides (current.nuc1 and prev.nuc1) must be on the same strand.
         * 3. The "to" nucleotides (current.nuc2 and prev.nuc2) must be on the same strand (but not necessarily the same strand as nuc1)
         * 4. The "from" nucleotide of current must be positioned immediately after that of "prev". (i.e. current.nuc1.index == prev.nuc1.index + 1)
         * 5. The "to" nucleotide of current must be positioned immediately after OR before that of "prev". (i.e. current.nuc2.index == prev.nuc1.index +/- 1)
         * (whether it is +1 or -1 will depend on the orientation of the second in the helix strand relative to the first.)
         * 6. If the orientation of the second strand to the first is known (which will always be the case when at least two adjacent bonds have been identified)
         * then the orientation between this bond and the last must match the one previously determined (e.g. from the preceding two bonds).
         */
        private static boolean bondIsNextInHelix(Bond prev, Bond current, ObjTools.RefInt orientation) {
            // Verify that the "nuc1" nucleotides of this bond and the last one are on the same strand.
            //    and that the "nuc2" nucleotides are on the same strand as each other (not necessarily the same as nuc1)
            if (current.n5.strand != prev.n5.strand || current.n3.strand != prev.n3.strand) return false;
            // Verify that current.nuc1 is exactly 1 position after ref.nuc1 (i.e. current.nuc1.index == ref.nuc1.index + 1)
            if (prev.n5.getNext() != current.n5) return false;
            // Verify that current.nuc1 is exactly 1 position before or after ref.nuc1 (i.e. current.nuc2.index == ref.nuc2.index +/- 1)
            int pos2 = current.n3.indexInScene - prev.n3.indexInScene; // should equal orientation, unless orientation==0, in which case it should be either 1 or -1.

            // The orientation can be either 1 (parallel) or -1 (anti-parallel).
            if (orientation.value == 0 && (pos2 == 1 || pos2 == -1)) {
                // The orientation was NOT previously established (because prior to this only ONE known bond existed)
                // But now that we have a second bond, it determines the orientation.
                orientation.value = pos2;
            }
            // If the orientation was previously known, the value determined here must match it exactly.
            // (if the orientation was set in the block above, this is redundant, but that is not the common case).
            return pos2 == orientation.value;
        }

//        /**
//         * An alternate method to get a helix. This traverses a strand nucleotide by nucleotide.
//         * It is slower than the bond approach, but perhaps easier to understand and requires less code.
//         */
//        public Motif.Helix getHelixAlt(final Nuc nuc) {
//            if (nuc.pair == null)
//                return null;
//            // `start` and `end` represent the nucleotides on the same side of the helix at the 5' and 3' ends respectively.
//            // n is the subsequent nucleotide and p is its pair.  pp is the nuc paired to the previous value of `n`.
//            Nuc n, p, pp, start = nuc, end = nuc;
//            int orientation = 0;
//
//            pp = start.getPaired();
//            while ((n = start.getPrev()) != null && (p = n.getPaired()) != null) {
//                if (orientation == 0)
//                    orientation = p == pp.getNext() ? 1 : p == pp.getPrev() ? -1 : 0;
//                if (p != (pp = pp.getNext(-orientation)))
//                    break;
//                start = n;
//            }
//            pp = end.getPaired();
//            while ((n = end.getNext()) != null && (p = n.getPaired()) != null) {
//                if (orientation == 0)
//                    orientation = p == pp.getNext() ? -1 : p == pp.getPrev() ? 1 : 0;
//                if (p != (pp = pp.getNext(-orientation)))
//                    break;
//                end = n;
//            }
//            return new Helix(start.getPairBond(), end.getPairBond());
//        }

        public boolean isParallel() {
            if (pairCount == 1) return false; //unable to determine
            Bond b1 = getPair(0), b2 = getPair(1);
            return b1.n5.indexInScene - b2.n5.indexInScene == b1.n3.indexInScene - b2.n3.indexInScene;
        }
        public Bond lastPair() {
            return getPair(pairCount-1);
        }
    }

    /**
     * A Loop represents a series of adjacent unpaired nucleotides.
     * Loops can be categorized as follows:
     * <dl>
     *     <dt>Hairpin</dt>
     *     <dd>A single-stranded region that connects one side of a helix to the other without any intervening helices.
     *     e.g.:  A:B-C-D-E-F-(A)
     *           (Notation:   X:Y represents a base-pair;  X-Y represents a 3'-to-5' covalent bond between nucleotides X and Y.
     *            (A) means the same "A" nucleotide that appeared earlier in the sequence.)
     *     </dd>
     *     <dt>Internal Loop</dt>
     *     <dd>A single-stranded region that connects one side of a helix to the same side of another helix. (Same side means without any intervening helices.
     *     e.g.:  A:B-C-D-E-F-G:H-[I-J-...]-(A)  --There can be no intervening helices between H and A
     *     </dd>
     *     <dt>Multibranch Loop</dt>
     *     <dd>A series of single-stranded regions that connect multiple helices together.
     *     e.g.:  A:B-C-D-E-F-G:H-I-J-K:L-M-N-(A)
     *     </dd>
     *     <dt>Terminal Loop</dt>
     *     <dd>A single-stranded regions that is only connected to a helix on one side.
     *     e.g.:  A-B-C-D-E-F-G:H
     *     </dd>
     * </dl>
     */
    public static class Loop extends Segment implements IMotif {
        private LoopType type;
        public LoopType getType() {
            if (type == null)
                type = calcLoopType(this);
            return type;
        }

        public enum LoopType {
            /** The nucleotides at both ends of the loop are paired to each other. */
            Hairpin,
            /** Traversal from the end of the loop and back to the start of the loop encounters only two helices (i.e. one at the start of the loop and one at the end). */
            Internal,
            /** Traversal from the end of the loop and back to the start of the loop encounters more than two helices. */
            Multibranch,
            /** The loop is part of a pseudoknot */
            Pseudo,
            /** The loop is exterior loop -- traversal from the end does not connect back to the start. (i.e. this loop is not included in any domain) */
            Exterior,
            /** The loop is an {@link LoopType#Exterior} one at the start or end of the strand. */
            Terminal
        }
        public Loop() {}
        // public Loop(final Strand s, final int start, final int end) {super(s, start, end); }
        public Loop(final Nuc start, final Nuc end) { super(start, end); }
        // public Loop(final Nuc start, final int length) {  super(start, length);  }
        public Loop(final Nuc start, final Nuc end, LoopType type) {
            super(start, end);
            this.type = type;
        }

        /**
         * Return the single-stranded region (loop) that contains the given nuc or null if nuc is not in a loop (i.e. it paired)
         */
        public static Loop getLoop(final Nuc nuc) {
            if (nuc.isPaired(false))
                return null;
            Nuc n, n5 = nuc, n3 = nuc;
            while ((n = n3.getNext()) != null && !n.isPaired(false))
                n3 = n;
            while ((n = n5.getPrev()) != null && !n.isPaired(false))
                n5 = n;
            return new Loop(n5, n3);
        }

        /**
         * Find all distinct single-stranded regions in specified list of strands.
         *
         * @return a list of single-stranded regions (aka loops).
         */
        public static List<Loop> findAll(RnaScene scene) {
            List<Loop> found = new ArrayList<>();
            int start;

            for (Strand s : scene.strands) {
                start = -1;
                for (int i = 0; i < s.size(); i++) {
                    if (s.get(i).isPaired(false)) {
                        if (start != -1) {
                            found.add(new Loop(s.get(start), s.get(i - 1)));
                            start = -1;
                        }
                    } else {
                        if (start == -1)
                            start = i;
                    }
                }
                if (start != -1)
                    found.add(new Loop(s.get(start), s.get(s.size() - 1)));
            }
//            for (Loop l : found)
//                l.type = Loop.calcLoopType(l);

            return found;
        }

        /**
         * Retrieve a list of all bonds that are part of a multi-branch loop.
         * (i.e. the bases of all helices that extend from the multi-branch loop.
         * The list is built by starting with n and traversing the strand in a forward direction,
         * jumping across base-pairs until either n is reached again or the strand terminates.
         * If the strand terminates before n is reached, the n nucleotide was NOT in
         * a multi-branch loop, and NULL is returned.
         *
         * @param n                          Any nucleotide in the multi-branch loop.
         *                                   Note that if n is paired, it is technically possible for it to be simultaneously in two
         *                                   adjacent multi-branch loops (assuming it belongs to a "helix" containing only a single base-pair).
         *                                   For this reason, it is better to pass an unpaired nucleotide in for n, to avoid any ambiguity.
         * @param includeExteriorLoop        if false, null is returned if the given nucleotide is part of an exterior loop.
         *                                   Otherwise all the branches starting at or after n will be included in the returned list
         *                                   (but those BEFORE n will not be.)
         * @param basePairRole Only applies when n is paired and the bases in either direction are unpaired (or paired
         *                     to bases in different helices than n).
         *                     (i.e. n is in a helix comprised of a single pair).
         *                     This parameter imfluences how the basepair is interpreted -- as an exit or entrance pair.
         * @return a list of bonds (if nStart was in fact in a multi-branch loop) or null otherwise.
         */
        public static List<Bond> getMultiBranchBasePairs(Nuc n, boolean includeExteriorLoop, Bond.MultiLoopRole basePairRole) {
            List<Bond> bonds = new ArrayList<>();
            if (n.isPaired(false)) {
                // determine which side of the bond is in the loop.
                Bond.MultiLoopRole calculatedRole = n.pair.getMultiLoopRole(null);
                if (calculatedRole == Bond.MultiLoopRole.Unknown) calculatedRole = basePairRole;
                switch (calculatedRole) {
                    case Entrance:
                        // This causes the loop to be transversed from n's  5' side (n's 3' pair is the LAST bond).
                        // This assumes that n is the 5'-most pair of a branch that is a "Trunk" entering the MB-loop.
                        n = n.getPaired();
                        break;
                    case Exit:
                        // Leaving n as-is will cause the loop to be transversed from n's  3' side (n's 5' pair is the LAST bond).
                        // This assumes that n's pair is the 5'-most pair of a branch that is a "Limb" exiting from the MB-loop
                        break;
                    case None:
                        throw new IllegalArgumentException("The nucleotide passed to getMultiBranchBasePairs is internal to a helix -- not part of a branched loop.");
                    default:
                        // We could not determine
                        // Assume this is an entrance branch.
                        n = n.getPaired();
                        break;
                }
            }
            Nuc nStart = n;
            do {
                if (n.isPaired(false)) { // prevent pseudoknot bonds
                    n = n.getPaired();
                    if (n == nStart)
                        return bonds;
                    bonds.add(n.getPairBond());
                }
                n = n.getNext();
            } while (n != null && n != nStart);
            return (n == nStart || includeExteriorLoop) ? bonds : null;
        }

        private static LoopType calcLoopType(final Loop l) {
            Nuc start = l.getNuc5().getPrev(); // get the first paired base IMMEDIATELY BEFORE the loop (or null if the strand starts at the beginning of the loop)
            if (start == null) return LoopType.Terminal;

            Nuc end = l.getNuc3(); // get last nuc in the loop
            int helices = 0; //number of intervening helices

            //System.out.println("calc loop: ");
            while (true) {
                //System.out.println("end: " + end);
                if (end == null)
                    return helices == 0 ? LoopType.Terminal : LoopType.Exterior;
                Nuc pair = end.getPaired(false);
                //System.out.println("pair: " + pair);
                //check for pseudoknot
                // check for differing strand.

                if (pair == null) {
                    end = end.getNext();
                } else if (pair.equals(start)) {
                    //System.out.println("reached start");
                    // we made it back to the base-pair at the start of the loop.
                    return helices == 0 ? LoopType.Hairpin : helices == 1 ? LoopType.Internal : LoopType.Multibranch;
                } else {
                    helices++;
                    end = pair.getNext();
                }
            }
        }
        public static List<Bond> getMultiBranchBasePairs(Nuc n) { return getMultiBranchBasePairs(n, false, Bond.MultiLoopRole.Entrance); }
        public static List<Bond> getMultiBranchBasePairs(Nuc n, boolean includeExteriorLoop) { return getMultiBranchBasePairs(n, includeExteriorLoop, Bond.MultiLoopRole.Entrance); }

        //        /** Returns true if the bond should end the series of nucleotides that compose a loop. */
//        protected static boolean isLoopEndingBond(Bond b) {
//            return b.type == BondType.Pair || b.type == BondType.Wobble;
//        }
        @Override
        public String toString() {
            return (type==null?"Loop ":type + " Loop ") + getNuc5() + ":" + getNuc3();
        }
    }

    public static class MultiLoop implements IMotif {
        public final List<Object> elements = new ArrayList<>();
        private Bond[] branchBonds;
        private Loop[] loops;
        private int loopNucs;

        protected MultiLoop() { }
        public MultiLoop(Nuc loopNuc) {
            if (loopNuc.isPaired()) throw new IllegalArgumentException("The nucleotide passed to a MultiLoop cannot be paired.");
            List<Bond> bonds = Loop.getMultiBranchBasePairs(loopNuc, true);
            if (bonds == null) return;
            init(bonds, false);
        }
        public MultiLoop(final List<Bond> bonds, final boolean preSorted) {
            init(bonds, preSorted);
        }
        protected void init(final List<Bond> bonds, boolean preSorted) {
            branchBonds = bonds.toArray(new Bond[bonds.size()]);
            if (!preSorted) Arrays.sort(branchBonds);
            int loopCount = 0;
            for(int i = 0; i < branchBonds.length; i++) {
                elements.add(branchBonds[i]);
                // first bond [0] is from the Trunk helix, so it's 5' side starts the loop.
                // in contrast, all the other bonds are from Limb helices, so their 3' sides continue the loop.
                Nuc nextNuc = i == 0 ? branchBonds[i].n5 : branchBonds[i].n3;
                Loop nextLoop = nextNuc.nextInScene(1, true).getLoop();
                if (nextLoop!=null) {
                    loopCount++;
                    elements.add(nextLoop);
                    loopNucs += nextLoop.size();
                }
            }
            loops = new Loop[loopCount];
            loopCount = 0; // used for counting in for-loop
            for (Object o : elements)
                if (o instanceof Loop)
                    loops[loopCount++] = (Loop)o;
        }
        @Override
        public Collection<Nuc> getBases() {
            Collection<Nuc> bases = new ArrayList<Nuc>(nucCount());
            for(Object o : elements)
                if (o instanceof Bond) {
                    bases.add(((Bond) o).n5);
                    bases.add(((Bond) o).n3);
                } else if (o instanceof Loop)
                    bases.addAll(((Loop) o).getBases());
            return bases;
        }
        public List<Bond> getBonds() {
            return Arrays.asList(branchBonds);
        }
        public List<Loop> getLoops() {
            return Arrays.asList(loops);
        }
        /** Get the "starting" base -- i.e. the 5'-most base (the base with the lowest scene-index) */
        public Nuc getN5() { return branchBonds[0].n5; }
        /** Get the "ending" base -- i.e. the 3'-most base (the base with the highest scene-index) */
        public Nuc getN3() { return branchBonds[0].n3; }
        public Bond getTrunk() { return branchBonds[0]; }

        public static List<MultiLoop> findAll(final RnaScene scene, final int minimumBranches, boolean includeExteriorLoop) {
            List<MultiLoop> list = new ArrayList<>();
            List<Bond> bonds = Loop.getMultiBranchBasePairs(scene.getNuc(0), true);
            for(Bond b : bonds) {
                if (includeExteriorLoop && bonds.size()>=minimumBranches)
                    list.add(new MultiLoop(bonds, true));
                findAll(list, b.getHelix().lastPair(), minimumBranches);
            }
            return list;
        }
        protected static void findAll(List<MultiLoop> list, Bond trunkBond, final int minimumBranches) {
            List<Bond> bonds = Loop.getMultiBranchBasePairs(trunkBond.n5);
            if (bonds.size()>=minimumBranches)
                list.add(new MultiLoop(bonds, true));
            // skip the first bond, which we know already is the far-end of the trunkBond
            // TODO: Optimize this by moving the getMultiBranchBasePairs to this class and
            // write an overload the accepts a reusable array of bonds and add a parameter to skip the trunk bond
            for(int i = 1; i < bonds.size(); i++)
                findAll(list, bonds.get(i).getHelix().lastPair(), minimumBranches);
        }
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("MultiLoop @ ").append(getN5()).append("[").append(elements.size()).append("]");
            int i = 0;
            for(Object o : elements)
                sb.append("  ").append(i++).append(":").append(o).append(";");
            return sb.toString();
        }
        public static MultiLoop getFor(final Nuc nuc, final int minimumBranches, final boolean includeExterior) {
            List<Bond> bonds = Loop.getMultiBranchBasePairs(nuc, includeExterior);
            if (bonds == null || bonds.size() < minimumBranches) return null;
            return new MultiLoop(bonds, false);
        }
        public int helixCount() {
            return branchBonds.length;
        }
        public int countLoops() {
            return loops.length;
        }
        public int loopNucCount() {
            return loopNucs;
        }
        public int nucCount() {
            return branchBonds.length*2+loopNucs;
        }
    }

    /**
     * Represents a segment of an RNA strand that extends from a nucleotide to its pair.
     */
    public static class Domain extends Segment implements IMotif {
        public Domain(){}
        public Domain(final Nuc n) {
            super(n, validatePair(n));
            orderNucs();
        }
        public Domain(final Bond b) { this(b.getNuc5());}
        private static Nuc validatePair(Nuc n) {
            Nuc p = n.getPaired();
            if (p == null) throw new IllegalArgumentException("A Domain can only be defined by a paired nucleotide.");
            if (p.strand != n.strand) throw new IllegalArgumentException("A Domain can only be defined by two nucleotide on the same strand.");
            return p;
        }

        /**
         * Gets the domain that starts at the specified nucleotide (and ends with its paired nucleotide)
         * or null if the nucleotide is unpaired or is paired with a nucleotide on a different strand.
         */
        public static Domain getDomain(final Nuc nuc) {
            Bond b = nuc.getPairBond();
            if (b == null || b.getNuc5().strand != b.getNuc3().strand)
                return null;
            return new Motif.Domain(b);
        }
    }

    /**
     * Represents a complete helix (the "base" of the Branch) along with all helix and loop regions on one side of it.
     * A more precise definition depends on whether or not the base helix is uni-molecular or bi-molecular
     * (i.e. composed of one or two strands).
     *
     * For a uni-molecular base helix, the corresponding Branch is the same as the Domain defined
     *     by the 5'-most nucleotide in the base helix. (This differs from a Domain in general because a Domain
     *     need not include all nucleotides in a helix, while a Branch by definition does.)
     *     Specifically, a uni-molecular Branch with base helix 'H' is the union of
     *     {@code getDomain(nuc)}  with  {@code getHelix(nuc)} for any nucleotide 'nuc' in 'H'.
     *
     * For a bimolecular base helix the definition becomes more complicated, and the choice of which side of the
     * helix to include becomes arbitrary, leading to two possible branches, both with the same base:
     *{@code
     *                      'N5' 'N3'
     *                        |'H'|
     *      S1         5'.....GCAUG.....3'
     *      S2     (left).....CGUAC.....(right)        Note: S2 can be aligned parallel or anti-parallel to S1
     *                        |   |                          so "left" can be 5' for parallel or 3' for anti-parallel.
     *                     'N5P' 'N3P'
     *      B1                --------->
     *      B2           <---------
     *}
     * Let 'H' be the base helix, 'S1' the strand with the lower index, and 'S2' the other strand.
     *
     * Let 'N5' be the 5'-most nucleotide on S1 that is in H and 'N3' the 3'-most nucleotide on 'S1' that is in 'H'.
     * Let 'N5P' be the nucleotide paired with N5 and N3P the one paired with N3.
     * Note that N5P and N3P are both on S2, but the alignment of S1 with respect to S1 could be anti-parallel
     * (which places N3P closer to the 5' end of S2 than N5P) or parallel (in which case N5P is on 5' side of N3P).
     *
     * The two possible branches based at H are:
     *    B1: All nucleotides in H and all nucleotides on the 3' side of N3 on S1 and to the right* of N3P on S2.
     *    B2: All nucleotides in H and all nucleotides on the 5' side of N5 on S1 and to the left* of N5P on S2.
     *      * The terms "right" and "left" depend on the alignment of S2 with respect to S1:
     *        For S2-parallel:  "right" means the 3' side whereas "left" means the 5' side.
     *        For S2-anti-parallel:  the above definitions are reversed -- "right" means 5' and "left" means 3'.
     */
    public static class Branch implements IMotif {
        private Helix base;
        private boolean isLowerBranch;

        public Branch(){}

        /**
         * Creates a unimolecular branch that extends from the specified nucleotide to its pair.
         * @param n A nucleotide in the base helix of a unimolecular branch.
         */
        public Branch(final Nuc n) {
            base = Helix.getHelix(n);
            if (base == null)
                throw new IllegalArgumentException("A branch cannot be created from an unpaired nucleotide.");
        }
        public Branch(Helix base, boolean isLowerBranch) {
            this.base = base;
            this.isLowerBranch = isLowerBranch;
        }

        public boolean isHybrid() {
            return base.isHybrid();
        }
        public boolean isParallel() {
            return base.isParallel();
        }
        public Helix getBaseHelix() {
            return base;
        }


        public void validate() {
            // Parallel:
            //  S1:       n5-ABC.........-n3   or   n5-.........XYZ-n3
            //  S2:       n5-ABC.......-n3            n5-.......XYZ-n3
            // Anti:                                                                         (in "ordinal" view.  like-letters are paired.)
            //  S1:       n5-ABC.........-n3   or   n5-.........XYZ-n3  ==>         n5-ABC.........-n3  or   n5-.........XYZ-n3
            //  S2:       n3-ABC.......-n5            ns-.......XYZ-n5       n5-.......CBA-n3                         n5-ZYX.......-n3

            // Parallel:       s1n1:s2n1  or  s1n2:s2n2
            // Anti-Parallel:  s1n1:s2n2[-H.length+1]  or  s1n2:s2n1[-H.length+1]
        }

        /**
         * Get a set of all nucleotides included in the group.
         * @return a Collection of nucleotides.
         */
        @Override
        public Collection<Nuc> getBases() {
            if (isHybrid()) {
                Collection<Nuc> list1, list2;
                Bond start = base.getPair(0), end = base.getPair(base.pairCount - 1);
                Strand s1 = start.getNuc5().getStrand(), s2 = start.getNuc3().getStrand();
                list1 = isLowerBranch ?
                        s1.subList(0, end.getNuc5().indexInStrand() + 1) :
                        s1.subList(start.getNuc5().indexInStrand(), s1.size());
                list2 = isLowerBranch == isParallel() ?
                        s2.subList(0, end.getNuc3().indexInStrand() + 1) :
                        s2.subList(start.getNuc3().indexInStrand(), s2.size());
                List<Nuc> all = new ArrayList<>(list1.size()+list2.size());
                all.addAll(list1);
                all.addAll(list2);
                return all;
            } else {
                Nuc n5 = base.getPair(0).getNuc5();
                return n5.getStrand().subList(n5.indexInStrand(), n5.getPaired().indexInStrand() + 1);
            }
        }

        /**
         * Return the smallest branch that contains nuc or null if nuc is
         * not contained in any branch (e.g. if it is in a terminal loop or part of a multibranch loop segment.)
         * @param nuc
         */
        public static Branch getBranch(final Nuc nuc) {
            // If the nuc is paired, it's own Helix defines the base of the branch.
            if (nuc.isPaired())
                return new Branch(nuc);

            // If the nuc is unpaired, search around the (hairpin, internal, or multibranch) loop until we find a
            // basepair that encloses the nucleotides in the loop.
            // If we reach the end of the strand, the nuc was in a terminal loop and no enclosing branch can exist, so return null.
            Nuc n = nuc.getPrev(); // search in reverse (arbitrary--fwd is just as good)
            while(n != null) {
                if (n.isPaired()) {
                    Bond p = n.pair;
                    if (p.n5.indexInScene < nuc.indexInScene && p.n3.indexInScene > nuc.indexInScene)
                        return new Branch(n); // This basepair encloses the loop, so return the branch defined by the helix that contains it.

                    // Otherwise, jump over the basepair
                    n = n.getPaired();
                }
                n = n.getPrev();
            }
            return null; // the nuc was in a terminal loop. (note that either 5' or 3' terminus will end here, even though we are searching in reverse, because any intervening bonds will be "jumped" over.
        }
    }

    /**
     * Represents a contiguous section of an RNA strand.
     * Note that in a Segment, the Nucs 'start' and 'end' are not necessarily in order (5' to 3' respectively)
     * (However for derived types, such as a Loop or Branch, start is always the 5'-most Nuc and end is always the 3'-most Nuc.)
     */
    public static class Segment implements INucGroup {
        public Segment(){}
        public Segment(Strand s, int start, int end) {
            this.start = s.get(start);
            this.end = s.get(end);
        }
        public Segment(Nuc start, Nuc end) {
            if (start.strand != end.strand)
                throw new IllegalArgumentException("The two nucleotides in a segment must both be from the same strand.");
            this.start = start;
            this.end = end;
        }
        public Segment(Nuc start, int length) {
            this(start, start.getNext(length-1));
        }

        public void orderNucs() {
            if (start.indexInScene > end.indexInScene)
                swapNucs();
        }
        public void swapNucs() {
            Nuc swap = start;
            start = end;
            end = swap;
        }

        /** The start of the segment (the 5' side) */
        public Nuc start;
        /** The end of the segment (the 3' side) */
        public Nuc end;

//        public int getStrand() {
//            return start.strandIndex;
//        }
//        public int getPos5() {
//            return Math.min(start.index, end.index);
//        }
//        public int getPos3() {
//            return Math.max(start.index, end.index);
//        }
        /** Get the number of nucleotides in the loop. */
        public int size() {
            return Math.abs(end.indexInScene - start.indexInScene) + 1;
        }
        public Nuc getNuc5() {
            return start.indexInScene <= end.indexInScene ? start : end;
        }
        public Nuc getNuc3() {
            return start.indexInScene <= end.indexInScene ? end : start;
        }
        public Nuc getNucAt(final int pos) {
            return getNuc5().getNext(pos);
        }

//        public int getStrand() {
//            return strand;
//        }
//        public int getStart() {
//            return start;
//        }
//        public int getLength() {
//            return length;
//        }
//        public Nuc getStartNuc() {
//            return scene.strands.get(strand).get(start);
//        }
//        public Nuc getEndNuc() {
//            return scene.strands.get(strand).get(start+length-1);
//        }
//        public Nuc getNuc(int position) {
//            return scene.strands.get(strand).get(start+position);
//        }


        /**
         * Get a set of all nucleotides included in the group.
         * @return a Collection of nucleotides.
         */
        @Override
        public Collection<Nuc> getBases() {
            int i = start.indexInStrand(), j = end.indexInStrand();
            if (i <= j)
                return start.getStrand().subList(i, j+1);
            return start.getStrand().subList(j, i+1);
        }
    }

    public static boolean hasPseudoKnots(RnaScene scene, boolean ignoreMarkedPseudoBonds) {
        final int length = scene.getNucCount();
        Bond[] lookup = new Bond[length];
        List<Bond> bonds = scene.getBonds(!ignoreMarkedPseudoBonds);
        for (Bond b : bonds)
            lookup[b.n5.indexInScene] = b;

        IntervalStack stack = new IntervalStack(Math.min(8,length / 4));
        stack.push(0, length - 1);
        Bond b = null;
        while (stack.pop()) {
            while(stack.i <= stack.j && null == (b = lookup[stack.i])) // backtrack stores returns the 3' position or -1
                stack.i++;
            if (stack.i > stack.j || b == null)
                continue; // pop next item from stack
            // here b != null and i < j
            int k = b.n3.indexInScene;
            if (k > stack.j) return true;
            if (stack.i + 1 <= k-1)
                stack.push(stack.i+1, k-1);
            if (k + 1 <= stack.j)
                stack.push(k+1, stack.j);
        }
        return false;
    }

    public static Set<Bond> findPseudoKnots(RnaScene scene) {
        Set<Bond> pseudo = new HashSet<>();
        findNonCrossingBonds(scene.getBonds(), pseudo);
        return pseudo;
    }
    public static Set<Bond> findNonCrossingBonds(RnaScene scene) { return findNonCrossingBonds(scene.getBonds()); }
    public static Set<Bond> findNonCrossingBonds(Collection<Bond> bonds) { return findNonCrossingBonds(bonds, null); }
    public static Set<Bond> findNonCrossingBonds(Collection<Bond> bonds, Collection<Bond> crossing) {
        if (bonds.size()==0) return Collections.emptySet();
        final int length = bonds.iterator().next().getScene().getNucCount();
        Bond[] lookup = new Bond[length];
        for (Bond b : bonds)
            lookup[b.n5.indexInScene] = lookup[b.n3.indexInScene] = b;

        short[][] outer = new short[length][length];
        short[][] backtrack = new short[length][length];

        for (int n = 1; n < length; n++) {
            for (int i = 0; i < length - n; i++) {
                int j = i + n;
                outer[i][j] = outer[i + 1][j];
                backtrack[i][j] = -1;
                if (lookup[i]!=null && lookup[i].n5.indexInScene==i) {
                    int k = lookup[i].n3.indexInScene;
                    if (k > i && k <= j) {
                        int tmp = 1;
                        if (i + 1 <= k - 1) {
                            tmp += outer[i + 1][k - 1];
                        }
                        if (k + 1 <= j) {
                            tmp += outer[k + 1][j];
                        }
                        if (tmp > outer[i][j]) {
                            outer[i][j] = (short) tmp;
                            backtrack[i][j] = (short) k;
                        }
                    }
                }
            }
        }
        HashSet<Bond> nonCrossing = new HashSet<>(bonds.size());

        // Backtrace
        IntervalStack stack = new IntervalStack(Math.min(8,length / 4));
        stack.push(0, length - 1);
        short k = -1;
        while (stack.pop()) {
            while(stack.i < stack.j && -1 == (k = backtrack[stack.i][stack.j])) // backtrack stores returns the 3' position or -1
                stack.i++;
            if (stack.i >= stack.j || k == -1)
                continue; // pop next item from stack
            // here k != -1 and i < j
            nonCrossing.add(lookup[k]);
            if (stack.i + 1 < k-1)
                stack.push(stack.i+1, k-1);
            if (k + 1 < stack.j)
                stack.push(k+1, stack.j);
        }

        if (crossing != null && nonCrossing.size() < bonds.size()) {
            for (Bond b : bonds)
                if (!nonCrossing.contains(b))
                    crossing.add(b);
        }

        return nonCrossing;

//        HashMap<Integer, ArrayList<Bond>> index2BPs = new HashMap<>();
//        for (Bond b : bonds) {
//            int i = b.n5.indexInScene;
//            if (!index2BPs.containsKey(i)) {
//                index2BPs.put(i, new ArrayList<Bond>());
//            }
//            index2BPs.get(i).add(b);
//        }
//
//        short[][] tab = new short[length][length];
//        short[][] backtrack = new short[length][length];
//        int theta = 3;
//
//        for (int i = 0; i < length; i++) {
//            for (int j = i; j < Math.min(i + theta, length); j++) {
//                tab[i][j] = 0;
//                backtrack[i][j] = -1;
//            }
//        }
//        for (int n = theta; n < length; n++) {
//            for (int i = 0; i < length - n; i++) {
//                int j = i + n;
//                tab[i][j] = tab[i + 1][j];
//                backtrack[i][j] = -1;
//                if (index2BPs.containsKey(i)) {
//                    ArrayList<Bond> vi = index2BPs.get(i);
//                    for (int numBP = 0; numBP < vi.size(); numBP++) {
//                        Bond mb = vi.get(numBP);
//                        int k = mb.n3.indexInScene;
//                        if ((k != -1) && (k <= j) && (i < k)) {
//                            int tmp = 1;
//                            if (i + 1 <= k - 1) {
//                                tmp += tab[i + 1][k - 1];
//                            }
//                            if (k + 1 <= j) {
//                                tmp += tab[k + 1][j];
//                            }
//                            if (tmp > tab[i][j]) {
//                                tab[i][j] = (short) tmp;
//                                backtrack[i][j] = (short) numBP;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        HashSet<Bond> planar = new HashSet<>();
//
//        // Backtracking
//        Stack<Point> intervals = new Stack<Point>();
//        intervals.add(new Point(0, length - 1));
//        while (!intervals.empty()) {
//            Point p = intervals.pop();
//            if (p.x <= p.y) {
//                if (backtrack[p.x][p.y] == -1) {
//                    intervals.push(new Point(p.x + 1, p.y));
//                } else {
//                    int i = p.x;
//                    int j = p.y;
//                    int nb = backtrack[p.x][p.y];
//                    Bond mb = index2BPs.get(i).get(nb);
//                    int k = mb.n3.indexInScene;
//                    planar.add(mb);
//                    intervals.push(new Point(i + 1, k - 1));
//                    intervals.push(new Point(k + 1, j));
//                }
//            }
//        }

//        // Remaining base pairs
//        for (int i : index2BPs.keySet()) {
//            ArrayList<Bond> vi = index2BPs.get(i);
//            for (Bond mb : vi) {
//                if (!planar.contains(mb)) {
//                    others.add(mb);
//                }
//            }
//        }
    }
    private static class IntervalStack {
        short[] stack;

        public IntervalStack() { this(10); }
        public IntervalStack(int capacity) {
            stack = new short[capacity*2];
        }

        /** Current values */
        public short i, j;
        int pos;
        public void push(int i, int j) {
            if (stack.length < pos +2)
                stack = Arrays.copyOf(stack, stack.length * 2);
            stack[pos++] = (short)i;
            stack[pos++] = (short)j;
        }
        public boolean pop() {
            if (pos==0) return false;
            j = stack[--pos]; // pop in reverse order or push (j, then i)
            i = stack[--pos];
            return true;
        }
    }

}
