package ur_rna.StructureEditor.models;

import ur_rna.StructureEditor.services.History;
import ur_rna.StructureEditor.services.SceneColorizer;
import ur_rna.StructureEditor.services.SceneDrawMode;

import java.util.*;

/**
 * Encompases all the information needed to an RNA "scene" which can include
 * one or more strands of RNA in a specific 2D layout.
 */
public class RnaScene implements Cloneable {
    public String title;
    public SceneDrawMode drawMode = SceneDrawMode.Standard;
    /** Whether this nucleotides were drawn flipped (Clockwise) or normal (Counter-Clockwise) */
    public boolean drawFlipped;

    private List<Nuc> nucs = new ArrayList<>();
    public StrandList strands = new StrandList(this, nucs);
    // public List<Bond> bonds = new ArrayList<>();
    public History<HistoryState> history = new History<>();


    public RnaScene() {}
    public Nuc getNuc(final int index) {
        return nucs.get(index);
    }

    /** Returns a read-only list of all nucleotides. */
    public List<Nuc> allNucs() { return strands.allNucs(); }

    public Strand getStrand(final int strandIndex) { return strands.get(strandIndex); }
    public Strand getStrand(final int strandIndex, final boolean createIfNotExisting) {
        if (createIfNotExisting)
            while(strands.size()<=strandIndex)
                strands.add();
        return strands.get(strandIndex);
    }
//    void setStyle(final Nuc n, final NucStyle style) {
//        if (style.index == -1)
//            addStyle(style);
//        else if (styles.get(style.index)!= style)
//            throw new NoSuchElementException();
//        n.style = style;
//    }
//    void setStyle(final Nuc n, final int styleIndex) {
//        if (styleIndex == -1)
//            n.style = null;
//        else
//            n.style = styles.get(styleIndex);
//    }

    /** Returns the list of helices in the RNA scene */
    public List<Motif.Helix> getHelices() {
        return Motif.Helix.findAll(this);
    }
    public List<Motif.Loop> getLoops() {
        return Motif.Loop.findAll(this);
    }
    public List<Motif.MultiLoop> getMultiLoops(int minimumBranches, boolean includeExterior) {
        return Motif.MultiLoop.findAll(this, minimumBranches, includeExterior);
    }
    public Motif.Helix getHelix(final Nuc nuc) {
        return Motif.Helix.getHelix(nuc);
    }
    public Motif.Loop getLoop(final Nuc nuc) {
        return Motif.Loop.getLoop(nuc);
    }
    public Motif.Domain getDomain(final Nuc nuc) {
        return Motif.Domain.getDomain(nuc);
    }
    public Motif.Branch getBranch(final Nuc nuc) { return Motif.Branch.getBranch(nuc); }

//    public void reIndex() {
//        strands.reIndexNucs();
//        int i = 0;
//        for(Nuc n : nucs)
//            n.pair = null;
//        for(Bond b : bonds) {
//            b.scene = this;
//            b.index = i++;
//            b.n5.pair = b;
//            b.n3.pair = b;
//        }
//        i = 0;
//        for(NucStyle s : styles)
//            s.index = i++;
//    }

    public int getNucCount() {
        return nucs.size();
//        int count = 0;
//        for(Strand s : strands)
//            count += s.size();
//        return count;
    }

//    public int getNucCount(int strands) {
//        int count = 0;
//        for(int i = 0; i < strands; i++)
//            count += this.strands.get(i).size();
//        return count;
//    }

    public Strand addStrand() { return strands.add(); }
//    public Strand addStrand(final String title) {
//        Strand s = new Strand();
//        s.title = title;
//        s.scene = this;
//        s.index = strands.size();
//        strands.add(s);
//        return s;
//    }

    public List<Bond> getBonds() { return getBonds(nucs, true); }
    public List<Bond> getBonds(boolean includePseudo) { return getBonds(nucs, includePseudo); }
    /** Get all bonds associated with nucleotides in the specified region, which can be a Strand or a Motif etc. */
    public List<Bond> getBonds(final Iterable<Nuc> region, boolean includePseudo) {
        List<Bond> list = new ArrayList<>();
        boolean[] added = new boolean[nucs.size()];
        for(Nuc n : region)
            if (n.pair!=null&&!added[n.pair.n5.indexInScene]&&(includePseudo||n.pair.isTruePair())) { // only add the bond if it has not already been added. Use Bond.n5 to track it.
                list.add(n.pair);
                added[n.pair.n5.indexInScene] = true;
            }
        return list;
    }

//    public List<Nuc> getBases() {
//        List<Nuc> all = new ArrayList<>(getNucCount());
//        for(Strand s : strands)
//            for (Nuc n : s)
//                all.add(n);
//        return all;
//    }
    public Nuc nucAtSceneIndex(int index) {
        return nucs.get(index);
//        int length = strands.get(0).size();
//        if (index < length)
//            return strands.get(0).get(index);
//        for (int i = 1; i < strands.size(); i++) {
//            index -= length;
//            length = strands.get(i).size();
//            if (index < length)
//                return strands.get(i).get(index);
//        }
//        return null;
    }
    public Bond addBond(Bond.BondInfo bond) { return addBond(bond.n1, bond.n2, bond.type); }
    public Bond addBond(final Nuc nuc1, final Nuc nuc2) { return addBond(nuc1, nuc2, BondType.Default); }
    public Bond addBond(final Nuc nuc1, final Nuc nuc2, final BondType type) { return addBond(new Bond(nuc1, nuc2, type)); }
    public Bond addBond(final int nuc1, final int nuc2, final BondType type) { return addBond(new Bond(nucs.get(nuc1), nucs.get(nuc2), type)); }
    private Bond addBond(Bond b) {
//        b.index = bonds.size();
//        b.scene = this;

        if (b.n5 == null || b.n3 == null) throw new IllegalStateException("The bond does not list valid nucleotides.");
        if (b.n5.pair != null || b.n3.pair != null)
            throw new IllegalStateException("A nucleotide referenced by a bond is already paired.");
        b.n5.pair = b;
        b.n3.pair = b;

        //bonds.add(b);
        return b;
    }

    public void breakBonds(final int ... nucIndices) {
        for (int i : nucIndices) {
            Nuc n = nucAtSceneIndex(i);
            if (n.isPaired())
                breakBond(n.pair);
        }
    }
    public void breakBond(Bond b) {
        // bonds.remove(b);
        b.n5.pair = null;
        b.n3.pair = null;
    }

//    public NucStyle addStyle() {
//        NucStyle s = new NucStyle();
//        s.index = styles.size();
//        styles.add(s);
//        return s;
//    }
//    public void addStyle(NucStyle style) {
//        style.index = styles.size();
//        styles.add(style);
//    }
    /** Creates a lookup-index where element[I] is the starting nuclotide index of the nucleotide belonging to strand[I].
     * For example: 3 strands with size 10, 5, and 8:
     *  {@code buildStrandIndex(0, 0) => int[] { 0, 10, 15 }} .. i.e. nucs: { S1{0-9}, S2{10-14}, S3{15-22}  }
     *  {@code buildStrandIndex(0, 3) => int[] { 0, 13, 21 }} .. i.e. nucs: { S1{0-9}, S2{13-17}, S3{21-28}  }
     *  {@code buildStrandIndex(1, 3) => int[] { 1, 14, 22 }} .. i.e. nucs: { S1{1-10}, S2{14-18}, S3{22-29}  }
     *
     * @param start The desired index of the first nucleotide in the first strand.
     * @param gapBetweenStrands The number of index positions to skip between strands (i.e. if inter-molecular linkers need to be inserted)
     * @return An integer array with each element[I] equal to the overall position of the first nucleotide in strand[I]
     */

    public int[] buildStrandIndex(int start, int gapBetweenStrands) {
        int[] index = new int[strands.size()];
        index[0] = start;
        for (int i = 1; i < index.length; i++) {
            index[i] = index[i-1] + strands.get(i-1).size()+gapBetweenStrands;
        }
        return index;
    }

    public void divideStrand(final Nuc n) {
        strands.divideStrand(n.strand, n.indexInStrand());
    }
    public void joinStrands(final Strand s1, final Strand s2) {
        strands.joinStrands(s1, s2);
    }


    public RnaSceneState getState() {
        RnaSceneState state = new RnaSceneState();
        state.strands = new RnaSceneState.NucState[strands.size()][];
        for (int si = 0; si < strands.size; si++) {
            Strand s = strands.get(si);
            RnaSceneState.NucState[] ss = state.strands[si] = new RnaSceneState.NucState[s.size()];
            for(int ni=0; ni < s.size; ni++) {
                Nuc n = s.get(ni);
                RnaSceneState.NucState ns = ss[ni] = new RnaSceneState.NucState();
                ns.number = n.number;
                ns.style = n.style==null?null:n.style.clone();
                ns.symbol = n.symbol;
                ns.X = n.location.x;
                ns.Y = n.location.y;
            }
        }
        List<Bond> bonds = getBonds();
        state.bonds = new Bond.BondInfo[bonds.size()];
        for (int i = 0; i < bonds.size(); i++) {
            Bond b = bonds.get(i);
            state.bonds[i] = new Bond.BondInfo(
                    b.n5.indexInScene,
                    b.n3.indexInScene,
                    b.type);
        }
        return state;
    }
    public void loadState(final RnaSceneState state) {
        int total = 0;
        for (int si = 0; si < state.strands.length; si++)
            total += state.strands[si].length;
        nucs = new ArrayList<>(total);
        strands = new StrandList(this, nucs);
        strands.suspendUpdate();
        for (int si = 0; si < state.strands.length; si++) {
            RnaSceneState.NucState[] ss = state.strands[si];
            Strand s = getStrand(si, true);
            for(int ni=0; ni < ss.length; ni++) {
                RnaSceneState.NucState ns = ss[ni];
                Nuc n = s.add(ns.symbol);
                n.number = ns.number;
                n.style = ns.style==null?null:ns.style.clone();
                n.location.x = ns.X;
                n.location.y = ns.Y;
            }
        }
        strands.resumeUpdate(true);
        for (Bond.BondInfo b : state.bonds)
            addBond(b);
    }

    public void clearBonds() {
        for (Nuc n : nucs)
            n.pair = null;
    }
    public boolean hasBonds() {
        for (Nuc n : nucs)
            if (n.pair != null)
                return true;
        return false;
    }
    public Strand firstStrand() {
        return strands.get(0);
    }
    public Strand secondStrand() {
        return strands.size>1?strands.get(1):null;
    }
    public void removeNucs(final Nuc[] nucs) {
        Arrays.sort(nucs);
        for (int i = nucs.length - 1; i >=0; i--) {
            if (nucs[i].isPaired())
                breakBond(nucs[i].pair);
            nucs[i].getStrand().remove(nucs[i]);
        }
    }
    public void colorize(final SceneColorizer c) {
        c.color(allNucs());
    }
    public void copyBonds(final RnaScene copyFrom) {
        for(Bond b : copyFrom.getBonds())
            addBond(b.getNuc5().indexInScene, b.getNuc3().indexInScene, b.type);
    }
    public void identifyPseudoknots() {
        Set<Bond> crossing = new HashSet<>();
        List<Bond> bonds = getBonds(false);
        Motif.findNonCrossingBonds(bonds, crossing);
        for(Bond b : bonds)
            if (crossing.contains(b))
                b.type = BondType.Pseudo;
    }
}
