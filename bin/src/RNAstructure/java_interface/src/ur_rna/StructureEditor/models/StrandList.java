package ur_rna.StructureEditor.models;

import java.util.*;

/**
 * A list of {@link ur_rna.StructureEditor.models.Strand} that coordinates strand divisions etc to ensure
 * the segmentation of nucleotides is consistent.
 */
public class StrandList extends AbstractList<Strand> {
    private final RnaScene scene;
    private List<Nuc> nucs, readOnlyNucs; // nucs is the primary list. all strands are sub-lists of this master list.
    /** Returns a read-only list of all nucleotides. */
    public List<Nuc> allNucs() { return readOnlyNucs; }
    public RnaScene getScene() { return scene; }

    int size; // start with a single strand. this will be the case most of the time.
    Strand s1, s2;
    List<Strand> others = null;

    public StrandList(RnaScene scene, List<Nuc> nucs) {
        this.scene = scene;
        this.nucs = nucs;
        this.readOnlyNucs = Collections.unmodifiableList(nucs);
        s1 = newStrand(0);
        size = 1;
    }
    //public void setNucs(List<Nuc> nucs) { this.nucs = nucs; }

//    public void reindex() {
//        for (int i = 0; i < s1.size; i++)
//            nucs.get(i).reIndex(scene, 0, i, i);
//        if (size == 1) return;
//        int total = s1.getEnd();
//        if (size > 1)
//            for (int i = 0; i < s2.size; i++)
//                nucs.get(total).reIndex(scene, 1, i, total++);
//        if (size > 2)
//            for (int j = 0; j <; j++) {
//
//            }
//    }

    private Strand newStrand(int index) { return newStrand(index, 0, 0); }
    private Strand newStrand(int index, int offset, int size) { return new Strand(scene, this, index, nucs, offset, size); }
//    {
//        return new InternalStrand(list, index, offset, size);
//    }

//    public Iterator<Nuc> nucIterator(boolean includeNullBetweenStrands) {
//        if (!includeNullBetweenStrands) return nucs.iterator();
//        return new Iterator<Nuc>() {
//            E strand = s1;
//            int nextStrand = strand.getEnd();
//            int nucID, endID = nucs.size();
//            boolean showingNull;
//
//            @Override
//            public boolean hasNext() {
//                return nucID < endID;
//            }
//            @Override
//            public Nuc next() {
//                // When we reach the end of a strand, do NOT return the current nucleotide. Instead return null. Then on the next iteration, return the current nucleotide (which is the first one in the second strand)
//                Nuc n = nucs.get(nucID);
//                if (showingNull)
//                    showingNull = false;
//                else if (++nucID == nextStrand) {
//                    showingNull = true;
//                    if (hasNext()) {
//                        strand = get(strand.getIndex() + 1);
//                        nextStrand = strand.getEnd();
//                    } else
//                        strand = null;
//                    return null;
//                }
//                return n;
//            }
//            @Override
//            public void remove() {
//                if (showingNull)
//                    throw new NoSuchElementException("Cannot remove null because it isn't in the actual list.");
//                nucs.remove(nucID);
//                endID--;
//            }
//        };
//    }

    //<editor-fold desc="List<Strand> Implementation">
    private UnsupportedOperationException noChangeError() { return new UnsupportedOperationException("Strands cannot be added or removed. Use the divide methods instead."); }
    @Override
    @Deprecated
    public boolean add(final Strand strand) { throw noChangeError(); }
    @Override
    @Deprecated
    public boolean remove(final Object o) { throw noChangeError(); }
    @Override
    @Deprecated
    public boolean addAll(final Collection<? extends Strand> c) { throw noChangeError(); }
    @Override
    @Deprecated
    public boolean addAll(final int index, final Collection<? extends Strand> c) { throw noChangeError(); }
    @Override
    @Deprecated
    public boolean removeAll(final Collection<?> c) { throw noChangeError(); }
    @Override
    @Deprecated
    public boolean retainAll(final Collection<?> c) { throw noChangeError(); }
    @Override
    @Deprecated
    public void clear() { throw noChangeError(); }
    @Override
    @Deprecated
    public Strand set(final int index, final Strand element) { throw noChangeError(); }
    @Override
    @Deprecated
    public void add(final int index, final Strand element) { throw noChangeError(); }
    @Override
    @Deprecated
    public Strand remove(final int index) { throw noChangeError(); }

    private void rangeCheck(int index) {
        if (index < 0 || index >= size)
            throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
    }

    private String outOfBoundsMsg(int index) {
        return "Invalid Strand index: " + index + (size==0?
                ". The list is empty.":
                ", Valid indices are " + 0 + " to " + (size - 1) + "."
        );
    }

    @Override
    public int size() {
        return size;
    }
    @Override
    public boolean isEmpty() {
        return false;
    }
    @Override
    public boolean contains(final Object o) {
        return o != null && s1 == o || s2 == o || others != null && others.contains(o);
    }
    @Override
    public Object[] toArray() {
        return size == 1 ? new Object[]{s1} :
                size == 2 ? new Object[]{s1, s2} :
                        toArray(new Object[size]);
    }
    @Override
    public Strand get(final int index) {
        rangeCheck(index);
        return index == 0 ? s1 : index == 1 ? s2 : others.get(index - 2);
    }

    @Override
    public int indexOf(final Object o) {
        return o == null ? -1 : o == s1 ? 0 : o == s2 ? 1 : others == null ? -1 : others.indexOf(o) + 2;
    }
    @Override
    public int lastIndexOf(final Object o) {
        return indexOf(o); // there will never be a duplicate.
    }
    //</editor-fold>

    public Strand first() { return s1; }
    public Strand second() { return s2; }
    public Strand last() { return size == 1 ? s1 : size == 2 ? s2 : others.get(size - 3); }
    public Strand add() { return add(size); }
    public Strand add(int insertAtIndex) {
        Strand ret;
        if (insertAtIndex == 0) {
            if (size > 1) others().add(0, s2);
            s2 = s1;
            ret = s1 = newStrand(0);
        } else if (insertAtIndex == 1) {
            if (size > 1) others().add(0, s2);
            ret = s2 = newStrand(1, s1.size, 0);
        } else {
            int offset = get(insertAtIndex - 1).getEnd();
            others().add(insertAtIndex - 2, ret = newStrand(insertAtIndex, offset, 0));
        }
        size++;
        for (int i = 0; i < size; i++)
            get(i).index = i;
        return ret;
    }
    public Strand removeDirect(int index) {
        Strand ret = get(index);
        if (index == 0) {
            if (size > 1)
                s1 = s2;
            s2 = size > 2 ? others.remove(0) : null;
        } else if (index == 1) {
            s2 = size > 2 ? others.remove(0) : null;
        } else {
            others().remove(index - 2);
        }
        size--;
        for (int i = 0; i < size; i++)
            get(i).index = i;
        return ret;
    }

    private List<Strand> others() {
        if (others == null) others = new ArrayList<>(4);
        return others;
    }

    /**
     * Divides one strand into two separate strands.
     * The division occurs immediately before nucIndexInStrand, so after the division,
     * the Nuc at nucIndexInStrand-1 will be the last Nuc in Strand, while the Nuc previously at
     * nucIndexInStrand will be the first Nuc in the newly created subsequent strand.
     * */
    public void divideStrand(Strand strand, int nucIndexInStrand) {
        int index = indexOf(strand);
        if (index == -1) throw new IllegalArgumentException("Strand does not belong to this list.");
        divideStrand(index, nucIndexInStrand);
    }
    /**
     * Divides one strand into two separate strands.
     * The division occurs immediately before nucIndexInStrand, so after the division,
     * the Nuc at nucIndexInStrand-1 will be the last Nuc in Strand, while the Nuc previously at
     * nucIndexInStrand will be the first Nuc in the newly created subsequent strand.
     * */
    public void divideStrand(int strandIndex, int nucIndexInStrand) {
        Strand strand1 = get(strandIndex);
        Strand strand2 = add(strandIndex + 1);
        strand2.offset = strand1.offset + nucIndexInStrand;
        strand2.setEnd(strand1.getEnd());
        strand1.setEnd(strand2.offset);
        reIndexNucs(strand2.offset);
    }

    public void joinStrands(Strand s1, Strand s2) {
        int distance = Math.abs(s1.index-s2.index);
        int first = Math.min(s1.index, s2.index);
        if (distance==0)
            throw new IllegalArgumentException("Cannot join a strand to itself.");
        if (distance!=1)
            throw new IllegalArgumentException("Only sequential strands can be joined."); // TODO: Allow re-arrangements.
        // make sure s1 < s2
        s1 = get(first); s2 = get(first+1);
        s1.size+=s2.size;
        s2.size=0;
        removeDirect(s2.index);
        reIndexNucs(s2.offset);
    }

    static final int NO_REINDEX = Integer.MAX_VALUE;
    private int suppressUpdateCounter = 0;
    private int lowestUpdateNucIndex = NO_REINDEX;
    public void suspendUpdate() { suppressUpdateCounter++; }
    public void resumeUpdate() { resumeUpdate(false); }
    public void resumeUpdate(boolean performUpdate) {
        suppressUpdateCounter--;
        if (performUpdate && suppressUpdateCounter == 0)
            updateStrands(0);
    }
    public void updateStrands(int updateNucIndex) {
        lowestUpdateNucIndex = Math.min(lowestUpdateNucIndex, updateNucIndex);
        if (suppressUpdateCounter != 0) return;

        int last = s1.size;
        if (size > 1) {
            s2.offset = last;
            last += s2.size;
        }
        if (size > 2)
            for (int i = 0; i < size -2; i++) {
                Strand s = others.get(i);
                s.offset = last;
                last += s.size;
            }
        last().setEnd(nucs.size());

        reIndexNucs(Math.max(0, lowestUpdateNucIndex));
        lowestUpdateNucIndex = NO_REINDEX;
    }
    private void reIndexNucs(int start) {
        if (start >= nucs.size()) return;
        for (int s = 0; s < size; s++) {
            Strand strand = get(s);
            if (start < strand.size)
                for (int i = start < 0 ? 0 : start; i < strand.size; i++) {
                    Nuc n = strand.get(i);
                    n.strand = strand;
                    n.indexInScene = strand.offset+i;
                }
            start -= strand.size;
        }
    }



    //    private boolean updateAfter(Supplier<Boolean> operation) {
//        boolean ret;
//        try {
//            suspendUpdate();
//            ret = operation.get();
//        } catch (Exception e) {
//            resumeUpdate();
//            throw e;
//        }
//        resumeUpdate();
//        if (ret) updateStrands();
//        return ret;
//    }


}
