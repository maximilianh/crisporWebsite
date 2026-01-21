package ur_rna.StructureEditor.models;

import ur_rna.Utilities.Strings;

import java.util.*;

import static ur_rna.StructureEditor.models.StrandList.NO_REINDEX;

/**
 * A Strand represents a single RNA strand.
 * An RNA drawing could potentially contain multiple strands, and
 * the strands can each contain inter-molecular connections (e.g. base pairs)
 * to other strands.
  */
public class Strand extends AbstractList<Nuc> implements INucGroup {
    private RnaScene scene;
    // Package-private constructor. Only to be created by a StrandList.
    Strand(RnaScene scene, StrandList owner, final int index, final List<Nuc> nucs, final int offset, final int size) {
        this.owner = owner;
        this.index = index;
        this.source = nucs;
        this.offset = offset;
        this.size = size;
        this.scene = scene;
    }

    public RnaScene getScene() { return scene; }

    /** Get the index (relative to all nucleotides in the scene) of the first nucleotide in this Strand */
    public int getSceneStart() { return offset; }
    /** Get the index (relative to all nucleotides in the scene) immediately AFTER last nucleotide in this Strand */
    public int getSceneEnd() { return offset + size; }

    public int strandToSceneIndex(int nucIndexInStrand) {
        return nucIndexInStrand + offset;
    }
    public int sceneToStrandIndex(int nucIndexInScene) {
        return nucIndexInScene - offset;
    }

    /** Implements {@link INucGroup#getBases()} */
    public Collection<Nuc> getBases() { return this; }

    public Nuc getNext(final Nuc nuc) { return getNext(nuc, 1); }
    public Nuc getPrev(final Nuc nuc) { return getNext(nuc, -1); }
    /** Returns the next nucleotide after the specified one,
     * regardless of whether it is in this strand or another one.
     * If the resulting index is greater than the number of nucleotides,
     * null will be returned, unless loopToBeginning is true, in which case
     * the scene is considered to be a circuit and numbering
     * continues at the beginning (i.e. the next Nuc after the last Nuc is the first Nuc etc).
     * */
    public Nuc getNextInScene(final Nuc nuc, int steps, boolean loopToBeginning) {
        int next = nuc.indexInScene+steps;
        final int size = source.size();
        if (next < 0 || next >= size) {
            if (!loopToBeginning) return null;
            next %= size; if (next < 0) next += size; // get modulus
        }
        return source.get(next);
    }
    public Nuc getNext(final Nuc nuc, int steps) {
        int next = nuc.indexInScene + steps;
        return next >= offset && next < offset + size ? source.get(next) : null;
    }
    public Nuc add(final String symbol) {
        Nuc n = new Nuc(symbol);
        add(n);
        return n;
    }
    public Nuc add() { return add(""); }

    protected final StrandList owner;
    protected final List<Nuc> source;
    int index; // index of this strand in parent List.
    int offset;
    int size;

    /** Get the index of this strand in the list of RnaScene strands. */
    public int getIndex() { return index; }
    /** Get the scene index of the last nucleotide in this strand plus one. I.e. the index just beyond this strand: {@code offset+size} */
    int getEnd() { return offset + size; }
    void setEnd(final int end) { this.size = end - offset; }

    public Nuc first() { return size == 0 ? null : source.get(offset); }
    public Nuc last() { return size == 0 ? null : source.get(offset + size - 1); }

    //<editor-fold desc="List<Strand> Implementation">
    protected void rangeCheck(int index) {
        if (index < 0 || index >= size)
            throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
    }

    protected void rangeCheckForAdd(int index) {
        if (index < 0 || index > size)
            throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
    }

    protected String outOfBoundsMsg(int index) {
        return "Invalid Nucleotide index: " + index + ", Valid indices are " + offset + " to " + (offset + size - 1) + ".";
    }

    @Override
    public Nuc set(int index, Nuc element) {
        rangeCheck(index);
        return source.set(index + offset, element);
    }
    private UnsupportedOperationException unsupported() {
        return new UnsupportedOperationException("This modification is not permitted on strands. Use the RnaScene or the StrandList to modify the nucleotides in a Strand.");
    }

    @Override
    public Nuc get(int index) {
        rangeCheck(index);
        return source.get(index + offset);
    }

    @Override
    public int size() {
        return size;
    }

    @Override
    public void add(int index, Nuc element) {
        rangeCheckForAdd(index);
        source.add(index + offset, element);
        element.indexInScene = offset + index;
        element.strand = this;
        size++;
        owner.updateStrands(offset + index + 1);
    }

    @Override
    public Nuc remove(final int index) {
        rangeCheck(index);
        Nuc n = source.remove(offset + index);
        size--;
        owner.updateStrands(offset+index);
        return n;
    }

    @Override
    public boolean removeAll(final Collection<?> c) {
        boolean ret;
        owner.suspendUpdate(); // don't get updates for each individual remove.
        try {
            ret = super.removeAll(c); // super checks to see whether each item is contained in THIS Strand, then calls remove, so size is already decremented.
        } finally {
            owner.resumeUpdate();
        }
        if (ret) owner.updateStrands(NO_REINDEX);
        return ret;
    }
    @Override
    public boolean retainAll(final Collection<?> c) {
        boolean ret;
        owner.suspendUpdate(); // don't get updates for each individual remove.
        try {
            ret = super.retainAll(c); // super checks to see whether each item is contained in THIS Strand, then calls THIS.remove, so size is already decremented.
        } finally {
            owner.resumeUpdate();
        }
        if (ret) owner.updateStrands(NO_REINDEX);
        return ret;
    }

    @Override
    public void clear() {
        if (size == 0) return;
        owner.suspendUpdate(); // don't get updates for each individual remove.
        try {
            super.clear(); // super checks to see whether each item is contained in THIS Strand, then calls THIS.remove, so size is already decremented.
        } finally {
            owner.resumeUpdate();
        }
        owner.updateStrands(offset);
    }

    //    public void reIndex(int strandID) {
    //        int nucID = offset;
    //        for (int i = 0; i < size; i++)
    //            nucs.get(nucID).reIndex(scene, this, nucID++);
    //    }

    public List<Nuc> addPlaceholders(int index, int count) { return addAll(index, Strings.fromChar('\0', count)); }
    public List<Nuc> addPlaceholders(int count) { return addPlaceholders(size, count); }
    public List<Nuc> addAll(CharSequence sequence) { return addAll(size, sequence); }
    public List<Nuc> addAll(int index, CharSequence sequence) {
        Nuc[] all = new Nuc[sequence.length()];
        for (int i = 0; i < all.length; i++)
            all[i] = new Nuc(Character.toString(sequence.charAt(i)));
        List<Nuc> list = Arrays.asList(all);
        addAll(index, list);
        return list;
    }

    @Override
    public boolean addAll(Collection<? extends Nuc> c) {
        return addAll(size, c);
    }

    @Override
    public boolean addAll(int index, Collection<? extends Nuc> c) {
        rangeCheckForAdd(index);
        int cSize = c.size();
        if (cSize == 0)
            return false;
        owner.suspendUpdate();
        try {
            source.addAll(offset + index, c); // NB: addAll does NOT call add !!
        } finally {
            owner.resumeUpdate();
        }
        size += cSize;
        owner.updateStrands(offset + index); // need to update all these nucs.
        return true;
    }

    //<editor-fold desc="Iterator<Nuc> and ListIterator<Nuc>">
    @Override
    public Iterator<Nuc> iterator() {
        return listIterator(0);
    }

    @Override
    public ListIterator<Nuc> listIterator(final int startIndex) {
        rangeCheckForAdd(startIndex);
        return new ListIterator<Nuc>() {
            private int index = startIndex;
            private int lastRet = -1; // stores the last index used in a call to next() or previous(). Used to remove the last item returned. reset to -1 to indicate remove has been called.

            public boolean hasNext() {
                return index < size;
            }
            public boolean hasPrevious() {
                return index > 0;
            }

            public Nuc next() {
                if (hasNext())
                    return source.get(lastRet = (offset + index++));
                else
                    throw new NoSuchElementException();
            }

            public Nuc previous() {
                if (hasPrevious())
                    return source.get(lastRet = (offset + --index));
                else
                    throw new NoSuchElementException();
            }
            public int nextIndex() {
                return index;
            }
            public int previousIndex() { return index - 1; }
            public void remove() {
                if (lastRet == -1)
                    throw new IllegalStateException();
                Strand.this.remove(lastRet);
                if (lastRet < index)
                    index--;
                lastRet = -1;
            }
            public void set(Nuc e) {
                if (lastRet == -1)
                    throw new IllegalStateException();
                Strand.this.set(lastRet, e);
            }
            public void add(Nuc e) {
                Strand.this.add(index, e);
                lastRet = -1;
                index++; // skip past the added item.
            }
        };
    }
    //</editor-fold>
    //</editor-fold>
}
