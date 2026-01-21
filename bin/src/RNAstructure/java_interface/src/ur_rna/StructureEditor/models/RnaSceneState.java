package ur_rna.StructureEditor.models;

/**
 * Represents a compressed version of an RnaScene that can be stored quickly and
 * restored to form an identical RnaScene that is completely separate, in terms of mutable Object references
 * from the original.
 */
public class RnaSceneState {
    public String title;
    public NucState[][] strands;
    public Bond.BondInfo[] bonds;
    //public StyleState[] styles;

    public static class NucState {
        public float X;
        public float Y;
        public String symbol;
        public NucStyle style;
        public int number;
    }

//    public Store getState() {
//        Store st = new Store(symbol, location);
//        st.style = this.styleIndex;
//        st.number = this.number;
//        return st;
//    }
//    public static Nuc fromState(Store st) {
//        Nuc n = new Nuc(st.symbol, st.locationX, st.locationY);
//        n.number = st.number; //historic numbering
//        return n;
//    }
//    public void reIndex(final RnaScene scene, Strand strand, final int indexInScene) {
//        this.indexInScene = indexInScene;
//        this.scene = scene;
//    }
//public static class Store {
//    public float locationX;
//    public float locationY;
//    public String symbol;
//    public int style;
//    public int number;
//    public Store(final String symbol, final Point2D.Float location) {
//        this.symbol = symbol;
//        if (location != null) {
//            this.locationX = location.x;
//            this.locationY = location.y;
//        }
//    }
//}

//    public Store getState() {
//        Store st = new Store();
//        st.N1 = N1;
//        st.N2 = N2;
//        st.S1 = S1;
//        st.S2 = S2;
//        st.type = type;
//        st.settings = style;
//        return st;
//    }
//    public static Bond fromState(final Bond.Store st) {
//        Bond b = new Bond();
//        b.N1 = st.N1;
//        b.N2 = st.N2;
//        b.S1 = st.S1;
//        b.S2 = st.S2;
//        b.type = st.type;
//        b.style = st.settings;
//        return b;
//    }
//    public static class Store {
//        public int S1, N1, S2, N2;
//        public BondType type;
//        public int settings;
//    }

//    public static class StyleState {
//        public Color fillColor;
//        public Color textColor;
//        public Color outlineColor;
//        public Font font;
//    }

//    /**
//     * Represents a trimmed-down version of an RnaScene useful for serialization
//     */
//    public static class Store {
//        public Strand.Store[] strands;
//        public Bond.Store[] bonds;
//    }
//
//    public RnaScene.Store getState() {
//        Store st = new Store();
//        st.strands = new Strand.Store[strands.size()];
//        st.bonds = new Bond.Store[bonds.size()];
//        int i = 0;
//        for(Strand s : strands)
//            st.strands[i++] = s.getState();
//        i = 0;
//        for(Bond b : bonds)
//            st.bonds[i++] = b.getState();
//        return st;
//    }
//
//
//    public void updatePairIndices() {
//        for (Strand s : strands)
//            for (Nuc n : s)
//                n.pairIndex = -1;
//        for (Bond b : bonds) {
//            b.getNuc5().pairIndex = b.index;
//            b.getNuc3().pairIndex = b.index;
//        }
//    }
//
//    public RnaScene clone() {
//        Store s = getState();
//        RnaScene copy = new RnaScene(); //sceneIndex
//        copy.setState(s);
//        return copy;
//    }
////    class MyList<T> extends ArrayList<T> {
////        public void setSize(int newSize) {
////            int size = size();
////            if (size < newSize) {
////                ensureCapacity(newSize);
////                while (size++ < newSize)
////                    add(null);
////            }
////
////            removeRange();
////        }
////    }
//
//    public void setState(Store st) {
//        bonds = new ArrayList<>(st.bonds.length);
//        for(int i = 0; i < st.bonds.length; i++)
//            bonds.add(Bond.fromState(st.bonds[i]));
//
//        strands = new ArrayList<>(st.strands.length);
//        for(int i = 0; i < st.strands.length; i++)
//            strands.add(Strand.fromState(st.strands[i]));
//
//        reIndex();
//    }
}
