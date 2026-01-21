package ur_rna.StructureEditor.services.fileIO;

import ur_rna.RNAstructure.Nucleobase;
import ur_rna.RNAstructure.RnaBackendException;
import ur_rna.RNAstructure.backend.*;
import ur_rna.StructureEditor.FileType;
import ur_rna.StructureEditor.models.*;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.PathTools;
import ur_rna.Utilities.StopWatch;
import ur_rna.Utilities.geom.Vec2D;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.regex.Pattern;

import static ur_rna.Utilities.Strings.isEmpty;

/**
 * RnaFileIO is used to load and save RNA files (CT, SEQ etc)
 * <p>
 * This is a further abstraction of the RNA Backend class. It is a mixture of
 * backend (c++) and Java-based functions.
 * Ideally, we would go the opposite route--integrate the backend classes into
 * the java program to make the distinction as seamless as possible, but
 * unfortunately the backend c++ is too brittle, and adding the necessary features
 * to the c++ classes might have a performance impact for the command-line interfaces.
 */
public class RnaFileIO {
    public static final int MAX_SEQUENCE_LENGTH = 100000000;
    private final static AppLog log = AppLog.getDefault();
    private static String DATA_PATH;
    public final static char InterMolecularBreakChar = 'I';
    public final static String InterMolecularBreakStr = "III"; // TODO: change "I" to another character, e.g. "_"
    public final static int InterMolecularBreakLength = InterMolecularBreakStr.length();
    public final static String ALPHABET_RNA = "rna";

//    public final static int FILE_CT = 1;
//    public final static int FILE_SEQ = 2;
//    public final static int FILE_PFS = 3;
//    public final static int FILE_SAV = 4;
//    public final static int FILE_DBN = 11;

    public final static int SEQUENCE_FMT_TEXT = 0;
    public final static int SEQUENCE_FMT_SEQ = 1;
    public final static int SEQUENCE_FMT_FASTA = 2;

//    public RnaScene getScene(int structureNumber) { return getScene(structureNumber, true, true); }
//    private RnaScene getScene(int structureNumber, boolean createIfNotLoaded, boolean doLayout) {
//        RnaScene r = structures[structureNumber - 1];
//        if (createIfNotLoaded && r == null) {
//            r = structures[structureNumber - 1] = new RnaScene(structureNumber);
//            createStrand(r, structureNumber);
//            updateSceneBonds(structureNumber);
//            if (doLayout) performNucLayout(structureNumber);
//        }
//        return r;
//    }
//
//    public void writeDrawingFile(RnaSceneGroup scenes, String path, boolean append, boolean extendedCT) throws IOException, RnaBackendException {
//        DrawingFileIO parser = new DrawingFileIO();
//        BufferedWriter writer = null;
//        if (extendedCT) throw new RuntimeException("Not implemented.");
//        try {
//            writer = new BufferedWriter(new FileWriter(path, append));
//            parser.write(scenes, writer, 2, false);
//        } catch (Exception ex) {
//            throw new IOException("Error writing Drawing file: "+ex.getMessage()+"\nPath: \"" + path + "\"", ex);
//        } finally {
//            try {
//                if (writer != null)
//                    writer.close();
//            } catch (IOException ex) {
//                ex.printStackTrace();
//            }
//        }
//    }
//
//    public RnaSceneGroup readDrawingFile(String path, boolean embeddedInCT) throws IOException, RnaBackendException {
//        BufferedReader reader = null;
//        DrawingFileIO parser = new DrawingFileIO();
//        try {
//            reader = new BufferedReader(new FileReader(path));
//            String line; int lineNumber = 0;
//            parser.parseStart();
//            while (null != (line = reader.readLine())) {
//                lineNumber++;
//                if (embeddedInCT) {
//                    // trim #! from the beginning of each line.
//                    if (line.length() < 2) continue;
//                    if (line.charAt(0) != '#' || line.charAt(1) != '!') continue;
//                    parser.parseContentLine(line.substring(2), lineNumber, 3);
//                } else {
//                    // read the content directly.
//                    parser.parseContentLine(line, lineNumber, 0);
//                }
//            }
//            parser.parseEnd();
//            return parser.getParsed();
//        } catch (Exception ex) {
//            File tmp = new File(path);
//            throw new IOException(fmt("Failed to load %s \"%s\": %s\nFull Path: \"%s\"", FileType.CteDrawing.getDescription(), tmp.getName(), ex.getMessage(), getCanonicalPath(tmp)), ex);
//        } finally {
//            try {
//                if (reader != null)
//                    reader.close();
//            } catch (IOException ex) {
//                ex.printStackTrace();
//            }
//        }
//    }

    private static void checkErrors(RNA rna, String operation) throws RnaBackendException {
        int err = rna.GetErrorCode();
        if (err != 0)
            throw new RnaBackendException("Error " + operation + ": " + rna.GetFullErrorMessage(), err);
    }

    public static RNA createRNAfromScene(RnaScene scene, boolean includePairs, boolean allowUnknownBases, boolean skipThermo)
            throws RnaBackendException {
        verifyDataPath();
        RNA rna = new RNA(getRNASequence(scene), RNABackend.SEQUENCE_STRING, ALPHABET_RNA, allowUnknownBases, skipThermo);
        checkErrors(rna, "building RNA structure");
        rna.EnsureStructureCapcacity(1); // indicate that there is a single structure.
        if (scene.title != null) {
            rna.GetStructure().SetSequenceLabel(scene.title);
            setStructureTitle(rna, 1, scene.title);
        }
        //rna.AddComment(scene.title, 1);
        if (includePairs)
            updateBackendBonds(scene, rna, 1);

        return rna;
    }
    public static RNA writeRNAstructureGroup(RnaSceneGroup scenes, boolean allowUnknownBases, boolean skipThermo)
            throws RnaBackendException {
        //DEBUG_TIMING: StopWatch sw = new StopWatch(true);
        verifyDataPath();
        //DEBUG_TIMING: sw.println("verifyDataPath").restart();
        String seq = getRNASequence(scenes.get(0));
        //DEBUG_TIMING: sw.println("getRNASequence").restart();
        RNA rna = new RNA(seq, RNABackend.SEQUENCE_STRING, ALPHABET_RNA, allowUnknownBases, skipThermo);
        //DEBUG_TIMING: sw.println("RNA").restart();
        if (!isEmpty(scenes.getTitle()))
            rna.GetStructure().SetSequenceLabel(scenes.getTitle());
        rna.EnsureStructureCapcacity(scenes.size());
        //DEBUG_TIMING: sw.println("EnsureStructureCapcacity").restart();
        for (int i = 0; i < scenes.size(); i++) {
            RnaScene scene = scenes.get(i);
            if (!seq.equals(getRNASequence(scene)))
                throw new RnaBackendException("RNA sequences differ between structures."); //  (Multiple RNA structures can only be saved together in a single file if all of them share the same sequence.)
            //DEBUG_TIMING: sw.println("verifySeq" + i).restart();
            int structureIndex = i + 1;
            if (scene.title != null)
                setStructureTitle(rna, structureIndex, scene.title);
            //DEBUG_TIMING: sw.println("setStructureTitle" + i).restart();
            updateBackendBonds(scene, rna, structureIndex);
            //DEBUG_TIMING: sw.println("updateBackendBonds" + i).restart();
        }
        return rna;
    }

    private static int backendPairIndex(Nuc n) {
        return n.indexInScene() + n.getStrand().getIndex() * InterMolecularBreakLength + 1;
    }

    public static String getRNASequence(RnaScene scene) {
        StringBuilder sb = new StringBuilder(scene.getNucCount() + (scene.strands.size() - 1) * getInterStrandLinker().length());
        for (int i = 0; i < scene.strands.size(); i++) {
            Strand s = scene.strands.get(i);
            if (i != 0)
                sb.append(getInterStrandLinker()); // indicates a break between nucleic acid molecules
            for (Nuc n : s)
                sb.append(Nucleobase.getVerifiedSymbol(n.symbol));
        }
        return sb.toString();
    }

//    private static Strand getSingleStrand(RnaScene scene) {
//        if (scene.strands.size() != 1)
//            throw new IllegalStateException("Currently only a single strand per scene is supported.");
//        return scene.strands.get(0);
//    }

//    private void updateBackendBonds(int structureIndex) {
//        updateBackendBonds(getSingleStrand(getScene(structureIndex)), this.rna, structureIndex);
//    }

    /**
     * Clear the RNA structure's bonds and re-read them from the RnaScene.
     */
    protected static void updateBackendBonds(RnaScene scene, RNA rna, int structureIndex) {
        List<Bond> bonds = scene.getBonds();
        rna.RemovePairs(structureIndex, false);
        rna.RemoveConstraints();
        for (Bond b : bonds) {
            //if (b.getNuc5().getStrand() == strand && b.getNuc5().getStrand() == strand)
            int p1 = backendPairIndex(b.getNuc5());
            int p2 = backendPairIndex(b.getNuc3());
            switch (b.type) {
                case Pseudo:
                    // DO NOTHING.
                    // rna.SpecifyPair(p1, p2, structureIndex);
                    break; // do not add as a pair.
                case Prohibited:
                    rna.ForceProhibitPair(p1, p2);
                    break;
                case Forced:
                    rna.ForcePair(p1, p2);
                    // fall-through to default.
                case Special: // fall-through to default.
                case Standard:// fall-through to default.
                default:
                    rna.SpecifyPair(p1, p2, structureIndex);
                    break;
            }
        }
    }

    public static void checkForErrors(RNA rna) throws RnaBackendException {
        if (rna.GetErrorCode() != 0)
            throw new RnaBackendException(rna.GetFullErrorMessage(), rna.GetErrorCode());
    }

    /** Create new RnaScene from an RNA structure. */
    public static RnaScene createSceneFromRNA(RNA rna, int structureIndex, boolean includeBonds)
            throws RnaBackendException {
        checkForErrors(rna);
        RnaScene scene = new RnaScene();
        scene.title = getStructureTitle(rna, structureIndex); //) rna.GetCommentString
        //if (!isEmpty(scene.title)) scene.title = scene.title.trim(); // e.g. remove trailing LF
        Strand strand = scene.strands.first();
        strand.addAll(rna.GetSequence());
        splitStrandsAtLinkers(scene);

//        for (char c : rna.GetSequence().toCharArray()) {
//            if (isInterStrandLinker(c)) {
//                if (strand.size() != 0)
//                    strand = scene.addStrand();
//                continue;
//            }
//            strand.add(Character.toString(c));
//        }

//        if (scene.strands.size() == 1)
//            scene.strands.get(0).title = scene.title;
//        else for (int i = 0; i < scene.strands.size(); i++) {
//            scene.strands.get(i).title = scene.title + "-- Strand " + i;
//        }

        if (includeBonds)
            updateSceneBonds(scene, rna, structureIndex);
        return scene;
    }

    public static void splitStrandsAtLinkers(final RnaScene scene) {
        boolean inLinker = false;
        List<Nuc> linkers = null, strandStarts = null;
        for (Strand strand : scene.strands) {
            for (Nuc n : strand)
                if (!isEmpty(n.symbol) && n.symbol.charAt(0) == RnaFileIO.InterMolecularBreakChar) {
                    inLinker = true;
                    if (linkers == null) {
                        linkers = new ArrayList<>(6);
                        strandStarts = new ArrayList<>(3);
                    }
                } else if (inLinker) {
                    strandStarts.add(n);
                    inLinker = false;
                }
            inLinker = false; // just in case the linker was at the end of the sequence.
            if (linkers != null && linkers.size() != 0)
                strand.removeAll(linkers);
        }
        if (strandStarts != null && strandStarts.size() != 0)
            for (Nuc n : strandStarts)
                scene.strands.divideStrand(n.getStrand(), n.indexInStrand());
    }

//    /**
//     * This function attempts to deal with the RNAstructure convention of inserting
//     * three 'I' characters into a sequence to indicate the break between two
//     * subsequent RNA strands (between which there may be inter-molecular pairs)
//     * e.g.
//     *   RNA sequence: AGUGAGGUCUIIIGUATTUGAUU  (two 10-base sequences joined by III)
//     *   RnaScene:  S1=AGUGAGGUCU S2=GUATTUGAUU
//     *   (Note that RNA indexing is 1-based while RnaScene indexing is 0-based)
//     *   So the RNA index  14 corresponds to the first base in S2, which has index 0 in S2 and 10 overall.
//     */
//    private Nuc getSceneNucFromBackendIndex(RnaScene scene, int backendPos, int[] index) {
//        int pos = backendPos - 1; // account for backend index being 1-based
//        int strands = scene.strands.size();
//        if (strands == 1)
//            return scene.strands.get(0).get(pos);
//
//        if ((strands-1)*2 != index.length)
//            throw new IllegalStateException("The backend strand index does not match the current scene.");
//
//        // int[end1,start2]  -- 2 strands
//        // int[end1,start2,end2,start3]  -- 3 strands
//        int strand = 0;
//        while(strand < strands) {
//            if (index[strand*2] > pos) { // if the end_# > pos
//                if (strand != 0) pos -= index[(strand - 1) * 2 + 1]; // subtract the start# preceeding this position
//                if (pos < 0)
//                    throw new IndexOutOfBoundsException("Invalid nucleobase index: " + backendPos + " (corresponds to intermolecular linker).");
//                return scene.strands.get(strand).get(pos);
//            }
//            strand++;
//        }
//        throw new IndexOutOfBoundsException("Invalid nucleobase index: " + backendPos + " (beyond end of sequence).");
//    }
//    public void updateSceneBonds(int structureIndex) {
//        updateSceneBonds(getSingleStrand(getScene(structureIndex)), rna, structureIndex);
//    }

    /** Clear the RnaScene's bonds and re-write them from the RNA structure. */
    protected static void updateSceneBonds(RnaScene scene, RNA rna, int structureIndex) {
        //List<Bond> newBonds = new ArrayList<>();
        int numBases = rna.GetSequenceLength();
        //verifySequenceLength(strand, numBases);
        //int strandIndex = strand.getIndex();
        // Try to preserve information about existing bonds instead of clearing it.
        //breakBonds(strand);
        NucIndexConverter nucIndex = NucIndexConverter.create(rna.GetSequence());
        List<Bond> oldBonds = scene.getBonds();
        scene.clearBonds();
        for (int i = 1; i <= numBases; i++) {
            int pair = rna.GetPair(i, structureIndex);
            if (pair > i) {
                scene.addBond(nucIndex.getSceneNuc(scene, i), nucIndex.getSceneNuc(scene, pair));
            }
        }

        // Attempt to preserve bond types.
        for (Bond bond : oldBonds) {
            if (bond.n3.isPaired())
                bond.n3.getPairBond().type = bond.type;
            if (bond.n5.isPaired())
                bond.n5.getPairBond().type = bond.type;
        }
    }

    public static boolean isInterStrandLinker(char c) { return c == InterMolecularBreakChar; }
    public static String getInterStrandLinker() { return InterMolecularBreakStr; }

    /** Calculates the free energy of the structure in kcal/mol */
    public static double calculateEnergy(final RnaScene scene) throws RnaBackendException {
        RNA rna = createRNAfromScene(scene, true, false, false);
        breakBackendPseudoknots(scene, rna, 1);
        double energy = rna.CalculateFreeEnergy(1, true);
        if (rna.GetErrorCode() != 0)
            throw new RnaBackendException(rna.GetFullErrorMessage(), rna.GetErrorCode());
        return energy;
    }

    /**
     * The backend uses "III" (or some other character in extended nucleobase alphabets) to
     * indicate an inter-molecular linker. I.e. a single sequence actually encodes two (or potentially more)
     * RNA molecules. In this Java interface, each molecule is stored in a separate "Strand" in an RnaScene.
     * This class builds an index that helps convert between the two representations by storing the position
     * at which an "I" is encountered (which ends the first strand) as well as the position when the next valid
     * base is encountered (which starts the second strand).
     * Usually, there are three consecutive "I" characters, but this function assumes that that number is subject
     * to change, which is why both the start and end positions of each linker are stored.
     */
    abstract static class NucIndexConverter {
        public abstract Nuc getSceneNuc(RnaScene scene, int backendIndex);
        public abstract int getBackendIndex(final Nuc sceneNuc);
        public static NucIndexConverter create(final String sequence) {
            char[] seq = sequence.toCharArray();
            int count = countLinkers(seq);
            if (count == 0) return SingleStrand;
            int[] linkers = new int[count * 2];
            readLinkers(seq, linkers);
            if (count == 1) return new DoubleStrand(linkers[0], linkers[1]);
            return new MultiStrand(linkers);
        }
        public static final NucIndexConverter SingleStrand = new NucIndexConverter() {
            @Override
            public Nuc getSceneNuc(final RnaScene scene, final int backendIndex) {
                return scene.strands.get(0).get(backendIndex - 1);
            }
            @Override
            public int getBackendIndex(final Nuc sceneNuc) {
                return sceneNuc.indexInScene() + 1;
            }
        };

        public static class DoubleStrand extends NucIndexConverter {
            final int linkerStart, linkerEnd;
            public DoubleStrand(final int linkerStart, final int linkerEnd) {
                this.linkerStart = linkerStart;
                this.linkerEnd = linkerEnd;
            }
            @Override
            public Nuc getSceneNuc(final RnaScene scene, final int backendIndex) {
                int pos = backendIndex - 1;
                if (pos < linkerStart)
                    return scene.strands.get(0).get(pos);
                return scene.strands.get(1).get(pos - linkerEnd);
            }
            @Override
            public int getBackendIndex(final Nuc sceneNuc) {
                if (sceneNuc.getStrand().getIndex() == 0)
                    return sceneNuc.indexInStrand() + 1;
                return sceneNuc.indexInStrand() + 1 + linkerEnd;
            }
        }

        public static class MultiStrand extends NucIndexConverter {
            final int[] linkers;
            public MultiStrand(final int[] linkers) { this.linkers = linkers; }
            @Override
            public Nuc getSceneNuc(RnaScene scene, int backendIndex) {
                int pos = backendIndex - 1;
                if (pos < linkers[0])
                    return scene.strands.get(0).get(pos); // this takes care of all nucs in the first strand.

                int strand = 0, lastStrand = linkers.length / 2;
                while (++strand < lastStrand) // starting index=1, which means it is actually the *second* strand.
                    if (pos < linkers[strand * 2]) // this handles all nucs in strands 2 through N-1
                        return scene.strands.get(strand).get(pos - linkers[strand * 2 - 1]); // subtract the end-position of the previous linker.
                // Otherwise it is in the last strand.
                return scene.strands.get(lastStrand).get(pos - linkers[lastStrand * 2 - 1]);
            }
            @Override
            public int getBackendIndex(final Nuc sceneNuc) {
                if (sceneNuc.getStrand().getIndex() == 0)
                    return sceneNuc.indexInStrand() + 1;
                return sceneNuc.indexInStrand() + 1 + linkers[sceneNuc.getStrand().getIndex() * 2 - 1];
            }
        }

        private static int countLinkers(final char[] sequence) {
            boolean prev = false;
            int i = 0;
            for (char c : sequence) {
                boolean current = isInterStrandLinker(c);
                if (current && !prev) i++;
                prev = current;
            }
            return i;
        }
        private static void readLinkers(final char[] sequence, final int[] linkers) {
            // int[end1, start2, end2, start3] ...
            int i = 0, pos = 0;
            boolean prev = false;
            for (char c : sequence) {
                boolean current = isInterStrandLinker(c);
                if (prev != current) linkers[i++] = pos;
                prev = current;
                pos++;
            }
        }
    }

//    /**
//     * Remove all bonds in an RnaScene that involve the specified Strand.
//     */
//    private static void breakBonds(final Strand strand) {
//        java.util.List<Bond> bonds = strand.getScene().bonds;
//        int strandIndex = strand.getIndex();
//        if (bonds.size() == 0)
//            return;
//        Bond[] arr = bonds.toArray(new Bond[bonds.size()]);
//        for (int i = arr.length - 1; i > -1; i--) {
//            Bond b = arr[i];
//            if (b.S1 == strandIndex || b.S2 == strandIndex)
//                bonds.remove(i);
//        }
//    }

    public static void breakBackendPseudoknots(final RnaScene scene, final RNA rna, final int structureNumber) {
        Set<Bond> pseudo = Motif.findPseudoKnots(scene);
        for (Bond b : pseudo)
            rna.RemoveBasePair(backendPairIndex(b.n5), structureNumber);
    }

    public static boolean getBestOrientation(double[] oldVals, double newVals[], double[] params, double allowedError) {
        // for all sets of two distinct points that are equi-distant in both models,
        // determine the rotation and translation that would make them congruent, then calculate how many other points match.
        allowedError = allowedError * allowedError; // square it so we can use the squared distances.
        Vec2D v1 = new Vec2D();
        Vec2D v2 = new Vec2D();
        for (int i = 0; i < oldVals.length; i += 2)
            for (int j = 0; j < oldVals.length; j += 2) {
                v1.setTo(oldVals[j] - oldVals[i], oldVals[j + 1] - oldVals[i + 1]);
                v2.setTo(newVals[j] - newVals[i], newVals[j + 1] - newVals[i + 1]);
                if (v1.lengthSqr() == v2.lengthSqr()) {
                    double angle = v1.angleTo(v2);
                }
            }
        return true;
    }

//    private static void verifySequenceLength(final Strand strand, int numBases) {
//        if (strand.size() != numBases)
//            throw new IllegalStateException("A strand has the wrong number of nucleotides.");
//    }
//    private void verifySequenceLength(final Strand strand) {
//        verifySequenceLength(strand, this.numBases);
//    }

    public interface AsyncTask<TResult> {
        Object getTag();
        int getProgress();
        String getStatus();
        boolean isDone();
        default boolean isCanceled() {
            return false;
        }
        default boolean canCancel() {
            return false;
        }
        default void cancel() {
            throw new UnsupportedOperationException("This Task cannot be canceled.");
        }
        TResult getResult();
        Exception getError();
        default boolean hadError() {
            return getError() != null;
        }
    }

    /**
     * A task that is run in the RNAstructure Backend (JNI).
     * These jobs update a (c++) ProgressHandler, but do not source
     * update events.
     * @param <TResult> The result that is produced by the calculation.
     */
    public static abstract class BackendCalc<TResult> extends BackgroundWork<TResult> {
        protected ProgressHandler backendProgress = new SimpleProgressHandler();
        protected int _stepWork=100, _workComplete;

        @Override
        public int getProgress() { return (int)(_workComplete+backendProgress.progress()/100f*_stepWork); }
        @Override
        public boolean isCanceled() {
            return backendProgress.canceled();
        }
        @Override
        public void cancel() {
            backendProgress.cancel();
        }
        protected void nextStep(int workInStep) {
            if (_stepWork!=100)
                _workComplete+=_stepWork;
            _stepWork=workInStep;
            backendProgress.update(0);
            notifyUpdate();
        }
        protected void nextStep(int workInStep, String status, boolean canCancel) {
            _canCancel=canCancel;
            _status = status;
            nextStep(workInStep);
        }
    }

    public abstract static class BackgroundWork<TResult>
     implements AsyncTask<TResult> {
        protected boolean _isDone, _isCanceled, _canCancel;
        protected int _progress = 0;
        protected Object _tag; String _status;
        protected TResult _result; Exception _err;
        protected Consumer<? super BackgroundWork<TResult>> _update, _complete;

        public Object getTag() { return _tag; }
        public int getProgress() { return _progress; }
        public String getStatus() { return _status; }
        public boolean isDone() { return _isDone; }
        public boolean isCanceled() { return _isCanceled; }
        public boolean canCancel() { return _canCancel; }
        public void cancel() {
            if (!_canCancel)
                throw new UnsupportedOperationException("This Task cannot be canceled.");
            _isCanceled = true;
            notifyUpdate();
        }
        public TResult getResult() { return _result; }
        public Exception getError() { return _err; }

        protected void setError(String message) { setError(new Exception(message)); }
        protected void setError(Exception ex) {
            _err = ex;
            notifyUpdate();
        }

        protected void update(int progress, String status) { _status = status; update(progress); }
        protected void update(int progress) {
            _progress = progress;
            notifyUpdate();
        }

        protected abstract TResult calcResult() throws Exception;

        public void start() {
            SwingWorker w = new SwingWorker() {
                @Override
                protected Object doInBackground() throws Exception {
                    try {
                        _result = calcResult();
                    } catch (Exception ex) {
                        _result = null;
                        setError(ex);
                    }
                    _isDone = true;
                    notifyUpdate();
                    return null;
                }
                @Override
                protected void done() {
                    _complete.accept(BackgroundWork.this);
                }
            };
            w.execute();
        }
        protected void notifyUpdate() {
            if (_update==null) return;
            if (SwingUtilities.isEventDispatchThread())
                _update.accept(this);
            else
                SwingUtilities.invokeLater( ()->_update.accept(this));
        }

        public void whenDone(Consumer<? super BackgroundWork<TResult>> actionOnComplete) {
            _complete = actionOnComplete;
        }
        public void onUpdate(Consumer<? super BackgroundWork<TResult>> actionOnUpdate) {
            _update = actionOnUpdate;
        }
    }

    public static BackendCalc<RnaSceneGroup> foldSeq(String sequence,int maxStructures) {
        //final int maxStructures = 10;
//        Strand s = getSingleStrand(getScene(structureIndex));
//        RnaFileIO copy = fromScene(getScene(structureIndex), 'U');
        return new BackendCalc<RnaSceneGroup>() {
            @Override
            protected RnaSceneGroup calcResult() throws Exception {
                nextStep(5, "Loading Sequence...", true);
                if (isCanceled()) return null;
                RNA rna = new RNA(sequence, true); //createRNAfromScene(scene, false, false, false);
                if (rna.GetErrorCode()!=0)
                    throw new RnaBackendException("Failed to Fold RNA Sequence:  " + rna.GetFullErrorMessage());
                rna.SetProgress(backendProgress);
                nextStep(95, "Folding Sequence...", true);
                if (isCanceled()) return null;
                int errCode = rna.FoldSingleStrand(20, maxStructures, defaultWindowSize(rna),
                        null, 30, false, true, false);
                if (errCode != 0)
                    throw new RnaBackendException("Failed to Fold RNA Sequence:  " + RNA.GetErrorMessage(errCode), errCode);
                return getRNAstructureGroup(rna);
            }
        };
    }

//    /** destructive operation that re-folds the RNA sequence and creates all new structures. */
//    protected static void refold_Internal(RNA rna, int maxStructures, PropertyChangeListener progressListener) throws
//            RnaBackendException {
//        //verifyDataPath();  now the datapath needs to be set BEFORE the RNA object is created.
////        TProgressDialog progress = new TProgressDialog();
////        rna.SetProgress(progress);
////        Timer t = new Timer(250, new ActionListener() {
////            int prev = -1;
////            @Override
////            public void actionPerformed(final ActionEvent e) {
////                int current = progress.progress();
////                if (current != prev)
////                    progressListener.propertyChange(new PropertyChangeEvent(this, "Progress",prev, current));
////            }
////        });
//
//        // int FoldSingleStrand(const float percent=20, const int maximumstructures=20, const int window=5, const char savefile[]="", const int maxinternalloopsize = 30, bool mfeonly=false, bool simple_iloops = true);
//        int errCode = rna.FoldSingleStrand(20, maxStructures, defaultWindowSize(rna), null, 30, false, true, false);
//        if (errCode != 0)
//            throw new RnaBackendException("Failed to Refold RNA Sequence:  " + RNA.GetErrorMessage(errCode), errCode);
//    }

//    private void setDefaultPositions(Strand strand) {
//        int length = strand.size();
//
//        final int groupCount = 10;
//        final int lineCount = 40;
//        final int maxSpaces = length < 20 ? lineCount : lineCount + (int) Math.ceil(lineCount / (float) groupCount) - 1;
//        int lines = 0;
//        int spaces = 0;
//
//        float x, y;
//        for (int i = 0; i < length; i++) {
//            Nuc n = strand.get(i);
//            y = nucRadius + lines * nucRadius * 2.6f;
//            int pos = (lines % 2) == 0 ? spaces : maxSpaces - spaces - 1;
//            x = nucRadius + pos * nucRadius * 2.4f;
//            n.location.setLocation(x, y);
//            if ((i + 1) % lineCount == 0) {
//                lines++;
//                spaces = 0;
//            } else {
//                spaces++;
//                if (length >= 20 && 0 == (i + 1) % groupCount)
//                    spaces++; //add an extra gap for groups
//            }
//        }
//    }

    private static int defaultWindowSize(RNA rna) {
        // Get the sequence length.
        int length = rna.GetSequenceLength();

        // Return a window size based on the length.
        if (length > 1200) { return 20; } else if (length > 800) { return 15; } else if (length > 500) {
            return 11;
        } else if (length > 300) { return 7; } else if (length > 120) { return 5; } else if (length > 50) {
            return 3;
        } else { return 2; }
    }

    private static void verifyDataPath() throws RnaBackendException {
        String path = getDataPath();
        if (path == null)
            throw new RnaBackendException("The thermodynamic data tables could not be found.\nPlease set the DATAPATH environment variable to the directory where they are stored.\n\t(e.g. <PATH-TO-RNAstructure>/data_tables).");
        RNABackend.setDataPath(path);
    }
    private static String getDataPath() {
        if (DATA_PATH != null) return DATA_PATH;
        if (verifyDataFiles(RNABackend.getDataPath()))
            return DATA_PATH;
        // DATAPATH was either NULL or not found. In either case, search for one in standard directories.
        String paths[] = {".", "..", "./resources", "../Resources"};
        String fullPaths[] = {"{path}/", "{path}/data_tables", "{app}/{path}/", "{app}/{path}/data_tables"};
        String appPath = PathTools.getAppPath(RnaFileIO.class);
        if (appPath == null) appPath = ".";
        for (String path : paths) {
            for (String full : fullPaths) {
                String fullPath = full.replace("{app}", appPath).replace("{path}", path);
                if (verifyDataFiles(fullPath))
                    return DATA_PATH;
            }
        }
        return null;
    }
    private static boolean verifyDataFiles(String path) {
        final String filter = "*.dat";
        if (path == null) return false;
        File dir = new File(path);
        if (dir.exists())
            try {
                String realpath = dir.getCanonicalPath();
                if (PathTools.listFiles(realpath, filter).size() != 0) {
                    DATA_PATH = realpath;
                    log.debug("The thermodynamic data tables were located at \"" + realpath + "\".");
                    return true;
                }
            } catch (Exception ex) {
                //Do nothing. continue the loop.
            }
        return false;
    }

    public static int getRNAType(FileType type) {
        switch (type) {
            case Structure:
            case CT:
                return RNABackend.FILE_CT;

            case DotBracket:
                return RNABackend.FILE_DBN;
            case FoldingSav:
                return RNABackend.FILE_SAV;
            case PartitionSav:
                return RNABackend.FILE_PFS;

            case Sequence:
            case Seq:
            case Fasta:
            case SeqText:
                return RNABackend.FILE_SEQ;
            default:
                return -1;
        }
    }

    public static RnaSceneGroup readFileAsScene(String path, FileType type) throws RnaBackendException {
        verifyDataPath();
        int backendType = getRNAType(type);
        RNA rna = new RNA(path, backendType, ALPHABET_RNA, true, true);
        rna.EnsureStructureCapcacity(1); // if it's a SEQ file, RNA doesn't create a structure by default.
        if (rna.GetErrorCode() != 0)
            throw new RnaBackendException(rna.GetFullErrorMessage(), rna.GetErrorCode(), path);
        // fmt("Failed to load %s File \"%s\": %s\nFull Path: \"%s\"\nInternal Error Number: %s", type.name(), tmp.getName(), message, getCanonicalPath(tmp), err)
        System.out.println("Done readRNAFile @" + System.currentTimeMillis());
        return getRNAstructureGroup(rna);
    }

    public static RnaSceneGroup getRNAstructureGroup(RNA rna) throws RnaBackendException {
        RnaSceneGroup group = new RnaSceneGroup();
        group.setTitle(rna.GetStructure().GetSequenceLabel());
        int count = rna.GetStructureNumber();
        for (int i = 1; i <= count; i++) {
            // copy the newly folded structure(s) from the RNA copy to the current strand.
            group.add(createSceneFromRNA(rna, i, true));
        }
        // updateSceneBonds(s, copy.rna, 1);
        // updateBackendBonds(structureIndex);
        return group;
    }

    public static boolean writeSequenceFile(RnaScene scene, final String path, final FileType type, boolean append)
            throws RnaBackendException, IOException {
        int fmt;
        switch (type) {
            case Fasta:
                fmt = SEQUENCE_FMT_FASTA;
                break;
            case Seq:
                fmt = SEQUENCE_FMT_SEQ;
                break;
            case SeqText:
                fmt = SEQUENCE_FMT_TEXT;
                break;
            default:
                throw new IllegalArgumentException("Invalid Sequence file format: " + type);
        }
        RNA rna = createRNAfromScene(scene, false, true, true);
        structure s = rna.GetStructure();
        if (s.writeseq(path, fmt, append) == 0) {
            String error = s.GetErrorDetails();
            throw new IOException(isEmpty(error) ? "Error writing file." : error);
        }
        return true;
    }
    public static boolean writeStructureFile(RnaSceneGroup scenes, final String path, final FileType type, boolean append)
            throws IOException, RnaBackendException {

        RNA rna = writeRNAstructureGroup(scenes, true, true);
        switch (type) {
            case Structure:
            case DotBracket:
                rna.WriteDotBracket(path); // append
                break;
            case CT:
                rna.WriteCt(path); // append
                break;
            default:
                throw new IllegalArgumentException("Unknown file type: " + type);
        }
        int err = rna.GetErrorCode();
        if (err != 0)
            throw new RnaBackendException(rna.GetFullErrorMessage(), err, path);
        return true;
    }

//    private void writeSeq(RNA rna, String path, FileType type, boolean append) throws IOException, RnaBackendException {
//        int fmt;
//        switch (type) {
//            case Fasta: fmt = SEQUENCE_FMT_FASTA; break;
//            case Seq: fmt = SEQUENCE_FMT_SEQ;  break;
//            case SeqText: fmt = SEQUENCE_FMT_TEXT;  break;
//            default: throw new IllegalArgumentException("Invalid Sequence file format: " + type);
//        }
//        structure s = rna.GetStructure();
//        if (s.writeseq(path, fmt, append) == 0) {
//            String error = s.GetErrorDetails();
//            throw new IOException(isEmpty(error)?"Error writing file.":error);
//        }
//            //throw new RnaBackendException(fmt("Failed to save %s File \"%s\": %s\nFull Path: \"%s\"", type.name(), outFile.getName(), s.GetErrorDetails(), getCanonicalPath(outFile)), 0, outFile.toString());
//    }
//

//    public static FileType guessFileType(final String path) {
//        String ext = PathTools.getExt(path).toUpperCase();
//        FileType type;
//        if (ext.equals("." + Program.DrawingFileExtension.toUpperCase()))
//            type = FileType.Drawing;
//        else
//            switch (ext) {
//                case ".CT":
//                case ".DBN":
//                case ".BRACKET":
//                    type = FileType.Structure;
//                    break;
//                case ".SEQUENCE":
//                case ".SEQ":
//                case ".FASTA":
//                    type = FileType.Sequence;
//                    break;
//                default:
//                    type = null; //could not be determined.
//            }
//        return type;
//    }

    public static String getStructureTitle(RNA rna, final int structureIndex) {
        String s = rna.GetCommentString(structureIndex);
        if (s.endsWith("\n"))
            return s.substring(0, s.length() - 1);
        return s;
    }
    public static void setStructureTitle(RNA rna, final int structureIndex, String title) {
        if (title == null)
            title = "";
        String s = rna.GetCommentString(structureIndex);
        if (s.endsWith("\n") && !title.endsWith("\n"))
            title += "\n";
        rna.GetStructure().SetCtLabel(title, structureIndex);
    }

    public static void identifyPseudoKnotsBackend(final RnaScene scene) throws RnaBackendException {
        StopWatch sw = new StopWatch(true);
        RNA rna = createRNAfromScene(scene, true, true, false);
        sw.println("createRNAfromScene").restart();
        rna.BreakPseudoknot(false, 1);
        sw.println("rna.BreakPseudoknot").restart();
        RnaScene updated = createSceneFromRNA(rna, 1, true);
        sw.println("createSceneFromRNA").restart();
        List<Nuc> updatedNucs = updated.allNucs();
        for (Nuc n : scene.allNucs()) {
            if (n.isPaired()) {
                boolean paired = updatedNucs.get(n.indexInScene()).isPaired(); // is this nuc paired in the backend version (with no pseudoknots)
                if (n.getPairBond().type == BondType.Default && !paired)
                    n.getPairBond().type = BondType.Pseudo; // set to pseudo if it is a normal bond in the current scene, but missing from the backend.
                else if (n.getPairBond().type == BondType.Pseudo && paired)
                    n.getPairBond().type = BondType.Default; // set to Normal if it is a pseudo  bond in the current scene, but found (as a standard bond) in the backend.
            }
        }
        sw.println("UpdatedNucs").restart();
    }
    static Pattern stripEnergyTitleRegex = Pattern.compile("^\\s*ENERGY\\s*=\\s*-?[\\d.]+\\s*", Pattern.CASE_INSENSITIVE);
    public static String stripEnergyTitle(String title) {
        if (title == null) return null;
        return stripEnergyTitleRegex.matcher(title).replaceFirst("");
    }
}
