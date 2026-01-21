package ur_rna.StructureEditor.services.fileIO;

import ur_rna.StructureEditor.FileType;
import ur_rna.StructureEditor.models.*;
import ur_rna.StructureEditor.services.SceneDrawMode;
import ur_rna.Utilities.*;
import ur_rna.Utilities.SimpleDataSerializer.SSONList;
import ur_rna.Utilities.SimpleDataSerializer.SSONMap;

import java.awt.*;
import java.io.*;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.regex.MatchResult;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static ur_rna.Utilities.Strings.*;

/**
 * There are two types of native drawing file formats.
 *
 * CTE "Extended CT File"               - Similar to CT, but with additional columns.
 * NSD "Nucleotide Structure Drawing"   - Hierarchical Structured Object-based format that stores
 *      Nucleotides, Bonds, and Styles in a format that is richer and more extensible than CTE
 *
 */
public class DrawingFileIO {

    //<editor-fold desc="Previous code for stripping NSD information from CT files">
    //    StringBuilder fullText;
    //    int[] sourceLineInfo;
    //    int lineCount;
    //    Object parsed;
    //
    //    public void parseStart() {
    //        fullText = new StringBuilder();
    //        sourceLineInfo = new int[300];
    //        lineCount = 0;
    //        //fullText.append('['); // assume data represents an array of RNAScenes
    //    }
    //
    //    // Each content line can be continued onto the next line by appending a '\' character before the newline.
    //    public void parseContentLine(String line, int sourceLineNumber, int sourceColumnNumber) throws SyntaxErrorException {
    //        if (sourceLineInfo.length < (lineCount+1) * 3) {
    //            sourceLineInfo = Arrays.copyOf(sourceLineInfo, sourceLineInfo.length * 2);
    //        }
    //        sourceLineInfo[lineCount * 3] = fullText.length();
    //        sourceLineInfo[lineCount * 3 + 1] = sourceLineNumber;
    //        sourceLineInfo[lineCount * 3 + 2] = sourceColumnNumber;
    //        fullText.append(line).append('\n');
    //        lineCount++;
    //    }
    //
    //    public void parseEnd() throws SyntaxErrorException {
    //        int[] index = Arrays.copyOf(sourceLineInfo, lineCount * 3);
    //        sourceLineInfo = null;
    //        //fullText.append(']');
    //        SimpleDataSerializer p = new SimpleDataSerializer();
    //        parsed = p.parse(fullText, index);
    //    }
    //
    //    public RnaSceneGroup getParsed() throws IOException {
    //    }
    //</editor-fold>

    //<editor-fold desc="NSD Drawing Files">
    /** Constants for Nucleotide Structure Drawing (NSD) File IO */
    private interface NSD {
        String  VERSION = "Version",
                STRANDS = "Strands", BONDS = "Pairs", STYLES = "Styles",
                TITLE = "Title", DRAW_MODE = "Mode", DRAW_FLIP = "Flip", STYLE = "Style", STRAND_NUCS = "Strand", NUC_SYMBOL = "Base", NUC_NUMBER = "HNum",
                NUC_ID = "ID", BOND_NUCS = "Pair", BOND_TYPE = "Type",
                STYLE_TEXT = "Color", STYLE_LINE = "Line", STYLE_FILL = "Fill", STYLE_BOND = "Bond", STYLE_FONT = "Font", STYLE_NUMOFF="NumOffset",
                X="X",Y="Y";
        float LOC_ZERO = 0f;
        int UNSTYLED = -1;
    }
    public static RnaSceneGroup readNsdDrawingFile(String path) throws SyntaxErrorException, IOException {
        String fullText = new String(Files.readAllBytes(Paths.get(path)), StandardCharsets.UTF_8);
        Object parsed = new SimpleDataSerializer().parse(fullText, null);
        return readNsdDrawingFile(parsed);
    }
    public static RnaSceneGroup readNsdDrawingFile(Reader input) throws SyntaxErrorException, IOException {
        int charsRead = 0;
        char[] buffer = new char[1024];
        StringBuilder sb = new StringBuilder();
        while((charsRead=input.read(buffer))>0)
            sb.append(buffer, 0, charsRead);
        Object parsed = new SimpleDataSerializer().parse(sb, null);
        return readNsdDrawingFile(parsed);
    }
    private static RnaSceneGroup readNsdDrawingFile(Object parsedNSDData) throws SyntaxErrorException {
        if (parsedNSDData == null) throw new IllegalStateException("Document content must be parsed before applying it to an RNA structure.");
        String location = "start of data";
        try {
            SSONList sdsScenes = expectList(parsedNSDData, "Expected a list of RNA scenes.");
            SSONMap mInfo = expectMap(sdsScenes.isEmpty()?null:sdsScenes.get(0), "Expected a Drawing file descriptor, e.g. { Version: ... }");
            Version version = new Version(Convert.toString(mInfo.get(NSD.VERSION), ""));

            if (version.isEmpty()) throw new IllegalArgumentException("The version descriptor was not found.");
            RnaSceneGroup scenes = new RnaSceneGroup();
            scenes.setTitle(mInfo.get(NSD.TITLE,"Untitled"));
            int sceneNum = 0;
            // Each file can contain a list of Scenes, each of which can contain one or more Nucleic acid sequences.
            for (Object objScene : sdsScenes.subList(1, sdsScenes.size())) {
                ++sceneNum;
                RnaScene scene = new RnaScene();
                HashMap<String, Nuc> idMap = new HashMap<>();
                HashMap<String, String> pairMap = new HashMap<>();
                scenes.add(scene);
                String sceneLoc = location = "scene " + sceneNum; // for error messages.
                SSONMap mScene = expectMap(objScene, "Expected an RNA 'scene' descriptor, e.g. { Title: ... Strands: ... Bonds: ... }");
                scene.title = mScene.get(NSD.TITLE, "Structure" + sceneNum);
                scene.drawMode =  SceneDrawMode.valueOf(mScene.get(NSD.DRAW_MODE, SceneDrawMode.Standard.name()));
                scene.drawFlipped =  mScene.get(NSD.DRAW_FLIP, false);

//                location = sceneLoc; //for error messages
//                // Load the list of styles (aka formats)
//                SSONList listStyles = expectList(mScene.get(NSD.STYLES), "Expected a list of Styles");
//                int styleNum = 0;
//                for (Object objStyle : listStyles) {
//                    ++styleNum;
//                    location = sceneLoc +", style #"+styleNum;
//                    SSONMap mStyle = expectMap(objStyle, "Expected a Style descriptor, e.g. { Color: ...  }");
//                    NucStyle s = scene.addStyle();
//                    s.fillColor = decodeColor(mStyle.get(NSD.STYLE_FILL), null);
//                    s.textColor = decodeColor(mStyle.get(NSD.STYLE_TEXT), null);
//                    s.outlineColor = decodeColor(mStyle.get(NSD.STYLE_LINE), null);
//                    s.font = decodeFont(mStyle.get(NSD.STYLE_FONT), null);
//                }

                // Load the list of strands (aka sequences)
                SSONList listStrands = expectList(mScene.get(NSD.STRANDS), "Expected a list of Strands");
                int strandNum = 1;
                for (Object objStrand : listStrands) {
                    location = sceneLoc +", strand #"+strandNum;
                    SSONMap mStrand = expectMap(objStrand, "Expected an RNA strand descriptor, e.g. { Title: ... Bases: ... }");
                    Strand strand = scene.getStrand(strandNum-1, true); //addStrand(mStrand.get(NSD.TITLE, "Sequence" + strandNum));
                    SSONList listNucs = expectList(mStrand.get(NSD.STRAND_NUCS), "Expected a list of Nucleotides");
                    // Load the list of nucleotides in each strand
                    for (Object objNuc : listNucs) {
                        location = String.format("scene %s, strand #%s, position %s.", sceneNum, strandNum, strand.size() + 1);
                        SSONMap mNuc = expectMap(objNuc, "Expected a Nucleotide (e.g. { Base: G, X: 3, Y: 5 ...  })");
                        Nuc n = new Nuc(mNuc.get(NSD.NUC_SYMBOL, "?"), mNuc.get(NSD.X, NSD.LOC_ZERO), mNuc.get(NSD.Y, NSD.LOC_ZERO));
                        strand.add(n);
                        n.number = mNuc.get(NSD.NUC_NUMBER, -1);
                        String id;
                        Object oid = mNuc.get(NSD.NUC_ID);
                        if (oid == null)
                            id = ""+(n.indexInScene()+1);
                        else if (oid instanceof Number)
                            id = fmtNumber(((Number)oid).doubleValue());
                        else
                            id = oid.toString();

                        if (idMap.containsKey(id))
                            throw new SyntaxErrorException("The ID "+id+" was used for more than one base.");

                        idMap.put(id, n);
                        //String pair = mNuc.get(NSD.BOND_NUCS, "");
                        //scene.setStyle(n, mNuc.get(NSD.STYLE, NSD.UNSTYLED));
                        if (mNuc.containsAnyKey(NSD.STYLE_FILL, NSD.STYLE_TEXT, NSD.STYLE_LINE, NSD.STYLE_FONT,NSD.STYLE_BOND, NSD.STYLE_NUMOFF)) {
                            n.style = new NucStyle(
                                    decodeColor(mNuc.get(NSD.STYLE_FILL), null),
                                    decodeColor(mNuc.get(NSD.STYLE_TEXT), null),
                                    decodeColor(mNuc.get(NSD.STYLE_LINE), null),
                                    decodeColor(mNuc.get(NSD.STYLE_BOND), null),
                                    decodeFont(mNuc.get(NSD.STYLE_FONT), null),
                                    mNuc.get(NSD.STYLE_NUMOFF, 0.0F)
                            );
                        }
                    }
                    strandNum++;
                }
                location = sceneLoc; //for error messages
                // load the list of Bonds
                SSONList listBonds = expectList(mScene.get(NSD.BONDS), "Expected a list of Bonds" + sceneNum);
                int bondNum = 0;
                for (Object objBond : listBonds) {
                    ++bondNum;
                    location = sceneLoc +", bond #"+bondNum;
                    SSONMap mBond = expectMap(objBond, "Expected a Bond descriptor, e.g. { Pair: \"5:12\" Type: BP }");
                    BondType type = BondType.fromAbbrev(mBond.get(NSD.BOND_TYPE,""), BondType.Default);
                    String[] parts = mBond.get(NSD.BOND_NUCS, "").split(":");
                    Nuc n1 = idMap.get(parts[0].trim()), n2 = idMap.get(parts[1].trim());
                    if (n1 == null) throw new SyntaxErrorException("The ID \""+parts[0]+"\" is not associated with any nucleotide.");
                    if (n2 == null) throw new SyntaxErrorException("The ID \""+parts[1]+"\" is not associated with any nucleotide.");
                    scene.addBond(n1, n2, type);
                    //b.style = mBond.get(NSD.STYLE, NSD.UNSTYLED);
                }
            }
            return scenes;
        } catch (Exception ex) {
            throw new SyntaxErrorException("Error reading Drawing data: while parsing " + location + ": " + ex.getMessage(), ex);
        }
    }
    private static final DecimalFormat _dblFormat = new DecimalFormat("#.##########");
    private static String fmtNumber(final double value) {
        long rounded = Math.round(value);
        if (value == rounded)
            return Long.toString(rounded);
        if (Math.abs(value) > 0.0001 && Math.abs(value) < 100000)
            return _dblFormat.format(value);
        return Double.toString(value);
    }

    private static SSONMap expectMap(Object o, String errorMessage) {
        if (o instanceof SSONMap)
            return (SSONMap)o;
        throw new IllegalArgumentException(errorMessage);
    }
    private static SSONList expectList(Object o, String errorMessage) {
        if (o instanceof SSONList)
            return (SSONList)o;
        throw new IllegalArgumentException(errorMessage);
    }

    public static void writeNsdDrawingFile(RnaSceneGroup scenes, String path, boolean append) throws IOException, FormatterException {
        writeNsdDrawingFile(scenes, path, append, 2, false);
    }
    public static void writeNsdDrawingFile(RnaSceneGroup scenes, String path, boolean append, int indent, boolean strict) throws IOException, FormatterException {
        try(Writer writer = newUTF8Writer(path, append)) {
            writeNsdDrawingFile(scenes, writer, indent, strict);
        }
    }
    public static void writeNsdDrawingFile(RnaSceneGroup scenes, Writer out, int indent, boolean strict) throws IOException, FormatterException {
        SSONMap mInfo = new SSONMap();
        SSONList listScenes = new SSONList();
        // Write File-level information.
        mInfo.put(NSD.VERSION, "1.0");
        mInfo.put(NSD.TITLE, scenes.getTitle());
        listScenes.add(mInfo);

        // Write each scene
        for(RnaScene scene : scenes) {
            SSONMap mScene = listScenes.addMap();
            mScene.put(NSD.TITLE, scene.title);
            mScene.put(NSD.DRAW_MODE, scene.drawMode.name());
            mScene.put(NSD.DRAW_FLIP, scene.drawFlipped);
            SSONList listStrands = mScene.putList(NSD.STRANDS) ,
                     listBonds = mScene.putList(NSD.BONDS)/*,
                     listStyles = mScene.putList(NSD.STYLES) */;
            // Write each strand
            for(Strand s : scene.strands) {
                SSONMap mStrand = listStrands.addMap();
                //mStrand.put(NSD.TITLE, s.title);
                SSONList listNucs = mStrand.putList(NSD.STRAND_NUCS);
                // Write each nucleotide
                for(Nuc n : s) {
                    SSONMap mNuc = listNucs.addMap();
                    mNuc.put(NSD.NUC_ID, n.indexInScene()+1);
                    mNuc.put(NSD.NUC_SYMBOL, n.symbol);
                    if (n.number!=-1) mNuc.put(NSD.NUC_NUMBER, n.number);
//                    if (n.isPaired()) {
//                        Bond b = n.getPairBond();
//                        mNuc.put(NSD.BOND_NUCS, b.getOther(n).indexInScene()+1);
//                        if (b.type != BondType.Default)
//                            mNuc.put(NSD.BOND_TYPE, b.type.getAbbrev());
//                    }
                    mNuc.put("X", roundLoc(n.location.x));
                    mNuc.put("Y", roundLoc(n.location.y));
                    if (n.style != null) {
                        //mNuc.put(NSD.STYLE, n.getStyle().getIndex());
                        NucStyle ns = n.style; SSONMap m = mNuc;
                        if (ns.textColor != null) m.put(NSD.STYLE_TEXT, Colors.getName(ns.textColor));
                        if (ns.fillColor != null) m.put(NSD.STYLE_FILL, Colors.getName(ns.fillColor));
                        if (ns.outlineColor != null) m.put(NSD.STYLE_LINE, Colors.getName(ns.outlineColor));
                        if (ns.bondColor != null) m.put(NSD.STYLE_BOND, Colors.getName(ns.bondColor));
                        if (ns.font != null) m.put(NSD.STYLE_FONT, ns.font.toString());
                        if (ns.numberOffset!=0) m.put(NSD.STYLE_NUMOFF, ns.numberOffset);
                    }
                }
            }

            // Write each bond
            for(Bond b : scene.getBonds()) {
                SSONMap mBond = listBonds.addMap();
                mBond.put(NSD.BOND_NUCS, (b.getNuc5().indexInScene()+1) + ":" + (b.getNuc3().indexInScene()+1));
                mBond.put(NSD.BOND_TYPE, b.type.getAbbrev());
            }
//
//            // Write each style
//            for(NucStyle s : scene.styles) {
//                SSONMap mStyle = listStyles.addMap();
//                if (s.textColor != null) mStyle.put(NSD.STYLE_TEXT, Colors.getName(s.textColor));
//                if (s.outlineColor != null) mStyle.put(NSD.STYLE_LINE, Colors.getName(s.outlineColor));
//                if (s.fillColor != null) mStyle.put(NSD.STYLE_FILL, Colors.getName(s.fillColor));
//                if (s.font != null) mStyle.put(NSD.STYLE_FONT, s.font.toString());
//            }
        }

        SimpleDataSerializer sds = new SimpleDataSerializer();
        SimpleDataSerializer.OutputFormatSettings settings = new SimpleDataSerializer.OutputFormatSettings(indent, strict);
        sds.write(listScenes, out, settings);
    }
    //</editor-fold>

    //<editor-fold desc="Common Drawing File Tools">

    private static Font decodeFont(final Object oFont, final Font defaultValue) {
        if (oFont == null || oFont.equals("")) return defaultValue;
        if (oFont instanceof Font) return (Font)oFont;
        if (oFont instanceof String) return Font.decode((String)oFont);
        return defaultValue;
    }
    private static Color decodeColor(Object oColor, Color defaultValue) {
        if (oColor == null || oColor.equals("")) return defaultValue;
        if (oColor instanceof Color) return (Color)oColor;
        if (oColor instanceof String) return Colors.getColor((String)oColor, defaultValue);
        if (oColor instanceof Number) return new Color(((Number)oColor).intValue());
        return defaultValue;
    }
    /** round locations to the nearest 0.1 */
    private static float roundLoc(float value) { return Math.round(10*value)/10f;}

    /**
     * Simplified version of {@link Files#newBufferedWriter(Path, Charset, OpenOption...)} that
     * uses {@link StandardCharsets#UTF_8} as the charset and accepts a boolean "append" instead of
     * an {@link OpenOption} argument list.
     */
    private static Writer newUTF8Writer(String path, boolean append) throws IOException {
        OpenOption[] options = append ? new OpenOption[] { StandardOpenOption.APPEND } : new OpenOption[0];
        return Files.newBufferedWriter(Paths.get(path), StandardCharsets.UTF_8, options);
    }

    //</editor-fold>

    //<editor-fold desc="Extended CT File (CTE) IO">
    /** Constants for CTE File IO */
    private interface CTE {
        String X = "X", Y = "Y",
            TEXT_COLOR = "TC", FILL_COLOR = "FC", LINE_COLOR = "LC", BOND_COLOR = "BC",
            FONT = "FS", BOND_TYPE = "B", LBL_OFFSET="LO";

        int PAD_TOTAL = 5, PAD_INDEX = 4, PAD_SYMBOL = -4, PAD_COLOR = 7, PAD_LOC = 6, PAD_FONT = 10;

        String PROP_SEPARATOR = ";!"; // separates standard CT columns from CTE columns
    }

    /**
     * Simplifies writing tabular text data in which the table formatting is achieved by padding values with spaces.
     * This class contains several "write" methods that accept an integer "padding" value. These pad spaces onto
     * the content. The padding will be on the left or right depending on whether the sign of the padding value is
     * positive or negative, respectively. (i.e. positive = LEFT-pad, negative = RIGHT-pad).
     */
    public static class TextWriter {
        protected final Writer writer;
        public TextWriter(final Writer writer) { this.writer = writer; }
        public TextWriter write(String s) { try { writer.write(s); } catch (IOException ex){ ex.printStackTrace(); } return this; }
        public TextWriter write(char c) { return write(Character.toString(c)); }
        public TextWriter write(int value) { return write(Integer.toString(value)); }
        public TextWriter write(Number value) { return write(value.toString()); }
        /** Write a value padded with spaces on the left (if padding is positive) or right (if padding is negative). */
        public TextWriter write(int padding, String value) { return write(pad(-padding, value)); } // Note: the Strings.pad function uses the opposite convention. Pad-characters are appended to the right side if padding is positive, or to the left if padding is negative.
        public TextWriter write(int padding, char c) { return write(padding, Character.toString(c)); }
        public TextWriter write(int padding, int value) { return write(padding, Integer.toString(value)); }
        public TextWriter write(int padding, Number value) { return write(padding, value.toString()); }
        /** Write a single space character. */
        public TextWriter space() { return write(' '); }
        /** Write the specified number of space characters. */
        public TextWriter space(int count) { return write(fromChar(' ', count)); }
        /** Write a system-specific line separator */
        public void writeln() { write(System.lineSeparator()); }
    }
    private static class CTEWriter extends TextWriter {
        public CTEWriter(final Writer writer) { super(writer); }
        //public CTEWriter write(String s) { return (CTEWriter) super.write(s);}
        public CTEWriter space() { return (CTEWriter) super.space();} // same as TextWriter.space, but return a CTEWriter
        /** Write a property value */
        public CTEWriter writep(String propName, int valuePadding, String value) { return (CTEWriter)space().write(propName).write(":").write(valuePadding, value); }
        public CTEWriter writep(String propName, int valuePadding, int value) { return writep(propName, valuePadding, Integer.toString(value)); }
        //public CTEWriter writep(String propName, int valuePadding, Object value) { return writep(propName, valuePadding, value.toString()); }
    }

    /**
     * Write an Extended CT (CTE) Drawing files.
     * A CTE is a CT file with additional information encoded in the columns past those in a normal CT.
     */
    public static void writeCteDrawingFile(RnaSceneGroup scenes, String path, boolean append) throws IOException {
        try(Writer writer = newUTF8Writer(path, append)) {
            writeCteDrawingFile(scenes, writer, true);
        }
    }
    public static void writeCtFile(RnaSceneGroup scenes, String path, boolean append) throws IOException {
        try(Writer writer = newUTF8Writer(path, append)) {
            writeCteDrawingFile(scenes, writer, false);
        }
    }
    /**
     * Write an Extended CT (CTE) Drawing files.
     * A CTE is a CT file with additional information encoded in the columns past those in a normal CT.
     */
    public static void writeCteDrawingFile(RnaSceneGroup scenes, Writer output, boolean includeExtendedProperties) {
            CTEWriter w = new CTEWriter(output);
            for (RnaScene scene : scenes) {
                int strands = scene.strands.size();
                int total = scene.getNucCount() + 3*(strands-1);
                // first line header
                w.write(CTE.PAD_TOTAL, total).space(2).write(scene.title).writeln(); // write the total number of bases and the structure label.
                // write nucleobases
                int index = 1; // first base is #1
                int[] strandIndex = scene.buildStrandIndex(1, RnaFileIO.InterMolecularBreakLength);
                for(Strand strand : scene.strands) {
                    if (index!=1) // If this is the second (or subsequent) strand, insert the inter-molecular linker (e.g. "III")
                        for (int i = 0; i < RnaFileIO.InterMolecularBreakStr.length(); i++) {
                            // Since this is a linker, we know it will never be the first or last nucleotide, so no index checks are necessary.
                            w.write(CTE.PAD_INDEX+1, index).space().write(CTE.PAD_SYMBOL, RnaFileIO.InterMolecularBreakChar); // write the index and symbol
                            w.space().write(CTE.PAD_INDEX, index-1).space().write(CTE.PAD_INDEX, index+1);  // write the previous and next index
                            w.space().write(CTE.PAD_INDEX, 0).space().write(CTE.PAD_INDEX, index); // write the pair and historical/annotative index
                            w.writeln();
                            index++;
                        }

                    for (Nuc n : strand) {
                        w.write(CTE.PAD_INDEX+1, index);
                        w.space().write(CTE.PAD_SYMBOL, isEmpty(n.symbol)?'?':n.symbol.charAt(0));
                        w.space().write(CTE.PAD_INDEX, index==1?0:index-1); // prev index
                        w.space().write(CTE.PAD_INDEX, index==total?0:index+1); // next index
                        w.space().write(CTE.PAD_INDEX, n.isPaired()?n.getPaired().indexInScene(strandIndex):0); // pair index
                        w.space().write(CTE.PAD_INDEX, n.number==-1?index:n.number); // historial index

                        if (includeExtendedProperties) {
                            w.space().write(CTE.PROP_SEPARATOR); // write comment character to indicate start of extended properties.
                            // Write "extended" nucleobase properties
                            w.writep(CTE.X, CTE.PAD_LOC, Math.round(n.location.x)); //write x location
                            w.writep(CTE.Y, CTE.PAD_LOC, Math.round(n.location.y));
                            Bond b = n.getPairBond();
                            if (b != null && b.type != BondType.Default) w.writep(CTE.BOND_TYPE, 1, b.type.getAbbrev());
                            NucStyle s = n.style;
                            if (s != null && !s.isEmpty()) {
                                if (s.textColor != null)
                                    w.writep(CTE.TEXT_COLOR, CTE.PAD_COLOR, Colors.getName(s.textColor));
                                if (s.fillColor != null)
                                    w.writep(CTE.FILL_COLOR, CTE.PAD_COLOR, Colors.getName(s.fillColor));
                                if (s.outlineColor != null)
                                    w.writep(CTE.LINE_COLOR, CTE.PAD_COLOR, Colors.getName(s.outlineColor));
                                if (s.bondColor != null)
                                    w.writep(CTE.BOND_COLOR, CTE.PAD_COLOR, Colors.getName(s.bondColor));
                                if (s.font != null)
                                    w.writep(CTE.FONT, CTE.PAD_FONT, s.font.toString());
                                if (s.numberOffset != 0F)
                                    w.writep(CTE.LBL_OFFSET, CTE.PAD_LOC, Float.toString(s.numberOffset));
                            }
                        }
                        w.writeln();
                        index++;
                    }
                }
            }
    }

    /**
     * This replacement for Scanner automatically skips whitespace, but does NOT skip over newlines unless
     * {@link #nextLine()} is called.
     *
     * The primary reason that this class was created (and that Scanner was unacceptable) was that
     * if Scanner's delimiter were set to just spaces and tabs, but not newlines, it would fail to recognize
     * valid content at the end of a line. for example getNextInt() would fail for "3\n".
     * This class allows testing (and extracting the results of) Patterns directly against the
     * content of the current line (at the current read position) without regard to delimiters, whereas Scanner
     * seems to only allow the testing of Patterns against pre-delimited text. For example:
     * Text: "Hello there!"
     * {@code Scanner(text).hasNext("Hello There") => false} (because the test is only applied to { "Hello", "There!" }[0]
     * {@code SimpleScanner(text).hasNext("Hello There") => true} (although whitespace BEFORE a match has already been skipped,
     * it can still be matched in the subsequent content.)
     * In general, Scanner's hasNextXXX functions are replaced by atXXX. {@link Scanner#hasNext()} is replaced by
     * a (negated) call to {@link #atEndOfLine()}.
     * This function also keeps track of the current line and column numbers so users can be provided with more detailed
     * information about where a potential error was found.
     */
    private static class SimpleScanner {
        private static final String END_TOKEN = "(?=\\s|$)";
        private static final Pattern
                SPACE= Pattern.compile("[\\s\\t]+"),
                NUM = Pattern.compile("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?"+END_TOKEN),
                BOOL = Pattern.compile("true|false", Pattern.CASE_INSENSITIVE),
                //PROP = Pattern.compile("(\\w+)[\\s\\t]*:"),
                ANY = Pattern.compile("[^\\s]+"+END_TOKEN);

        protected final BufferedReader reader;
        private String line;
        protected int lineNum, pos;
        protected Matcher matcher;
        protected boolean eof;

        public SimpleScanner(String content) { this (new BufferedReader(new StringReader(content)));  }
        public SimpleScanner(Reader input) {
            reader = input instanceof BufferedReader ? (BufferedReader)input : new BufferedReader(input);
        }

        public int lineNumber() { return lineNum; }
        public int colNumber() { return pos; }

        public boolean atNumber() { return matchNext(NUM);  }
        public boolean atBool() { return matchNext(BOOL);  }
        public boolean atEndOfLine() { return !matchNext(ANY); }
        public boolean atEOF() { return eof; }
        //public boolean atProperty() { return matchNext(PROP); }
        public boolean at(Pattern pattern) { return matchNext(pattern); }

        /** Returns the text of the next token, regardless of type, but does not advance the read position. */
        public String peek() { return matchNext(ANY)?matcher.group():null; }
        /** Returns the text of the next token, regardless of type, and advance the read position. */
        public String next() { String ret = peek(); accept(); return ret; }
        //public String getProperty() { require(PROP); accept(); return matcher.group(1); }
        public String next(Pattern pattern) { require(pattern); String s = matcher.group(); accept(); return s; }
        public boolean nextBool() { return Boolean.parseBoolean(next(BOOL)); }
        public double nextDouble() { double d = Double.parseDouble(require(NUM)); accept(); return d; }
        public float nextFloat() { float f = Float.parseFloat(require(NUM)); accept(); return f; }
        public int nextInt() { int ret = Integer.parseInt(require(NUM)); accept(); return ret; }
        /** Returns a Matchresult if the pattern was matched, or null otherwise. */
        public MatchResult peekMatch(final Pattern pattern) { return matchNext(pattern)?matcher:null; }
        public MatchResult nextMatch(final Pattern pattern) {
            MatchResult m = peekMatch(pattern);
            accept(); // does nothing if m is null.
            return m;
        }

        /**
         * Returns the remaining content on the current line (i.e. the substring after the current read position).
         * Also advances the read position to the next line, so that subsequent queries (atXXX) and extractions (getXXX)
         * will be applied to the next line.
         */
        public String nextLine() {
            if (line == null && !eof)
                readLine();
            String ret = line == null || pos == 0 ? line : pos == line.length() ? "" : line.substring(pos);
            line = null; // will force readline on next request.
            return ret;
        }


        // Non-public members
        protected boolean readLine() {
            if (!eof) try {
                if (null != (line = reader.readLine())) {
                    pos=0; lineNum++;
                    return true;
                }
                eof = true; // i.e. if line==null
            } catch (IOException ex) { ex.printStackTrace(); }
            return false;
        }

        protected String require(final Pattern pattern) {
            if  (matchNext(pattern)) return matcher.group();
            throw new NoSuchElementException();
        }

        // advance the read position past the last match.
        protected void accept() {
            if (matcher != null) pos = matcher.end();
            matcher = null;
        }
        protected boolean matchNext(Pattern p) {
            // If this pattern was just matched (and is still cached) return true.
            if (matcher != null && matcher.pattern() == p) return true;
            if (line == null && !readLine()) return false; // read the next line if necessary

            // If the next thing is whitespace, skip it.
            matcher = SPACE.matcher(line).region(pos, line.length());
            if (matcher.lookingAt())
                pos = matcher.end();

            // Now search for the requested pattern
            matcher = p.matcher(line).region(pos, line.length());
            if (matcher.lookingAt())
                return true;
            // No match
            matcher = null;
            return false;
        }
    }
    public static RnaSceneGroup readCteDrawingFile(String path) throws IOException, SyntaxErrorException {
        try(BufferedReader reader = Files.newBufferedReader(Paths.get(path), StandardCharsets.UTF_8)) {
            return readCteDrawingFile(reader);
        }
    }
    public static RnaSceneGroup readCteDrawingFile(Reader input) throws SyntaxErrorException {
        // read a CT file with extended information.
        RnaSceneGroup scenes = new RnaSceneGroup();
        scenes.setSource("<input>", FileType.CteDrawing);
        Pattern LINE_COMMENT = Pattern.compile("[;#]"), // allow comments. reserve for future use.
                PROP_NAME = Pattern.compile("(\\w+)[\\s\\t]*:"),
                PROP_START = Pattern.compile(CTE.PROP_SEPARATOR);
        SimpleScanner scanner = new SimpleScanner(input);
        int totalBases = 0;
        RnaScene scene = null;
        Strand strand = null;
        NucStyle style = null;
        int count = 0;
        while (!scanner.atEOF()) {
            // here we are at the start of a new line.
            if (scanner.at(LINE_COMMENT)||scanner.atEndOfLine()) {
                scanner.nextLine(); // skip the current line.
                continue;
            }
            if (scene == null) {
                try {
                    // First line of actual CT content. e.g.   "300   Some RNA Title"
                    totalBases = scanner.nextInt();
                } catch (NoSuchElementException ex) {
                    throw parseError(scanner, "Expected the number of nucleobases in the structure, but found '%s'.", scanner.next());
                }
                if (totalBases<0||totalBases>RnaFileIO.MAX_SEQUENCE_LENGTH)
                    throw parseError(scanner, "Invalid number of nucleobases: '%d' (expected from 1 to %d).", totalBases, RnaFileIO.MAX_SEQUENCE_LENGTH);
                scene = new RnaScene();
                scenes.add(scene);
                scene.title = scanner.nextLine().trim();
                if (scenes.size() == 0)
                    scenes.setTitle(scene.title); // if this is the first scene, set the overall title also
                strand = scene.firstStrand(); //Note: ALL nucleotides are added to the first (and only) strand. The strand will be split later at intermolecular linkers.
                strand.addPlaceholders(totalBases); // pre-create and add all nucleotides for improved performance.
                count = 0;
            } else {
                // read the next nucleobase
                String expected = null;
                try {
                    Nuc nuc = strand.get(count++); // Note: count will be incremented for the following messages, so it will be 1-based there, but 0-based here.

                    Bond bond = null;
                    if (style==null) style = new NucStyle();

                    expected = "index" ;
                    int index = scanner.nextInt();
                    if (index != count)
                        throw parseError(scanner, "Invalid nucleobase index %d (expected %d)", index, count);

                    expected = "symbol" ;
                    nuc.symbol = scanner.next();

                    expected = "prev-index" ;
                    int prev = scanner.nextInt();
                    if (prev != 0 && prev != count - 1)
                        throw parseError(scanner, "Invalid nucleobase connection to previous base %d (expected 0 or %d)", prev, count - 1);
                    if (prev < 0 || prev > totalBases)
                        throw parseError(scanner, "Invalid nucleobase index %d (expected range: 1 to %d)", prev, totalBases - 1);

                    expected = "next-index" ;
                    int next = scanner.nextInt();
                    if (next != 0 && next != count + 1)
                        throw parseError(scanner, "Invalid nucleobase connection to next base %d (expected 0 or %d)", next, count + 1);

                    // TODO: add this test back in once we switch from exceptions to warnings.
//                    if (next < 0 || next > totalBases)
//                        throw parseError(scanner, "Invalid nucleobase index %d (expected range: 1 to %d)", next, totalBases);

                    expected = "pair-index" ;
                    int pair = scanner.nextInt();
                    if (pair < 0 || pair > totalBases)
                        throw parseError(scanner, "Invalid nucleobase index %d (expected range: 1 to %d)", pair, totalBases);

                    expected = "historic-index" ;
                    nuc.number = scanner.atNumber() ? scanner.nextInt() : index;

                    if (next==0 && count != totalBases)
                        strand = scene.strands.add();

                    if (pair != 0) {
                        if (pair > index) {
                            // this is a forward reference to a nuc that is defined later.
                            //pairs[pair-1] = bond = new Bond.BondInfo(index-1, pair-1, BondType.Default);
                            if (strand.get(pair-1).isPaired())
                                throw parseError(scanner.lineNumber(), -1, "Inconsistent base pairing: bases %d and %d are both paired to base %d: ", index, strand.get(pair-1).getPaired().indexInScene()+1, pair);
                            bond = scene.addBond(index -1 , pair-1, BondType.Default); // note that 'pair' could refer to a nucleobase on another strand. so for now, use '0' as the strand index and re-index the bonds later.
                        } else {
                            bond = nuc.getPairBond(); // the bond already exists. it was created by a previous nucleotide.
                            // pairs[index-1];
                            if (bond == null || bond.getNuc3() != nuc) // it should point to the already-created bond.
                                throw parseError(scanner.lineNumber(), -1, "Inconsistent base pairing information between bases %d and %d: ", pair, index);
                        }
                    } else if (nuc.isPaired()) // i.e. bond was created by a previous nuc, but this one has pair=0.
                        throw parseError(scanner.lineNumber(), -1, "Inconsistent base pairing information between bases %d and %d: ", nuc.getPaired().indexInScene()+1, index); // this nuc's index was referenced by N1, but this nuc's pair is 0.

                    // If we want to skip unexpected columns, we could use:  while (!scanner.at(PROP_START)) scanner.next();
                    // For now, we'll be more strict and require that either the props start here or the line ends.
                    if (scanner.at(PROP_START)) {
                        scanner.next(PROP_START);
                        while (scanner.at(PROP_NAME)) {
                            String prop = scanner.nextMatch(PROP_NAME).group(1);
                            switch (prop.toUpperCase()) {
                                case CTE.X:
                                    expected = "X-pos" ;
                                    nuc.location.x = scanner.nextFloat();
                                    break;
                                case CTE.Y:
                                    expected = "Y-pos" ;
                                    nuc.location.y = scanner.nextFloat();
                                    break;
                                case CTE.TEXT_COLOR: // C (same as TC)
                                    expected = "Color" ;
                                    style.textColor = Colors.getColor(scanner.next());
                                    break;
                                case CTE.FILL_COLOR: //"FC"
                                    expected = "Fill-Color" ;
                                    style.fillColor = Colors.getColor(scanner.next());
                                    break;
                                case CTE.LINE_COLOR: //"LC":
                                    expected = "Line-Color" ;
                                    style.outlineColor = Colors.getColor(scanner.next());
                                    break;
                                case CTE.BOND_COLOR: //"BC":
                                    expected = "Bond-Color" ;
                                    style.bondColor = Colors.getColor(scanner.next());
                                    break;
                                case CTE.FONT: //"FS":
                                    expected = "Font" ;
                                    style.font = Font.decode(scanner.next());
                                    break;
                                case CTE.LBL_OFFSET: //"LO":
                                    expected = "Label-Offset" ;
                                    style.numberOffset = scanner.nextFloat();
                                    break;
                                case CTE.BOND_TYPE: //"B":
                                    expected = "Bond-Type" ;
                                    if (bond == null) throw parseError(scanner, "Basepair properties are specified for base number "+index+", which is not paired.");
                                    bond.type = BondType.fromAbbrev(scanner.next());
                                    break;
                                default:
                                    System.err.println("Warning: unknown nucleobase property in CTE file: " + prop);
                                    break;
                            }
                        }
                        if (!style.isEmpty()) {
                            nuc.style = style;
//                            scene.addStyle(style);
//                            scene.setStyle(nuc, style);
                            style = null; //reset so a new one will be created for the next nuc
                        }
                    }
                    if (!scanner.atEndOfLine())
                        throw parseError(scanner, "Unexpected text after columns: \"%s\".", scanner.next());
                    scanner.nextLine();
                } catch (NoSuchElementException | NumberFormatException ex) {
                    throw parseError(scanner, "The nucleobase information is malformed. (Expected %s, but found '%s').", expected, scanner.atEndOfLine() ? "end-of-line" : scanner.next());
                }
                if (count == totalBases) {
                    // this is the end of the current structure
//                    for (Bond.BondInfo pair : pairs)
//                        if (pair != null)
//                            scene.addBond(pair);

                    // Restructure the scene to account for Intermolecular Linkers and/or strand-breaks
//                    if (strandBreaks.size()!= 0)
//                        for (Nuc n : strandBreaks)
//                            scene.divideStrand(n);
                    RnaFileIO.splitStrandsAtLinkers(scene);
                    scene = null; //reset so new scene will be created.
                }
            }
        } // while(!scanner.atEOF())
        if (scene!=null && scene.getNucCount()!=totalBases)
            throw parseError(scanner, "The file ended before the expected number of nucleobases were listed. (Expected %i bases, but found '%i').", totalBases, scene.getNucCount());
        return scenes;
    }

    private static SyntaxErrorException parseError(SimpleScanner scanner, String message, Object... args) { return parseError(scanner.lineNumber(), scanner.colNumber(), message, args); }
    private static SyntaxErrorException parseError(int line, int column, String message, Object... args) {
        String sb = "CTE File parsing error" +
                ": " +
                (args.length == 0 ? message : String.format(message, args));
        //        if (line!=-1)
//            sb.append(" on line ").append(line);
        return new SyntaxErrorException(sb, null, line, column);
    }
    //</editor-fold>

//    private StringBuilder currentLine = new StringBuilder();
//    private int sourceStartLine;
//
//    // Each content line can be continued onto the next line by appending a '\' character before the newline.
//    public void parseContentLine(String line, int sourceLineNumber) throws SyntaxErrorException {
//        line = line.trim();
//        if (currentLine.length() == 0) sourceStartLine = sourceLineNumber;
//        if (line.length() == 0) return;
//        if (line.charAt(line.length() - 1) == '\\') {
//            currentLine.append(line.substring(0, line.length() -1));
//            return;
//        } else if (currentLine.length() > 0) {
//            currentLine.append(line);
//            line = currentLine.toString();
//            currentLine.setLength(0);
//        }
//        parseFullLine(line, sourceStartLine, sourceLineNumber);
//    }
//    public void parseEnd() throws SyntaxErrorException {
//        if (currentLine.length() != 0)
//            throw new SyntaxErrorException("The command started on line " + sourceStartLine + " was never completed.", currentLine.toString(), sourceStartLine, 0);
//    }
//
//    private static Pattern NSD_LINE_HEADER_REGEX = Pattern.compile(
//            "\\s*(?<name>\\w+)\\s*(?:(?<str>\\d+)(?:\\.(?<seq>\\d+))?)?\\s*:\\s*",
//            Pattern.CASE_INSENSITIVE);
//
//    // This function parses "complete" lines -- i.e. not partial lines that are continued using '\'
//    public void parseFullLine(String line, int lineNumberStart, int lineNumberEnd) throws SyntaxErrorException {
//
//    }
}
