package ur_rna.StructureEditor;

import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.PathTools;

import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.util.*;

import static ur_rna.Utilities.Strings.isEmpty;

public enum FileType {
    // Misc RNAstructure file types
    /** Oligo report file. */
    Report("Report File", "rep"),
    /** SHAPE chemical modification file.     */
    SHAPE("SHAPE File", "shape"),
    /** Select an oligo list file. */
    OligoList("Sequence List File", "lis;olis"),
    /** Select a general OUT file. */
    OUT("OUT File", "out"),
    /** Select a partition function save file. */
    PartitionSav("Partition Function Save File", "pfs"),
    /** Select a Dynalign save file. */
    DynalignSav("Dynalign Save File", "dsv"),
    /** Select a folding save file. */
    FoldingSav("Folding Save File", "sav,save"),
    /** Alignment file. */
    Alignment("Alignment File", "ali"),
    /** Constraints file. */
    Constraints("Constraint File", "con"),
    /** Dot plot file. */
    DotPlot("Dot Plot File", "dp"),

    // Sequence types
    /** FASTA Sequence file */
    Fasta("FASTA File", "fasta;fa"),
    /** SEQ Sequence file */
    Seq("SEQ File", "seq"),
    /** Plain text Sequence file */
    SeqText("Sequence Text File", "txt"),
    /** GenBank Flat File Format  e.g.: https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html */
    SeqGenbank("Genbank Sequence File", "gen"),

    // Structure files
    /** CT (Connection Table) structure file. */
    CT("CT File", "ct"),
    /** Dot-Bracket file. */
    DotBracket("Bracket File", "dbn,bracket,dot"),
    /** Helix format. (encodes bonds, but not sequence.) */
    Helix("Helix File", "helix;hlx;htxt"), // probably should just be "txt", but I'm not sure if the Helix format matches the SeqText format, which is more deserving of the "txt" extension

    // Structure drawing formats.
    /** Native structure drawing format. Similar to JSON. Encodes Bases, Bonds, and Styles separately. */
    NsdDrawing("Structure Drawing File", Program.DrawingFileExtension),
    /** Extended CT format. Similar to CT, but with additional columns for nucleobase position, bond type, and color/font styles */
    CteDrawing("Extended CT Drawing File", "cte"),

    // Vector graphics files
    /** Scalable Vector Graphics (SVG) File */
    SVG("SVG - Scalable Vector Graphics File", "svg"),
    /** Postscript (PS) file. */
    Postscript("Postscript File", "ps"),

    // Image files
    PNG("PNG Image File", "png"),
    GIF("GIF Image File", "gif"),
    JPEG("JPEG Image File", "jpg,jpeg"),

    // Compound types:
    Sequence("Sequence File", Fasta, Seq, SeqText),
    //VectorImage ("Vector Image File", SVG, Postscript),
    VectorImage ("Vector Image File", SVG),
    Image ("Image File", PNG, GIF, JPEG),
    Structure("Structure File", CT, DotBracket),
    Drawing("Structure Drawing", NsdDrawing, CteDrawing),

    AnyOpenable("Drawing, Structure, or Sequence File", Drawing, Structure, Sequence, PartitionSav, FoldingSav),
    AnySaveable("Any Save Type", Drawing, Structure, Sequence, VectorImage, Image);

    public static final List<FileType> allTypes;
    public static final List<FileType> primaryTypes;
    public static final Map<String, FileType> extMap;
    static {
        FileType[] values = values();
        allTypes = Collections.unmodifiableList(Arrays.asList(values));
        List<FileType> primary = new ArrayList<>();
        Map<String, FileType> map = new HashMap<>();
        for (FileType f : allTypes) {
            if (f.subTypes().size()==0) {
                primary.add(f);
                for (String ext : f.getExtensions()) {
                    ext = ext.toLowerCase();
                    if (map.containsKey(ext))
                        // Each extension must be associated with exactly one PRIMARY FileType. However COMPOUND types can be composed of multiple PRIMARY types.
                        throw new IllegalStateException(String.format("The extension %1s is included by more than one primary %2s: %3s and %4s.",
                                ext, FileType.class.getSimpleName(),
                                map.get(ext).name(),
                                f.name()
                        ));
                    map.put(ext, f);
                }
            }
        }
        primaryTypes = Collections.unmodifiableList(primary);
        extMap = Collections.unmodifiableMap(map);
    }

    private final String description, extensions;
    private final Set<FileType> subTypes;

    FileType(String description, FileType... subTypes) {
        this(description, "", new LinkedHashSet<>(Arrays.asList(subTypes)));
    }
    FileType(String description, String extensions) {
        this(description, extensions.trim().replace(',', ';').replace(' ', ';'), Collections.emptySet());
    }
    FileType(String description, String extensions, Set<FileType> subTypes) {
        this.description = description;
        this.extensions = extensions;
        this.subTypes = subTypes;
        String invalidState = null;
        if (extensions!=null&&!extensions.isEmpty()&&!subTypes.isEmpty()) invalidState = "both";
        else if ((extensions==null||extensions.isEmpty())&&subTypes.isEmpty()) invalidState = "neither";
        if (invalidState != null)
            throw new IllegalStateException(String.format("A %s must either contain primary file extensions OR be composed of sub-types, but %s has %s, which is invalid.", FileType.class.getSimpleName(), this.name(), invalidState));
    }
    public String[] getExtensions() {
        return extensions==null || extensions.isEmpty() ? ObjTools.EMPTY_STRING_ARRAY : extensions.split(";");
    }
    public String firstExtension() {
        return extensions==null || extensions.isEmpty() ? "" : extensions.split(";")[0];
    }
    public String getFilterString() { return this.getFilterString(false); }
    public String getFilterString(boolean splitSubTypes) {
        /*
         Example:
            Fasta("FASTA File", "fasta"),
            Seq("SEQ File", "seq"),
            Sequence("Sequence File", Fasta, Seq),

            Fasta    => "FASTA File|fasta"         (regardless of splitSubTypes)
            Seq      => "SEQ File|seq"             (regardless of splitSubTypes)
            Sequence => "Sequence File|fasta;seq"       (splitSubTypes = false)
            Sequence => "FASTA File|fasta|SEQ File|seq" (splitSubTypes = true)

         */

        StringBuilder sb = new StringBuilder();
        if (splitSubTypes) {
            addFilterString(sb); // ends with an extra '|'
        } else {
            sb.append(description).append('|');
            addExtensions(sb); // ends with an extra ';'
        }
        sb.setLength(sb.length() - 1); //remove the last '|' or ';'
        return sb.toString();
    }

    public void addFilterString(StringBuilder sb) {
        // Add "description|exts|" but ONLY if it has at least one extension.
        if (!isEmpty(extensions))
            sb.append(description).append('|').append(extensions).append('|');
        if (subTypes.size() != 0) {
            for(FileType t : subTypes)
                t.addFilterString(sb);
        }
    }
    private void addExtensions(StringBuilder sb) {
        if (!isEmpty(extensions))
            sb.append(extensions).append(';');
        for(FileType t : subTypes)
            t.addExtensions(sb);
    }
    public String getDescription() {
        return description;
    }
    public Set<FileType> subTypes() { return subTypes; }

    public boolean includes(FileType childType) {
        if (this == childType || subTypes().contains(childType))
            return true;
        for (FileType sub : subTypes())
            if (sub.includes(childType))
                return true;
        return false;
    }
    public boolean isSubTypeOf(final FileType... parentTypes) {
        for(FileType t : parentTypes)
            if (t.includes(this))
                return true;
        return false;
    }
    public static FileType fromExt(String ext) {
        if (ext.isEmpty() || ext.equals(".")) return null;
        if (ext.charAt(0)=='.') ext = ext.substring(1);
        return extMap.get(ext.toLowerCase());
    }
    public static FileType fromFilter(final FileFilter filter) {
        if (filter instanceof FileNameExtensionFilter) {
            for(String ext : ((FileNameExtensionFilter)filter).getExtensions()) {
                FileType f = fromExt(ext);
                if (f != null) return f;
            }
        }
        return fromDesc(filter.getDescription());
    }
    public static FileType fromDesc(String description) {
        int pos = description.indexOf(" (*.");
        if (pos != -1)
            description = description.substring(0, pos);
        description = description.trim();
        for(FileType f : allTypes)
            if (f.getDescription().equals(description))
                return f;
        return null;
    }

    public static FileType fromName(final String typeName) {
        if (typeName != null && !typeName.isEmpty())
            for (FileType f : allTypes)
                if (f.name().equals(typeName))
                    return f;
        return null;
    }
    public static FileType fromAny(final String typeName, final String extension, final String description, final String filePath, final FileType defaultValue) {
        FileType f = fromName(typeName);
        if (f != null) return f;
        if (extension != null)
            f = fromExt(extension);
        if (f != null) return f;
        if (description!=null)
            f = fromDesc(description);
        if (f != null) return f;
        if (filePath!=null)
            f = fromExt(PathTools.getExt(filePath));
        if (f != null) return f;
        return defaultValue;
    }

    @Override
    public String toString() { return getDescription(); }
}
