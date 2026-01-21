package ur_rna.RNAstructureUI.utilities;

import ur_rna.RNAstructure.backend.RNABackend;
import ur_rna.Utilities.PathTools;

/**
 *
 */
public enum FileType {
    CT(RNABackend.FILE_CT),
    SEQ(RNABackend.FILE_SEQ),
    PFS(RNABackend.FILE_PFS),
    SAV(RNABackend.FILE_SAV),
    DBN(RNABackend.FILE_DBN),
    SHAPE(RNABackend.FILE_DBN+1);

    private final int value;
    FileType(int backendValue) {
        value = backendValue;
    }
    public int backendValue() { return value; }
    public boolean equals(int value) { return value==this.value; }

    public static FileType guessTypeFromFilePath(String path) {
        String ext = PathTools.getExt(path, false, false);
        ext = ext.toUpperCase();
        switch (ext) {
            case "CT": return CT;
            case "SEQ":
            case "FA":
            case "FASTA":return SEQ;
            case "PFS":return PFS;
            case "SAV":
            case "FSV":
                return SAV;
            case "DBN":
            case "DOT":
            case "BRACKET":
            case "BRK": return DBN;
            default: return null; // unknown
        }
    }
}
