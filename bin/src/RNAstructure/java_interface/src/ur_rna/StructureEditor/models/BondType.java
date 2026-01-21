package ur_rna.StructureEditor.models;

import java.util.HashMap;
import java.util.Map;

/**
 * These Bond types reflect additional user-input that could not be garnered directly from the
 * bases on either side of the bonds.
 * For example, a user could mark bonds as Forced, Standard, or Prohibited.
 *
 * Regarding certain operations, such as automatic layout or user interaction with helices etc, any bonds marked as
 * {@link #Prohibited} or {@link #Pseudo} may be ignored or handled differently than other bonds,
 * depending on the circumstances of the operation.
 *
 */
public enum BondType {
    /** This bond is likely to be paired in the given structure. */
    Standard("BP", "Standard Basepair"),
    /** User (or external evidence etc) indicates this bond MUST be paired. */
    Forced("CF", "Forced Pair (Constraint)"),
    /** User (or external evidence etc) indicates this bond must NOT be paired. */
    Prohibited("CP", "Prohibited Pair (Constraint)"),
    /** This bond is likely to be paired, but should be highlighted as special or unusual for some reason. */
    Special("NC", "Non-Canonical/Special Bond"),
    /** Pseduo-knot or quaternary bond. */
    Pseudo("PK", "Pseudo-Knot Bond");

    public static final BondType Default = Standard;
    private final String description;

    private String abbrev;
    private static Map<String, BondType> _lookup;
    BondType(String abbrev, String desc) {this.abbrev = abbrev; this.description = desc;}
    public String getAbbrev() { return  abbrev; }
    /**
     * Gets the BondType corresponding to an abbreviation.
     * @throws IllegalArgumentException if the abbreviation is empty or does not correspond to a known BondType.
     * */
    public static BondType fromAbbrev(String abbrev) {
        BondType found = fromAbbrev(abbrev, null);
        if (found == null)
            throw new IllegalArgumentException("Invalid BondType abbreviation: " + abbrev + ".");
        return found;
    }
    /**
     * Gets the BondType corresponding to an abbreviation.
     * Returns the value passed in as 'defaultIfInvalid' if the abbreviation is empty or does not correspond
     * to a known BondType.
     * */
    public static BondType fromAbbrev(String abbrev, BondType defaultIfInvalid) {
        if (abbrev == null || abbrev.isEmpty())
            return defaultIfInvalid;
        if (_lookup == null) {
            _lookup = new HashMap<>();
            for(BondType b : BondType.values()) {
                _lookup.put(b.abbrev, b);
            }
        }
        return _lookup.getOrDefault(abbrev, defaultIfInvalid);
    }
    public String getDescription() {
        return description;
    }
}
