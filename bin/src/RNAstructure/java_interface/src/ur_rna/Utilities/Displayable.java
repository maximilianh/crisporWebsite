package ur_rna.Utilities;

import ur_rna.Utilities.annotation.ToFriendlyString;

@ToFriendlyString(method="HasDisplayString")
public interface Displayable {
    /**
     * Returns an informative representation of the object, which should indicate its type if it is not obvious.
     * Do NOT return {@code null} and avoid non-printable characters if possible.
     * @return An informative representation of the object.
     */
    String toDisplayString();
}
