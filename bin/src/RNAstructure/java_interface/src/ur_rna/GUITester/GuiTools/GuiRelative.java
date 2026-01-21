package ur_rna.GUITester.GuiTools;

/**
 * Describes the relationship between two GUI Components.
 * Relationships are directional, so it is important to note the order of the components.
 * Given components A and B, the relationship {@code A --[RELATION]--> B} is interpreted as "B is a RELATION of A".
 * e.g. {@code A --[Child]--> B}  means that B is the Child of A. Similarly {@code A --[PrevSibling]--> B} means that B is the previous sibling of A.
 *
 * In regards to the {@link GuiItemRef#relationship } property of a {@link GuiItemRef}, it should be
 * interpreted as GuiItemRef (subject) is a [RELATION] of {@link GuiItemRef#relative }
 *
 * So for example, if C and P are GuiItemRef objects and C.relative=P and C.relationship={@link GuiRelative#Child}
 * it indicates that C is a Child of P.
 *
 */
public enum GuiRelative {
    Sibling,
    Parent,
    Child,
    Ancestor,
    Descendant,
    LabelTarget
}
