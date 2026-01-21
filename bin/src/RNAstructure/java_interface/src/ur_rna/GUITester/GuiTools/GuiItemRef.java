package ur_rna.GUITester.GuiTools;

import abbot.finder.Matcher;
import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.annotation.ToFriendlyString;

import java.awt.*;

@ToFriendlyString
public class GuiItemRef {
    public final Matcher matcher;
    public Component lastFound;
    private GuiItemRef relative;
    private GuiRelative relationship;
    public Component getLastFound() {
        return lastFound;
    }
    public GuiItemRef getRelative() {
        return relative;
    }
    public GuiRelative getRelationship() {
        return relationship;
    }
    public GuiItemRef(final Matcher matcher) {
        this.matcher = matcher;
    }

    /**
     * Returns a GuiItemRef that always matches a single known component.
     * @param identity The component to match
     */
    public GuiItemRef(final Component identity) {
        this.lastFound = identity;
        this.matcher = new ComponentMatcher(identity);
    }

    /**
     * Returns a Matcher that always matches a single known component.
     */
    public static class ComponentMatcher implements Matcher {
        public final Component identity;
        public ComponentMatcher(Component identity) {
            this.identity = identity;
        }
        @Override
        public boolean matches(final Component c) {
            return c == identity;
        }
        @Override
        public String toString() {
            return "{ Component Matcher: " + identity + " }";
        }

    }

    public GuiItemRef clearFound(final boolean clearRelatives) {
        lastFound = null;
        if (clearRelatives && relative != null)
            relative.clearFound(true);
        return this;
    }

    @Override
    public String toString() {
        return toString(true);
    }

    public String toString(final boolean includeFullRelativeDescription) {
        return String.format(
                "{ GUI-Reference Criteria: %s %s }",
                ObjTools.toStr(matcher, "none"),
                relativeToString(includeFullRelativeDescription)
        );
    }

    public String relativeToString(boolean includeFullRelativeDescription) {
        if (relationship == null)
            return "";
        return relationship.name() + " of " + (
                includeFullRelativeDescription ?
                        this.relative :
                        "another component"
        );
    }

    public void setRelative(final GuiItemRef relative, final GuiRelative relationship) {
        this.relative = relative;
        this.relationship = relationship;
    }
}
