package ur_rna.GUITester.GuiTools;

public class GuiItemNotFoundException extends ur_rna.Utilities.ExtensibleException {
    private static final long serialVersionUID = 1L;

    public GuiItemNotFoundException() { }
    public GuiItemNotFoundException(String message) {
        super(message);
    }

    public GuiItemNotFoundException(GuiItemRef ref) {
        super("GUI item not found: " + ref);
    }

    public GuiItemNotFoundException(Throwable cause) {
        super(cause);
    }
    public GuiItemNotFoundException(String message, Throwable cause) {
        super(message, cause);
    }

    public GuiItemNotFoundException addRelative(GuiItemRef ref) {
        if (ref.getRelationship() != null)

        this.setMessage(
                String.format("When searching for %s component:\n->",
                        ref.getRelationship().name()) +
                        this.getMessage()
                );
        return this;
    }
}
