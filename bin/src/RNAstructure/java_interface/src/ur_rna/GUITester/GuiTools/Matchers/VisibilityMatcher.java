package ur_rna.GUITester.GuiTools.Matchers;

import java.awt.*;

/**
 * Tests the visibility of the GUI component (whether it is showing or not)
 */
public class VisibilityMatcher extends ComposableBase {
    private boolean _vis = true;

    public VisibilityMatcher() { }
    public VisibilityMatcher(boolean visible) {
        _vis = visible;
    }

    @Override
    public boolean matches(Component c) { return c.isShowing() == _vis; }

    public boolean getVisible() { return _vis; }
    public void setVisible(boolean value) { _vis = value; }

    @Override
    public String toString() {
        return "{ VisibilityMatcher vis=" + (_vis ? "showing" : "hidden") + " }";
    }
}
