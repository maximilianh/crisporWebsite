package ur_rna.GUITester.GuiTools.Matchers;

import abbot.util.Regexp;
import ur_rna.Utilities.Convert;

import javax.accessibility.AccessibleContext;
import javax.swing.*;
import javax.swing.text.JTextComponent;
import java.awt.*;
import java.util.Map;

public class DescriptionMatcher extends ComposableBase {
    private String _find = null;
    private int scope;

    public final static int SEARCH_NONE = (0);
    public final static int SEARCH_TRUE_NAME = (1<<0);
    public final static int SEARCH_TRUE_INPUT = (1<<1);
    public final static int SEARCH_TRUE_CAPTION = (1<<2);
    public final static int SEARCH_TRUE_TEXT = SEARCH_TRUE_CAPTION | SEARCH_TRUE_INPUT;
    public final static int SEARCH_TOOLTIP = (1<<10);
    public final static int SEARCH_SELECTION = (1<<14);
    public final static int SEARCH_ACCESSIBLE_NAME = (1<<20);
    public final static int SEARCH_ACCESSIBLE_TEXT = (1<<21);
    public final static int SEARCH_ACCESSIBLE_DESC = (1<<22);

    public final static int SEARCH_ACCESSIBLE = SEARCH_ACCESSIBLE_NAME | SEARCH_ACCESSIBLE_TEXT | SEARCH_ACCESSIBLE_DESC;
    public final static int SEARCH_NAME = SEARCH_TRUE_NAME | SEARCH_ACCESSIBLE_NAME | SEARCH_ACCESSIBLE_DESC;
    public final static int SEARCH_INPUT = SEARCH_TRUE_INPUT | SEARCH_ACCESSIBLE_TEXT;
    public final static int SEARCH_CAPTION = SEARCH_TRUE_CAPTION | SEARCH_ACCESSIBLE_DESC;
    public final static int SEARCH_TEXT = SEARCH_INPUT | SEARCH_CAPTION;
    public static final int SEARCH_COMMENT = SEARCH_TOOLTIP | SEARCH_ACCESSIBLE_DESC;

    public final static int SEARCH_ALL = -1;

    public DescriptionMatcher(String text) {
    	this(text, SEARCH_ALL);
    }

    public DescriptionMatcher(String text, int searchScope) {
    	this._find = text;
        this.scope = searchScope;
    }    

    @Override
    public boolean matches(Component c) {
    	if (0 != (scope & SEARCH_TRUE_NAME) &&
                isMatch(c.getName())) return true;

        if (0 != (scope & SEARCH_TOOLTIP) &&
   		    c instanceof JComponent &&
                isMatch(((JComponent) c).getToolTipText())) return true;

        if (0 != (scope & SEARCH_TRUE_INPUT)) {
            if (c instanceof JTextComponent && isMatch(((JTextComponent) c).getText())) return true;
            if (c instanceof TextComponent && isMatch(((TextComponent) c).getText())) return true;
        }
        if (0 != (scope & SEARCH_TRUE_CAPTION)) {
            if (c instanceof AbstractButton && isMatch(((AbstractButton) c).getText())) return true; // includes JButton, JMenuItem
            if (c instanceof Button && isMatch(((Button) c).getLabel())) return true;
            if (c instanceof JLabel && isMatch(((JLabel) c).getText())) return true;
            if (c instanceof Label && isMatch(((Label) c).getText())) return true;
            if (c instanceof Frame && isMatch(((Frame) c).getTitle())) return true;
            if (c instanceof Dialog && isMatch(((Dialog) c).getTitle())) return true;
            if (c instanceof JMenuItem && isMatch(getMenuPath((JMenuItem) c))) return true;
        }

        if (0 != (scope & SEARCH_SELECTION)) {
            if (c instanceof JList && isMatch("" + ((JList) c).getSelectedValue())) return true;
            if (c instanceof JComboBox && isMatch("" + ((JComboBox) c).getSelectedItem())) return true;
        }

        if (0 != (scope &  SEARCH_ACCESSIBLE)){
            AccessibleContext ac = c.getAccessibleContext();
            if (ac != null) {
                if (0 != (scope & SEARCH_ACCESSIBLE_NAME) &&
                    isMatch(ac.getAccessibleName())) return true;

                if (0 != (scope & SEARCH_ACCESSIBLE_DESC) &&
                    isMatch(ac.getAccessibleDescription())) return true;

//                if (0 != (scope & SEARCH_ACCESSIBLE_TEXT) &&
//                    isMatch("" + ac.getAccessibleText())) return true;
            }
        }
		return false;
    }

    boolean isMatch(String actual) {
        return stringsMatch(_find, actual);
    }

    /**
     * Get a list of all text descriptions associated with this Component.
     * This includes the component's name as well as input text, captions,
     * labels, list-selections, tool-tips, etc.
     *
     * This function is not used during GUI component matching. Instead it
     * may be useful for tools that provide information to the user about
     * components (e.g. during GUI-robot script authoring).
     *
     * @param c
     * @param props
     */
    public static void getDescriptions(Component c, Map<String,String> props) {
        props.put("Name", c.getName());
        if (c instanceof JComponent) props.put("ToolTip", ((JComponent) c).getToolTipText());

        if (c instanceof JTextComponent) props.put("InputText", ((JTextComponent) c).getText());
        else if (c instanceof TextComponent) props.put("InputText", ((TextComponent) c).getText());

        if (c instanceof AbstractButton ) props.put("Caption", (((AbstractButton) c).getText())); // includes JButton, JMenuItem
        else if (c instanceof Button) props.put("Caption", (((Button) c).getLabel()));
        else if (c instanceof JLabel) props.put("Caption", (((JLabel) c).getText()));
        else if (c instanceof Label) props.put("Caption", (((Label) c).getText()));
        else if (c instanceof Frame) props.put("Caption", (((Frame) c).getTitle()));
        else if (c instanceof Dialog) props.put("Caption", (((Dialog) c).getTitle()));
        if (c instanceof JMenuItem) props.put("Menu Text", (getMenuPath((JMenuItem) c)));

        if (c instanceof JList) props.put("SelItem", ("" + ((JList) c).getSelectedValue()));
        else if (c instanceof JComboBox) props.put("SelItem", ("" + ((JComboBox) c).getSelectedItem()));

        AccessibleContext ac = c.getAccessibleContext();
        if (ac != null) {
            props.put("AccessName", ac.getAccessibleName());
            props.put("AccessDesc", ac.getAccessibleDescription());
        }
    }

//    static boolean stringsMatch(String expected, String actual) {
//        return ExtendedComparator.stringsMatch(expected, actual);
//    }

    /**
     * Performs a fuzzy string match that can be altered with various
     * prefix symbols, called modifiers. By default (no modifiers) this
     * performs a case insensitive string comparison that allows wildcards
     * (specifically POSIX Standard Shell globbing patterns) such as:
     * <ul>
     *     <li>{@code ?} (any character)</li>
     *     <li>{@code *} (zero or more characters)</li>
     *     <li>{@code [..]} (character range)</li>
     *     <li>{@code {..,..}} (alternatives)</li>
     *</ul>
     * If the standard string matching is not sufficient, a different mode can
     * be selected by prepending one of the following prefixes to the search string:
     * <dl>
     *     <dt>{@code ~:} (tilde, colon)</dt>
     *     <dd>Perform a regular expression match using {@link gnu.regexp.RE}.
     *       Multiline matches are enabled by ~:(?m)regexp
     *       Embedded newlines ("\n") in the match string will then match
     *       end-of-lines.
     *     </dd>
     *     <dt>{@code ==} (double equals)</dt>
     *     <dd>Perform strict match (case sensitive, no wildcards)</dd>
     *     <dt>{@code ~=} (tilde, equals)</dt>
     *     <dd>Perform a case-sensitive match allowing wildcards. (i.e the same
     *     type of fuzzy matching as the default, except that character case matters.
     *     </dd>
     *     <dt>{@code ~~} (double tilde)</dt>
     *     <dd>Perform a fuzzy match (the default type). Use this if
     *     the search expression could start with one of the other prefixes.
     *     e.g. {@code ~~==? } would perform a fuzzy search for two equals signs
     *     followed by any other single character. (i.e. {@code ==? })
     *     </dd>
     *</dl>
     */
    // Requiring exact matches eliminates the need to always include the start
    // "^" and finish "$" symbols.
    public static boolean stringsMatch(String pattern, String actual) {
        if (pattern == null)
            return actual == null;

        if (pattern.startsWith("~:")) // regular expression
            return Regexp.stringMatch(pattern.substring(2), actual);

        if (pattern.startsWith("==")) // strict equality
            return pattern.equals(actual);

        if (pattern.startsWith("~=")) // case-sensitive fuzzy match
            return Regexp.stringMatch(Convert.globToRegex(pattern), actual);

        if (pattern.startsWith("~~")) // case-INsensitive fuzzy match (just remove prefix)
            pattern = pattern.substring(2);

        // case-INsensitive fuzzy match (default if no prefix)
        return Regexp.stringMatch(Convert.globToRegex(pattern.toLowerCase()), actual==null ? "" : actual.toLowerCase());
    }

//    public boolean matchesMenu(JMenuItem mi) {
//        return  stringsMatch(_find, getMenuPath(mi));
//    }
//
    public static String getMenuPath(JMenuItem item) {
        Object parent = item.getParent();
        if(parent instanceof JPopupMenu) {
            parent = ((JPopupMenu)parent).getInvoker();
        }
        return parent instanceof JMenuItem ? getMenuPath((JMenuItem) parent) + "|" + item.getText() : item.getText();
    }
//
//    public static java.util.List splitMenuPath(String path) {
//        int lastFoundIndex = -1;
//        ArrayList selectionPath;
//        for(selectionPath = new ArrayList(); (lastFoundIndex = path.indexOf(124, lastFoundIndex)) != -1; ++lastFoundIndex) {
//            selectionPath.add(path.substring(0, lastFoundIndex));
//        }
//        selectionPath.add(path);
//        return selectionPath;
//    }

	public String toString() {
        return "[ Description Matcher:  \"" + _find + "\" (in " + scopeToString(scope) + ")]";
    }

    public static String scopeToString(final int scope) {
        switch (scope) {
            case SEARCH_NONE:
                return "nowhere";
            case SEARCH_ALL:
                return "any field";
            case SEARCH_TRUE_NAME:
                return "true name";
            case SEARCH_TRUE_INPUT:
                return "true input";
            case SEARCH_TRUE_CAPTION:
                return "true caption";
            case SEARCH_TRUE_TEXT:
                return "true text";
            case SEARCH_TOOLTIP:
                return "tooltip";
            case SEARCH_SELECTION:
                return "selection";
            case SEARCH_ACCESSIBLE_NAME:
                return "accessible name";
            case SEARCH_ACCESSIBLE_TEXT:
                return "accessible text";
            case SEARCH_ACCESSIBLE_DESC:
                return "accessible desc";
            case SEARCH_ACCESSIBLE:
                return "accessible";
            case SEARCH_NAME:
                return "name";
            case SEARCH_INPUT:
                return "input";
            case SEARCH_CAPTION:
                return "caption";
            case SEARCH_TEXT:
                return "text";
            case SEARCH_COMMENT:
                return "tooltip or desc";
            default:
                return Integer.toString(scope);
        }

    }
}
