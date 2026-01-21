package ur_rna.GUITester.GuiTools;

import abbot.finder.*;
import ur_rna.GUITester.GuiTools.Matchers.DescriptionMatcher;
import ur_rna.Utilities.AppLog;

import javax.swing.*;
import java.awt.*;

/**
 * @author Richard M. Watson
 */
public class GuiFinder extends BasicFinder {
    public GuiAppManager _owner;
    private AppLog _log;
    public GuiFinder(GuiAppManager owner) { _owner = owner; _log = _owner.getLog(); }
    /**
     * Find a component. Returns null if no match is found instead of throwing a ComponentSearchException.
     * @param m The matcher to use to find the component.
     * @return Returns the matching component if found, or null if no match is found.
     */
    public Component find(Matcher m) {
        return find((Container)null, m);
    }

    /**
     * Find a component using a {@link DescriptionMatcher}. Returns null if no match is found instead of throwing a ComponentSearchException.
     * @param descriptor The text to search for. Many different textual fields of components are compared with this. Encapsulate regular expressions inside slashes (e.g. /[a-z]+/ )
     * @return Returns the matching component if found, or null if no match is found.
     */
    public Component find(String descriptor, GuiRelative relationship, Component relative) throws ComponentSearchException  {
        return this.find(new DescriptionMatcher(descriptor), relationship, relative);
    }

    /**
     * Find a component with the built-in BasicFinder. Returns null if no match is found instead of throwing a ComponentSearchException.
     * @param relative a reference component that has a specific relationship with the desired target component.
     * @param matcher The matcher to use to find the component.
     * @return Returns the matching component if found, or null if no match is found.
     */
    public Component find(Matcher matcher, GuiRelative relationship, Component relative) {
        // if there is no relationship, find based only on the matcher.
        if (relationship == null && relative == null) return find(matcher);

        // if there is no relationship, but there IS a relative, use the default relationship.
        if (relationship == null) relationship = GuiRelative.Child; //  relationship==null, but relative != null, so set default relationship

        // There is a required relative, but it evaluated to null, so we cannot find the component.
        if (relative == null) // note that (relationship != null) or we would have returned already.
            return null; //reference component was not found, so we cannot search for its relative.

        switch (relationship) {
            case Child:
                return find((Container) relative, matcher); // TODO: Limit search scope to just the immediate children. (Currently this searches all descendants. May need to modify 'hierarchy' of the finder)
            case Parent:
                return relative.getParent();
            case Sibling:
                return find(relative.getParent(), matcher); // TODO: Limit search scope to just the immediate siblings. (Currently this searches all descendants of the parent.)
            case Ancestor:
                Container c;
                while (null != (c = relative.getParent()))
                    if (matcher.matches(c))
                        return c;
                return null;
            case Descendant:
                return find((Container) relative, matcher);
            case LabelTarget:
                if (relative instanceof JLabel)
                    return ((JLabel) relative).getLabelFor();
                return null;
            default:
                throw new IllegalArgumentException("Unknown GUI component relationship: " + relationship);
        }
    }

    @Override
    public Component find(Container root, Matcher matcher) {
        //if (root == null) root = _root;
        try {
            return super.find(root, matcher);
        } catch (ComponentNotFoundException ex) {
            return null;
        } catch (MultipleComponentsFoundException ex) {
            _log.error("Multiple Components Found Exception", ex);
            return null;
        }
    }

    public Component find(GuiItemRef ref, boolean required) throws GuiItemNotFoundException {
        if (ref.lastFound != null) return ref.lastFound;
        Component related = null;
        if (ref.getRelative() != null) {
            try {
                related = find(ref.getRelative(), required);
            } catch (GuiItemNotFoundException ex) {
                ex.addRelative(ref);
                throw ex;
            }
            if (related == null) return null; //relative was not found.
        }
        ref.lastFound = find(ref.matcher, ref.getRelationship(), related);
        if (required && ref.lastFound == null) {
            _owner.printHierarchyTrace();
            throw new GuiItemNotFoundException(ref);
        }
        return ref.lastFound;
    }

//	public Component findMainWindow() {
//		try {
//			Component f = _finder.find(null, ClassMatcher.of(java.awt.Frame.class));
//			_log.info("Frame: %s. IsWindow: %s",  f.getName(), f instanceof JWindow);
//			return f;
//		} catch (MultipleComponentsFoundException | ComponentNotFoundException e) {
//			e.printStackTrace();
//			return null;
//		}
//	}
}
