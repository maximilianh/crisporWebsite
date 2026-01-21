package ur_rna.Utilities.swing;

import ur_rna.Utilities.ObjTools;

import javax.swing.*;
import java.awt.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static ur_rna.Utilities.Strings.isEmpty;

/**
 *
 */
public class MergeManager  {
    public static final Component[] EMPTY_COMPONENT_ARRAY = new Component[0];
    private Map<Component, Integer> mergePosMap = new HashMap<>();
    private Map<Container, Component[]> originalComponents = new HashMap<>();
    private static final Integer DEFAULT_MERGE_POS = IMergeItem.MergePosDefault;

//    public <T> T setMergePos(T c, int mergePos) {
//        if (c == null)
//            throw new IllegalArgumentException("The component cannot be null.");
//        if (c instanceof IMergeItem) {
//            ((IMergeItem) c).setMergePos(mergePos);
//            return c;
//        }
////        if (c instanceof Component) {
////            mergePosMap.put((Component)c, mergePos);
////        }
//        throw new IllegalArgumentException("The component's merge position cannot be set.");
//        //return c;
//    }
//
//    public int getMergePos(Object c) {
//        if (c == null)
//            throw new IllegalArgumentException("The component cannot be null.");
//        if (c instanceof IMergeItem)
//            return ((IMergeItem) c).getMergePos();
////        if (c instanceof Component) {
////            return mergePosMap.getOrDefault(c, DEFAULT_MERGE_POS);
////        }
//        return IMergeItem.MergePosDefault;
//    }

    public void reset() {
        for (Map.Entry<Container, Component[]> entry : originalComponents.entrySet()) {
            setComponents(entry.getKey(), entry.getValue());
        }
        originalComponents.clear();
    }

    private String getMergeName(Object c) {
        if (c instanceof IMergeItem)
            return ((IMergeItem)c).getMergeName();
        if (c instanceof AbstractButton)
            return IMergeItem.getMergeName((Component) c);
        if (c instanceof IMenuItem)
            return ((IMenuItem)c).getItemName();
        return null;
    }

    public boolean merge(final Component p1, final Component p2) {
        Container c1 = Menus.getSubItemContainer(p1), c2 = Menus.getSubItemContainer(p2);
        if (c1 == null) throw new IllegalArgumentException("Cannot merge items into a Component that is not a Container and has no associated Container.");
        if (c1 == null || c2.getComponentCount() == 0) return false;
        Component[] original1 = c1.getComponents();
        Component[] original2 = c2.getComponents();
        Component[] result = new Component[original1.length + original2.length];
        boolean modified = mergeItems(original1, original2, result);
        if (modified) {
            storeComponents(c1, original1);
            storeComponents(c2, original2);
            setComponents(c1, result);
        }
        return modified;
    }

    public boolean merge(final Component parent, final Component[] addSubItems) {
        Container container = Menus.getSubItemContainer(parent);
        if (container == null) throw new IllegalArgumentException("Cannot merge items into a Component that is not a Container and has no associated Container.");
        if (addSubItems.length == 0) return false;
        Component[] original = container.getComponents();
        Component[] result = new Component[original.length + addSubItems.length];
        boolean modified = mergeItems(original, addSubItems, result);
        if (modified) {
            storeComponents(container, original);
            setComponents(container, result);
        }
        return modified;
    }

    private boolean mergeItems(final Component[] existing, final Component[] adding, final Component[] merged) {
        boolean modified = false;
        final String[] names = new String[existing.length];
        for (int i = 0; i < existing.length; i++)
            names[i] = getMergeName(existing[i]);
        System.arraycopy(existing, 0, merged, 0, existing.length);
        int pos = existing.length;
        // Container cAdd, cExist;
        for (Component add : adding) {
            //if ((cAdd = Menus.getSubItemContainer(add)) != null) {
            String name = getMergeName(add);
            int found = isEmpty(name) ? -1 : ObjTools.indexOf(names, name, true);
            // cExist = found == -1 ? null : Menus.getSubItemContainer(existing[found]);
            if (found != -1 && Menus.getSubItemContainer(existing[found]) != null) {
                merge(existing[found], add);
                //continue;
            } else {
                modified = true;
                merged[pos++] = add;
            }
            //}
        }
        //Arrays.sort(merged, 0, pos, IMergeItem::compareByMergePos);
        return modified;
    }

    public void sort(final Component parent) {
        int max = Integer.MIN_VALUE;
        Container container = Menus.getSubItemContainer(parent);
        if (container == null) throw new IllegalArgumentException("Cannot sort items in a Component that is not a Container and has no associated Container.");
        int length = container.getComponentCount();
        boolean sorted = true;
        for (int i = 0; i < length; i++) {
            Component c = container.getComponent(i);
            int pos = IMergeItem.getValidMergePos(c);
            if (pos < max) {
                sorted = false;
                break;
            }
            max = pos;
            if (Menus.hasSubItems(c))
                sort(c);
        }
        if (!sorted) {
            Component[] list = container.getComponents();
            Arrays.sort(list, IMergeItem::compareByMergePos);
            setComponents(container, list);
        }
    }

    private void storeComponents(final Container parent, final Component[] original) {
        if (!originalComponents.containsKey(parent)) // do not change it if already set.
            originalComponents.put(parent, original);
    }

    private void setComponents(final Container parent, final Component[] list) {
        parent.removeAll();
        for(Component c : list) {
            if (c != null)
                parent.add(c);
        }
    }
    public void clearCache(final Container component) {
        originalComponents.remove(component);
    }

//    private Component[] getComponents(final Component c) {
//        if (c instanceof IContainer)
//            return ((IContainer) c).getComponents();
//        if (c instanceof Container)
//            return ((Container) c).getComponents();
//        return EMPTY_COMPONENT_ARRAY;
//    }
//
//    private boolean arraysEqual(final Component[] original, final Component[] merged) {
//        if (merged.length > original.length && merged[original.length] != null)
//            return false;
//        for (int i = 0; i < original.length; i++) {
//            if (!Objects.equals(original[i], merged[i]))
//                return false;
//        }
//        return true;
//    }
}
