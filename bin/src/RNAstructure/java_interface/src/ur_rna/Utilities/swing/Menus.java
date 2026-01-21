package ur_rna.Utilities.swing;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

/**
 * A utility class to help with JMenus and JMenuItems
 */
public class Menus {
    public static boolean hasSubItems(Component menu) {
        Container c = getSubItemContainer(menu);
        return c != null && c.getComponentCount() != 0;
    }

    public static JComponent getSubItemContainer(Component menu) {
        if (menu instanceof JMenu)
            return ((JMenu)menu).getPopupMenu();
        if (menu instanceof JComponent)
            return (JComponent)menu;
        return null;
    }

    public static List<JMenu> getMenus(Component menu) {
        JComponent m = getSubItemContainer(menu);
        if (m == null)
            return Collections.emptyList();
        int length = m.getComponentCount();
        ArrayList<JMenu> list = new ArrayList<>(length);
        Component c;
        for(int i = 0; i < length; i++)
            if ((c = m.getComponent(i)) instanceof JMenu)
                list.add((JMenu)c);
        return list;
    }

    public static java.util.List<JMenuItem> getMenuItems(Component menu) {
        final JComponent m = getSubItemContainer(menu);
        if (m == null)
            return Collections.emptyList();
        int max = m.getComponentCount();
        ArrayList<JMenuItem> list = new ArrayList<>(max);
        Component c;
        for(int i = 0; i < max; i++)
            if ((c = m.getComponent(i)) instanceof JMenuItem)
                list.add((JMenuItem)c);
        return list;
    }

    public static JMenuItem findByTreePath(Component menu, String name) {
        return findByTreePath(IMenuItem.from(menu), name);
    }
    public static JMenuItem findByTreePath(IMenuItem parent, String name) {
        if (name == null) return null;
        String[] parts = name.split("->");
        int i = 0;
        while (true) {
            JMenuItem mi = findByName(parent, parts[i]);
            if (++i == parts.length)
                return mi;
            // If we are not at the end, then the found item MUST be a Menu
            //   in order for the next search to proceed.
            if (!(mi instanceof JMenu)) // hasSubItems verifies that it is a JMenu and that it has sub-items
                return null;
            parent = IMenuItem.from(mi);
        }
    }

    public static JMenuItem findByName(Component menu, String name) {
        return findByName(IMenuItem.from(menu), name);
    }

    public static JMenuItem findByName(IMenuItem parent, String name) {
        if (name == null) return null;
        if (name.contains("->"))
            return findByTreePath(parent, name);

        List<? extends JMenuItem> items = parent.getMenuItems();

        // perform search across, then down through hierarchy
        int i = 0; //track index
        for (JMenuItem mi : items) {
            if (name.equals(Components.getNameOrText(mi)))
                return mi;
            if (name.equals("[" + i + "]"))
                return mi;
            i++;
        }
        for (JMenuItem m : items)
            if (m instanceof JMenu) {
                JMenuItem mi = findByName(IMenuItem.from(m), name);
                if (mi != null)
                    return mi;
            }
        return null;
    }

    public static void setTextWithMnemonic(final AbstractButton item,  final String text, final Consumer<String> setTextDirectly) {
        int mn = KeyMnemonic.getMnemonic(text);
        if (mn == 0)
            setTextDirectly.accept(text);
        else {
            item.setMnemonic(mn);
            int mni = KeyMnemonic.getMnemonicIndex(text);
            setTextDirectly.accept(KeyMnemonic.stripMnemonics(text));
            if (mni != item.getDisplayedMnemonicIndex())
                item.setDisplayedMnemonicIndex(mni);
        }
    }

    public static JMenuItem findByAction(Component menu, Action action) {
        return findByAction(IMenuItem.from(menu), action);
    }
    public static JMenuItem findByAction(IMenuItem menu, Action action) {
        if (action == null) return null;
        // perform search across, then down through hierarchy
        for (JMenuItem mi : menu.getMenuItems())
            if (action.equals(mi.getAction()))
                return mi;
        for (JMenu m : menu.getMenus()) {
            JMenuItem mi = findByAction(IMenuItem.from(m), action);
            if (mi != null)
                return mi;
        }
        return null;
    }
    public static void setEnabled(IMenuItem parent, boolean enable) {
        for(JMenu m : parent.getMenus())
            setEnabled(IMenuItem.from(m), enable);
        for(JMenuItem mi : parent.getMenuItems())
            mi.setEnabled(enable);
    }

//
//    public static void mergeMenus(JComponent parent, List<? extends JMenu> menus) {
//        parent = Menus.getSubItemContainer(parent);
//
//        if (menus.size() == 0)
//            return Collections.emptyList();
//
//        List<JMenu> list = new ArrayList<>(menus.size());
//
//        Iterator<? extends JMenu> it = menus.iterator();
//
//        list.add(it.next()); // must have at least one element.
//        while(it.hasNext()) {
//            JMenu jm = it.next();
//            String name = Menus.getItemName(jm);
//            if (!isEmpty(name)) {
//                for (JMenu j1 : list)
//
//            }
//
//        }
//
//        for (int i = 2; i < lists.length; i++)
//            list = merge(list, lists[i].getSubItems());
//        return new MenuList(list);
//    }
//
//    private static void mergeSubMenus(final JMenu m) {
//        List<JMenu> list = Menus.getMenus(m);
//        int length = list.size();
//
//        loop:
//        boolean match;
//
//    }
//
//    public static List<Component> mergeSubItems(List<? extends Component> list1, List<? extends Component> list2) {
//        if (list1 == null) list1 = EMPTY_COMPONENT_LIST;
//        if (list2 == null) list2 = EMPTY_COMPONENT_LIST;
//        final Component[] items = new Component[list1.size()+list2.size()];
//        ObjTools.copyTo(list1, items);
//        ObjTools.copyTo(list2, items, list1.size());
//        Arrays.sort(items, IMergeItem::compareByMergePos);
//        int count = 0;
//        ArrayList<String> names = new ArrayList<>(items.length);
//        for (int i = 0; i < items.length; i++) {
//            IMergeItem menu =
//            if ()
//                if (items[i] instanceof MergeMenu) {
//                    MergeMenu m = (MergeMenu) items[i];
//                    String name = m.getName();
//                    int found = name.length() == 0 ? -1 : names.indexOf(name); //no need to search if ="", because we don't merge items without a non-empty mergeName
//                    if (found == -1) {
//                        items[count++] = m;
//                        names.add(name);
//                    } else {
//                        //merge the two items
//                        MergeMenu m1 = (MergeMenu) items[found];
//                        if (!m1.isMergeCopy) {
//                            items[found] = m1 = (MergeMenu) m1.cloneShallow(); //replace the original element with this clone so we can change properties and add new items.
//                            m1.isMergeCopy = true;
//                        }
//                        m1.items = merge(m1.items, m.items);
//                    }
//                } else
//                    count++;
//        }
//        return Arrays.asList(items).subList(0, count);
//    }

}
