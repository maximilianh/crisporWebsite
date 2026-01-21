package ur_rna.Utilities.swing;

import ur_rna.Utilities.annotation.NotNull;

import javax.swing.*;
import java.awt.*;
import java.util.List;

/**
 * Utility class for dealing with JMenuItem
 */
public interface IMenuItem {
    Component getComponent();
    default JMenuItem getMenuItem() {
        Component c;
        return (c = getComponent()) instanceof JMenuItem ? (JMenuItem)c : null;
    }
    /**
     * Gets the container that holds JMenuItem and JMenu children.
     * (if different than getComponent).
     * For example, a JMenu is a component, but its children are contained by a JPopupMenu,
     * not by its own container.
     */
    default JComponent getContainer() { return Menus.getSubItemContainer(getComponent());  }
    default boolean hasSubItems() { return Menus.hasSubItems(getContainer()); }
    @NotNull
    default String getItemName() { return Components.getNameOrText(getComponent());  }
    default List<? extends JMenu> getMenus() { return Menus.getMenus(getContainer()); }
    default List<? extends JMenuItem> getMenuItems() { return Menus.getMenuItems(getContainer()); }
    default JMenuItem findByTreePath(String name) { return Menus.findByTreePath(this, name); }
    default JMenuItem findByName(String name) { return Menus.findByName(this, name); }
    default JMenuItem findByAction(Action find) { return Menus.findByAction(this, find); }

//    static List<IMenuItem> getIMenuItems(Component menu) {
//        final JComponent m = Menus.getSubItemContainer(menu);
//        if (m == null)
//            return Collections.emptyList();
//        int length = m.getComponentCount();
//        ArrayList<IMenuItem> list = new ArrayList<>(length);
//        Component c;
//        for(int i = 0; i < length; i++)
//            if ((c = m.getComponent(i)) instanceof JMenuItem)
//                list.add(from((JMenuItem)c));
//        return list;
//    }

    static IMenuItem from(Component mi) {
        if (mi == null) throw new NullPointerException("Cannot create an IMenuItem from a NULL Component reference.");
        if (mi instanceof IMenuItem)
            return (IMenuItem)mi;
        if (mi instanceof JMenu)
            return new IMenuItem() {
                JMenu c = (JMenu)mi;
                @Override
                public Component getComponent() { return c; }
                public JComponent getContainer() { return c.getPopupMenu(); }
            };
        return new IMenuItem() {
            Component c = mi;
            @Override
            public Component getComponent() { return c; }
        };
    }
//    class MenuItemComponent implements IMenuItem {
//        private JComponent parent;
//        public MenuItemComponent(final JComponent parent) { this.parent=parent; }
//        public void setComponent(final @NotNull JComponent parent) {
//            if (parent == null) throw new NullPointerException("MenuItemComponent cannot have a null component reference.");
//            this.parent=parent;
//        }
//        @Override
//        public Component getComponent() {
//            return parent;
//        }
//    }

//    class MenuItemList<T extends Component> implements IMenuItem {
//        private final Iterable<T> list;
//        public MenuItemList(final Iterable<T> list) { this.list = list; }
//        @Override
//        public Component getComponent() {
//            return null;
//        }
//        @Override
//        public JMenuItem getMenuItem() {
//            return null;
//        }
//        @Override
//        public boolean hasSubItems() {
//            return list.iterator().hasNext();
//        }
//        @Override
//        public String getItemName() {
//            return null;
//        }
//        @Override
//        public List<JMenu> getMenus() {
//            ArrayList<JMenu> menus =new ArrayList<>();
//            for(Component c : list)
//                if (c instanceof JMenu) menus.add((JMenu)c);
//            return menus;
//        }
//        @Override
//        public List<JMenuItem> getMenuItems() {
//            ArrayList<JMenuItem> menus =new ArrayList<>();
//            for(Component c : list)
//                if (c instanceof JMenuItem) menus.add((JMenu)c);
//            return menus;
//        }
//    }
}
