/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.Utilities.swing;

import javax.swing.*;
import java.awt.*;

/**
 * A class that creates a menu.
 *
 * @author Richard M. Watson
 */
public class MergeMenu extends JMenu implements IMergeItem, IMenuItem {
	public int subItemMergePos = 0;
	//private String mergeName;
	private int mergePos = MergePosInitial;

	public MergeMenu(String text) { this(text, MergePosInitial); }
    public MergeMenu(Action a, int mergePos) { super(a); setMergePos(mergePos); }
    public MergeMenu(Action a, String text, int mergePos) { this(a, mergePos); setText(text); }
	public MergeMenu(String text, String name) { this(text, name, MergePosInitial);	}
	public MergeMenu(String text, String name, int mergePos) { this(text, mergePos); setName(name); }
    public MergeMenu(String text, int mergePos) { setText(text); setMergePos(mergePos); }

    public void addAll(Component... components) { for(Component c : components) add(c); }

	@Override
    public void setText(String text) { Menus.setTextWithMnemonic(this, text, super::setText); }

	/**
	 * MergeMenu Conversion Constructor.
	 *
	 * @param oldMenu The menu to convert.
	 */
	public MergeMenu(JMenu oldMenu) {
		this(oldMenu.getText());
		for (Component component : oldMenu.getComponents())
				add(component);
	}

	@Override
	public int getMergePos() { return mergePos;}
	@Override
	public MergeMenu setMergePos(final int value) {mergePos = value; return this; }

//	@NotNull
//	@Override
//	public String getName() {
//		return isEmpty(mergeName) ?
//				super.getName() :
//				mergeName;
//	}

//	public String getMergeName() { return this.mergeName; }
//	public void setMergeName(final String name) { this.mergeName = name; }

	/**
	 * Enable or Disable the menus in this list.
	 */
	public void setSubItemsEnabled(final boolean enabled) {
		for (JMenuItem mi : this.getMenuItems()) {
			if (mi != null)
				mi.setEnabled(enabled);
		}
	}

	/**
	 * Set visibility for all sub-items
	 */
	public void setSubItemsVisible(final boolean visible) {
		for (JMenuItem mi : this.getMenuItems()) {
			if (mi != null)
				mi.setVisible(visible);
		}
	}

//	public List<Component> merged() {
//		return merge(EMPTY_COMPONENT_LIST, getMenuItems());
//	}

//	public List<IMergeItem> getMergeItems() {return IMergeItem.asMergeItems(getMenuItems()); }

//	public void add(JComponent menuItem) {
//		items.add(new MergeInfo(menuItem, subItemMergePos));
//	}
//	public void insert(int index, JComponent menuItem) {
//		items.add(index, new MergeInfo(menuItem, subItemMergePos));
//	}
//	public void add(JComponent... menuItems) {
//		for (JComponent m : menuItems) add(new MergeInfo(m, subItemMergePos));
//	}
//	public void add(MergeInfo m) {
//		items.add(m);
//	}
//	public void insert(int index, MergeInfo m) {
//		items.add(index, m);
//	}
//	public void add(MergeInfo... menus) {
//		for (MergeInfo m : menus) add(m);
//	}

	/**
	 * Creates a new menu item attached to the specified
	 * <code>Action</code> object and appends it to the end of this menu.
	 *
	 * @param a the <code>Action</code> for the menu item to be added
	 * @see Action
	 */
	@Override
	public MergeItem add(final Action a) { return (MergeItem)add(new MergeItem(a, subItemMergePos)); }
	public CheckMergeItem addCheck(final Action a) { return (CheckMergeItem)add(new CheckMergeItem(a, subItemMergePos)); }
	public MergeItem add(final Action a, int mergePos) { return add(new MergeItem(a), mergePos); }
	public MergeItem add(final MergeItem mi, int mergePos) {
		super.add(mi);
		IMergeItem.setMergePos(mi, mergePos);
		return mi;
	}

	@Override
	public JMenuItem add(final JMenuItem mi) {
		if (mi instanceof IMergeItem)
			IMergeItem.setInitialMergePos((IMergeItem) mi, subItemMergePos);
		return super.add(mi);
	}
	private int checkMergePos(int mergePos) {
		return mergePos; //mergePos == MergePosDefault ? subItemMergePos :
	}

	/**
	 * Add a menu separator
	 */
	public void addSeparator() {  addSeparator(subItemMergePos); }

	/**
	 * Add a menu separator with the specified mergePos
	 */
	public void addSeparator(int mergePos) { add(new MergeItem.Separator(checkMergePos(mergePos))); }

	/**
	 * Increments by 10 the current mergePos at which sub-items will be added by default.
	 */
	public void insertMergeGap() {
		subItemMergePos += 10;
	}
	/**
	 * Increments (by gap) the current mergePos at which sub-items will be added by default.
	 */
	public void insertMergeGap(int gap) {
		subItemMergePos += gap;
	}
	/**
	 * Gets the current mergePos at which sub-items will be added by default.
	 * @return An integer representing the current merge position
	 *         at which sub-items will be added by default.
	 */
	public int getSubItemMergePos() {
		return subItemMergePos;
	}
	/**
	 * Sets the current mergePos at which sub-items will be added by default.
	 */
	public void setSubItemMergePos(final int subItemMergePos) {
		this.subItemMergePos = subItemMergePos;
	}

//	public boolean areMergeItemsSorted() {
//		int max = Integer.MIN_VALUE;
//		Container container = getPopupMenu();
//		int length = container.getComponentCount();
//		boolean sorted = true;
//		for (int i = 0; i < length; i++) {
//			Component c = container.getComponent(i);
//			int pos = IMergeItem.getMergePos(c);
//			if (pos < max) return false;
//			if (Menus.hasSubItems(c) && !)
//
//		}
//
//	}
//
//	public void sortMergeItems() {
//		if (areMergeItemsSorted()) return;
//
//		if (!sorted) {
//			Component[] list = container.getComponents();
//			Arrays.sort(list, IMergeItem::compareByMergePos);
//			setComponents(container, list);
//		}
//	}

	//	public MergeInfo add(JMenuItem m, int mergePos) {
//		MergeInfo mi = new MergeInfo(m, mergePos);
//		add(mi);
//		return mi;
//	}
//	public MergeInfo insert(int index, JMenuItem m, int mergePos) {
//		MergeInfo mi = new MergeInfo(m, mergePos);
//		insert(index, mi);
//		return mi;
//	}

//
//	public java.util.List<MergeInfo> getSubItems() {
//		return items;
//	}
//	public MergeMenu[] getSubMenus() {
//		int len = 0;
//		for (MergeInfo mi : items) {
//			if (mi.item instanceof MergeMenu)
//				len++;
//		}
//		MergeMenu[] arr = new MergeMenu[len];
//		len = 0;
//		for (MergeInfo mi : items) {
//			if (mi.item instanceof MergeMenu)
//				arr[len++] = (MergeMenu)mi.item;
//		}
//		return arr;
//	}
//

//    public List<MergeItem> merge(MergeItem menu) {
//        List<MergeItem> list = new ArrayList<>();
//        mergeInto(list);
//        return list;
//    }
//
//    public void mergeInto(List<MergeItem> mergeInto) {
//        if (!this.hasChildren())
//            return;
//        mergeInto(this.getSubItems(), )
//    }

//
//
//	public Object clone()  {
//		MergeMenu m = (MergeMenu)super.clone();
//		m.items = new ArrayList<>(m.items);
//		return m;
//	}
//	/**
//	 * Creates a shallow copy of this MergeMenu, unlike (@link #clone},
//	 * which creates a deeper copy in which the sub-item array is also cloned.
//	 */
//	protected final Object cloneShallow()  {
//			return super.clone();
//	}
//
//	/**
//	 * Add a check box menu item with tooltip text.
//	 *
//	 * @param text    The menu item text.
//	 * @param tooltip The tooltip text.
//	 */
//	public MergeItem addCheckItem(String text, String tooltip) {
//		return addCheckItem(text, tooltip, false);
//	}
//	public MergeItem addCheckItem(String text, String tooltip, boolean checked) {
//		MergeItem check = add(createCheckItem(text, tooltip));
//		check.setSelected(checked);
//		return check;
//	}
//	protected JCheckBoxMenuItem createCheckItem(String text, String tooltip) {
//		JCheckBoxMenuItem item = new JCheckBoxMenuItem();
//		initItem(item, text, tooltip, null, null);
//		return item;
//	}
//
//	/**
//	 * Add a plain text menu item with tooltip text.
//	 *
//	 * @param text    The menu item text.
//	 * @param tooltip The tooltip text.
//	 */
//	public MergeItem addItem(String text, String tooltip) {
//		return addItem(text, tooltip, (KeyStroke) null, null);
//	}
//
//	/**
//	 * Add a plain text menu item with tooltip text.
//	 *
//	 * @param text    The menu item text.
//	 * @param tooltip The tooltip text.
//	 */
//	public MergeItem addItem(String text, String tooltip, String shortcutKey) {
//		return addItem(text, tooltip, AcceleratorKey.getKey(shortcutKey), null);
//	}
//
//	/**
//	 * Add a menu item with tooltip text and a key stroke which
//	 * can activate it.
//	 *
//	 * @param text        The menu item text.
//	 * @param tooltip     The tooltip text.
//	 * @param shortcutKey The key stroke.
//	 */
//	public MergeItem addItem(String text, String tooltip, char shortcutKey) {
//		return addItem(text, tooltip, AcceleratorKey.getKey(shortcutKey), text);
//	}
//
//	/**
//	 * Add a plain text menu item with tooltip text and a key stroke which
//	 * can activate it.
//	 *
//	 * @param text    The menu item text.
//	 * @param tooltip The tooltip text.
//	 * @param key     The key stroke.
//	 */
//	public MergeItem addItem(String text, String tooltip, KeyStroke key) {
//		return addItem(text, tooltip, key, null);
//	}
//
//	/**
//	 * Add a plain text menu item with tooltip text.
//	 *
//	 * @param text    The menu item text.
//	 * @param tooltip The tooltip text.
//	 */
//	public MergeItem addItem(String text, String tooltip, String shortcutKey, String actionCommand) {
//		return addItem(text, tooltip, AcceleratorKey.getKey(shortcutKey), actionCommand);
//	}
//
//	/**
//	 * Add a plain text menu item with tooltip text.
//	 *
//	 * @param text    The menu item text.
//	 * @param tooltip The tooltip text.
//	 */
//	public MergeItem addItem(String text, String tooltip, char shortcutKey, String actionCommand) {
//		return addItem(text, tooltip, AcceleratorKey.getKey(shortcutKey), actionCommand);
//	}
//
//	/**
//	 * Add a menu item with tooltip text.
//	 *
//	 * @param text          The menu item text.
//	 * @param tooltip       The tooltip text.
//	 * @param actionCommand A string specifying the action identifier
//	 */
//	public MergeItem addItem(String text, String tooltip, KeyStroke keyStroke, String actionCommand) {
//		JMenuItem item = createMenuItem(text, tooltip, actionCommand, keyStroke);
//		return add(item);
//	}
//
//	protected JMenuItem createMenuItem(String text, String tooltip, String actionCommand, KeyStroke keyStroke) {
//		JMenuItem item = new JMenuItem();
//		initItem(item, text, tooltip, actionCommand, keyStroke);
//		return item;
//	}
//
//	/**
//	 * Perform initial setup of the JMenuItem
//	 * @param item 			The JMenuItem to initialize.
//	 * @param text 			The text for the item, which can contain an ampersand (&amp;)
//	 *                		in front of the mnemonic character.
//	 * @param tooltip 		Descriptive text for the item.
//	 * @param actionCommand The command this will invoke when clicked.
//     * @param keyStroke 	An accelerator shortcut key.
//     */
//	protected void initItem(JMenuItem item, String text, String tooltip, String actionCommand, KeyStroke keyStroke) {
//		if (text != null) {
//			int mn = KeyMnemonic.getMnemonic(text);
//			if (mn != 0)
//				text = KeyMnemonic.stripMnemonics(text);
//			if (isEmpty(actionCommand)) actionCommand = text;
//			actionCommand = KeyMnemonic.stripMnemonics(actionCommand);
//			item.setText(text);
//			if (mn != 0)
//				item.setMnemonic(mn);
//		}
//		if (actionCommand != null) {
//			item.setName(actionCommand);
//			item.setActionCommand(actionCommand);
//		}
//		if (tooltip != null)
//			item.setToolTipText(tooltip);
//		if (keyStroke != null)
//			item.setAccelerator(keyStroke);
//	}
//
//	@Override
//	public String toString() {
//		String type = (getMenuItem() == SEPARATOR) ? "(SEPARATOR)" : getMenuItem().getClass().getSimpleName();
//		return String.format("%s [ Name: \"%s\"  Text: \"%s\" Pos: %d ]", type, this.getName(), this.getMenuItem().getText(), this.getMergePos());
//	}

//	public JMenu toJMenu() {
//		JMenuItem jSrc = this.uiItem;
//		JMenu jDst = new JMenu();
//		CloneJMenuItems(jSrc, jDst);
//		for (MergeItem mi : getSubItems())
//			jDst.add(mi.uiItem);
//		return jDst;
//	}
//	public static void CloneJMenuItems(JMenuItem jSrc, JMenuItem jDst) {
//		jDst.setText(jSrc.getText());
//		jDst.setToolTipText(jSrc.getToolTipText());
//		jDst.setEnabled(jSrc.isEnabled());
//		jDst.setName(jSrc.getName());
//		jDst.setMnemonic(jSrc.getMnemonic());
//		// Safely clone Accelerator KeyStroke (JMenu does not throws an Error when calling setAccelerator, but no Exception is raised for getAccelerator
//		if (!(jDst instanceof JMenu))
//			jDst.setAccelerator(jSrc.getAccelerator());
//		jDst.setAction(jSrc.getAction());
//		jDst.setActionCommand(jSrc.getActionCommand());
//		jDst.setActionMap(jSrc.getActionMap());
//		jDst.setSelected(jSrc.isSelected());
//		jDst.setForeground(jSrc.getForeground());
//		for (ActionListener a : jSrc.getActionListeners())
//			jDst.addActionListener(a);
//	}
}