package ur_rna.RNAstructureUI.menus;

import ur_rna.RNAstructureUI.AppMainFrame;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.swing.*;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import static ur_rna.Utilities.Strings.isEmpty;

/**
 * A MainMenu is essentially an MergeMenu that is designed to be subclassed and to
 * handle the actions of its child menu items by attaching itself as an
 * ActionListener to all child menu items.
 */
public class MainMenu extends MergeMenu implements ActionListener {
    private static final long serialVersionUID = 20160315;
    private final AppMainFrame mainFrame = AppMainFrame.getAppRoot();

    /**
     * The ActionLister that will be used to delegate ActionEvents
     * produced by the menu items. The default is null, which indicates
     * that menu events are handled internally (i.e. in a sub-class)
     */
    private final ActionListener menuListener;

    public MainMenu(String title) {  this(title, MergePosInitial, null); }
    public MainMenu(String title, int mergePos) {  this(title, mergePos, null); }
    public MainMenu( String title, ActionListener menuListener) { this(title, MergePosInitial, menuListener); }
    public MainMenu( String title, int mergePos, ActionListener menuListener) {
        super(title, mergePos);
        this.menuListener = menuListener;
    }

    @Override  // from ActionListener
    public void actionPerformed( ActionEvent ev ) {
        if (menuListener == null)
            onMenuAction(ev.getActionCommand(), ev);
        else
            menuListener.actionPerformed(ev);
    }

    /**
     * Do any actions specific to this menu.
     *
     * @param command   The command that signifies a particular action.
     * @param ev
     */
    protected void onMenuAction(String command, final ActionEvent ev) {
        AppLog.getDefault().warn("Unhandled Menu Action: \"" + command + "\" in " + this.getClass().getName());
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

	/**
	 * Add a check box menu item with tooltip text.
	 *
	 * @param text    The menu item text.
	 * @param tooltip The tooltip text.
	 */
	public CheckMergeItem addCheckItem(String text, String tooltip) {
		return addCheckItem(text, tooltip, false);
	}
	public CheckMergeItem addCheckItem(String text, String tooltip, boolean checked) {
        CheckMergeItem check =createCheckItem(text, tooltip);
		check.setSelected(checked);
		add(check);
		return check;
	}
	protected CheckMergeItem createCheckItem(String text, String tooltip) {
        CheckMergeItem item = new CheckMergeItem();
		initItem(item, text, tooltip, null, null);
		return item;
	}

	/**
	 * Add a plain text menu item with tooltip text.
	 *
	 * @param text    The menu item text.
	 * @param tooltip The tooltip text.
	 */
	public MergeItem addItem(String text, String tooltip) {
		return addItem(text, tooltip, (KeyStroke) null, null);
	}

	/**
	 * Add a plain text menu item with tooltip text.
	 *
	 * @param text    The menu item text.
	 * @param tooltip The tooltip text.
	 */
	public MergeItem addItem(String text, String tooltip, String shortcutKey) {
		return addItem(text, tooltip, AcceleratorKey.getKey(shortcutKey), null);
	}

	/**
	 * Add a menu item with tooltip text and a key stroke which
	 * can activate it.
	 *
	 * @param text        The menu item text.
	 * @param tooltip     The tooltip text.
	 * @param shortcutKey The key stroke.
	 */
	public MergeItem addItem(String text, String tooltip, char shortcutKey) {
		return addItem(text, tooltip, AcceleratorKey.getKey(shortcutKey), text);
	}

	/**
	 * Add a plain text menu item with tooltip text and a key stroke which
	 * can activate it.
	 *
	 * @param text    The menu item text.
	 * @param tooltip The tooltip text.
	 * @param key     The key stroke.
	 */
	public MergeItem addItem(String text, String tooltip, KeyStroke key) {
		return addItem(text, tooltip, key, null);
	}

	/**
	 * Add a plain text menu item with tooltip text.
	 *
	 * @param text    The menu item text.
	 * @param tooltip The tooltip text.
	 */
	public MergeItem addItem(String text, String tooltip, String shortcutKey, String actionCommand) {
		return addItem(text, tooltip, AcceleratorKey.getKey(shortcutKey), actionCommand);
	}

	/**
	 * Add a plain text menu item with tooltip text.
	 *
	 * @param text    The menu item text.
	 * @param tooltip The tooltip text.
	 */
	public MergeItem addItem(String text, String tooltip, char shortcutKey, String actionCommand) {
		return addItem(text, tooltip, AcceleratorKey.getKey(shortcutKey), actionCommand);
	}

	/**
	 * Add a menu item with tooltip text.
	 *
	 * @param text          The menu item text.
	 * @param tooltip       The tooltip text.
	 * @param actionCommand A string specifying the action identifier
	 */
	public MergeItem addItem(String text, String tooltip, KeyStroke keyStroke, String actionCommand) {
        MergeItem item = createMenuItem(text, tooltip, actionCommand, keyStroke);
		add(item);
        return item;
	}

	protected MergeItem createMenuItem(String text, String tooltip, String actionCommand, KeyStroke keyStroke) {
        MergeItem item = new MergeItem();
		initItem(item, text, tooltip, actionCommand, keyStroke);
		return item;
	}

	/**
	 * Perform initial setup of the JMenuItem
	 * @param item 			The JMenuItem to initialize.
	 * @param text 			The text for the item, which can contain an ampersand (&amp;)
	 *                		in front of the mnemonic character.
	 * @param tooltip 		Descriptive text for the item.
	 * @param actionCommand The command this will invoke when clicked.
     * @param keyStroke 	An accelerator shortcut key.
     */
	protected void initItem(JMenuItem item, String text, String tooltip, String actionCommand, KeyStroke keyStroke) {
		if (text != null) {
			int mn = KeyMnemonic.getMnemonic(text);
			if (mn != 0)
				text = KeyMnemonic.stripMnemonics(text);
			if (isEmpty(actionCommand)) actionCommand = text;
			actionCommand = KeyMnemonic.stripMnemonics(actionCommand);
			item.setText(text);
			if (mn != 0)
				item.setMnemonic(mn);
		}
		if (actionCommand != null) {
			item.setName(actionCommand);
			item.setActionCommand(actionCommand);
		}
		if (tooltip != null)
			item.setToolTipText(tooltip);
		if (keyStroke != null)
			item.setAccelerator(keyStroke);
        item.addActionListener(this);
        if (item instanceof IMergeItem)
			((IMergeItem)item).setMergePos(subItemMergePos);
	}

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
