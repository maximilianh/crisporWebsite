package ur_rna.StructureEditor.menus;

import ur_rna.Utilities.swing.MergeItem;
import ur_rna.Utilities.swing.MergeMenu;

import javax.swing.*;

/**
 * A MainMenu is essentially an MergeMenu that is designed to be subclassed and to
 * handle the actions of its child menu items by attaching itself as an
 * ActionListener to all child menu items.
 */
public class MainMenu extends MergeMenu  { //implements ActionListener
    private static final long serialVersionUID = 20160315;
    //private final AppMainFrame mainFrame = AppMainFrame.getAppRoot();

    /**
     * The ActionLister that will be used to delegate ActionEvents
     * produced by the menu items. The default is null, which indicates
     * that menu events are handled internally (i.e. in a sub-class)
     */
//    private final ActionListener menuListener;
    // public MainMenu(String title) {  super(title, MergePosDefault, title); }
    public MainMenu(String title) {  super(title); }
//    public MainMenu(String title, int mergePos) {  this(title, title, mergePos, null); }
//    public MainMenu( String title, ActionListener menuListener) { this(title, title, MergePosDefault, menuListener); }
//    public MainMenu( String title, String mergeName, int mergePos, ActionListener menuListener) {
//        super(title, mergePos, KeyMnemonic.stripMnemonics(mergeName));
//        this.menuListener = menuListener;
//    }

//    @Override  // from ActionListener
//    public void actionPerformed( ActionEvent ev ) {
//        String cmd = ev.getActionCommand();
//        AppActions ac = AppActions.find(cmd);
//        if (ac != null)
//            Program.getInstance().invoke(ac);
//        else if (menuListener != null)
//            menuListener.actionPerformed(ev);
//        else
//            Program.getInstance().invoke(cmd);
//            //AppLog.getDefault().warn("Unhandled Menu Action: \"" + cmd + "\" in " + this.getClass().getName());
//            //onMenuAction(cmd, ev);
//    }

//    @Override
//    /**
//     * {@inheritDoc}
//     *
//     * This override adds the MainMenu as an ActionListener for the item.
//     */
//    protected void initItem(final JMenuItem item, final String text, final String tooltip, final String actionCommand, final KeyStroke keyStroke) {
//        super.initItem(item, text, tooltip, actionCommand, keyStroke);
//        item.addActionListener(this);
//    }

//    /**
//     * Do any actions specific to this menu.
//     *
//     * @param command   The command that signifies a particular action.
//     * @param ev
//     */
//    protected void onMenuAction(String command, final ActionEvent ev) {
//        Program.getInstance().invoke(command);
//        //
//    }

    public MergeItem addItem(final Action a) {
        return addItem(a, MergePosDefault);
    }
    public MergeItem addItem(final Action a, int mergePos) {
        MergeItem mi = new MergeItem(a, mergePos);
        add(mi);
        return mi;
        //return addItem(null, null, (KeyStroke)null, cmd, mergePos);
    }
//    public MergeItem addItem(String text, String tooltip, String keyStroke, Action cmd) {
//        return addItem(text, tooltip, getKey(keyStroke), cmd, MergePosDefault);
//    }
//
//    public MergeItem addItem(String text, String tooltip, String keyStroke, final Action cmd, int mergePos) {
//        return addItem(text, tooltip, getKey(keyStroke), cmd, mergePos);
//    }
//
//    public MergeItem addItem(String text, String tooltip, KeyStroke keyStroke, final Action cmd, int mergePos) {
//        JMenuItem item = new JMenuItem(cmd);
////        if (tooltip == null) tooltip = UiAction.getDesc(cmd);
////        if (keyStroke == null) keyStroke = UiAction.getKeyStroke(cmd);
////        if (text == null) text = UiAction.getText(cmd);
//        //UiAction.getName(cmd)
//        initItem(item, text, tooltip, null, keyStroke);
//        return add(item, mergePos);
//    }
//    public MergeItem addItem(String text, String tooltip, String keyStroke, final ActionListener action, int mergePos) {
//        MergeItem mi = addItem(text, tooltip, keyStroke);
//        replaceAction(mi.uiItem, this, action);
//        mi.setMergePos(mergePos);
//        add(mi);
//        return mi;
//    }
//    /**
//     * Search the registered item's list of ActionListeners for the specified listener ("find"),
//     * and if found, replace it with another listener ("replaceWith")
//     * @param item  The AbstractButton (Button, CheckBox, MenuItem etc) whose list of ActionListeners should be modified.
//     * @param find  The ActionListener to find. If found, it will be removed.
//     * @param replaceWith  The ActionListener to add if the "find" listener was found.
//     * @return Returns true if the list of ActionListeners was modified or false otherwise.
//     */
//    protected boolean replaceAction(final AbstractButton item, final ActionListener find, final ActionListener replaceWith) {
//        for(ActionListener l : item.getActionListeners())
//            if (l == find) {
//                item.removeActionListener(l);
//                item.addActionListener(replaceWith);
//                return true;
//            }
//        return false;
//    }
}
