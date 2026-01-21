package ur_rna.StructureEditor;

import ur_rna.StructureEditor.menus.FileMenu;
import ur_rna.StructureEditor.menus.SettingsMenu;
import ur_rna.StructureEditor.menus.WindowsMenu;
import ur_rna.StructureEditor.services.RecentFileList;
import ur_rna.StructureEditor.services.RecentFileMenuManager;
import ur_rna.StructureEditor.ui.SliderToolWindow;
import ur_rna.StructureEditor.windows.ChildFrame;
import ur_rna.StructureEditor.windows.DashboardFrame;
import ur_rna.StructureEditor.windows.DrawWindow;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.OSInfo;
import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.Strings;
import ur_rna.Utilities.swing.*;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.BevelBorder;
import javax.swing.event.AncestorEvent;
import javax.swing.event.AncestorListener;
import javax.swing.event.InternalFrameEvent;
import java.awt.*;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.prefs.Preferences;

import static ur_rna.Utilities.Strings.asBool;
import static ur_rna.Utilities.Strings.fmt;
import static ur_rna.Utilities.swing.AcceleratorKey.resetTabTraversalKeys;

/**
 * The application's main window, which hosts the menubar, toolbar, statusbar, and the
 * contains the desktop on which child windows are shown.
 */
public class MainFrame extends MdiParentFrame {
    private AppLog log = AppLog.getDefault();
    private final JMenuBar menuBar;
    private final JToolBar toolBar;
    private final JPanel statusBar;
    private final JLabel infoLabel;
    public final JDesktopPane desktop;
    private final Program program;
    // private MenuList mainMenus;
    private RecentFileList recentFiles;
    private RecentFileMenuManager recentFilesMenu;
    private DashboardFrame dashboard;
    private MergeManager menuMerger = new MergeManager(); // used to sort and merge the main menus and toolbar buttons with those from child windows.

    /**
     * Constructor.
     * Create the frame, set the icon and layout, add the toolbar, back panel,
     * and info label, add a focus listener that shows the default menus, then
     * set the size and location before setting the frame to be visible.
     */
    public MainFrame() {
        super(Program.TITLE);
        this.program = Program.getInstance();
//        String iconString = "Icon.gif";
//        setIconImage(ImageGrabber.tryGetImage(iconString));
        setLayout(new BorderLayout());
        desktop = buildDesktop();

        setJMenuBar(menuBar = new JMenuBar());
        createMainMenus();

        toolBar = buildToolBar();
        createMainToolButtons();

        program.settings().addChangeHandler(this::settingsChanged);
        recentFiles = new RecentFileList();
        dashboard = new DashboardFrame();

        add(toolBar, BorderLayout.PAGE_START);
        add(desktop, BorderLayout.CENTER);
        infoLabel = buildInfoLabel();
        statusBar = buildStatusBar(infoLabel);
        //statusBar.add(infoLabel, BorderLayout.SOUTH);
        resetInfoLabel();

        updateActiveFrameInfo();
        setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        setSize(1024, 730);
        setLocationRelativeTo(null);
        setVisible(true);

        loadRecentFiles();
        addMainActions();
        resetTabTraversalKeys(desktop);
    }
    private void settingsChanged() {
        Settings s = program.settings();
        if (!s.ShowDashboard && dashboard.isVisible())
            hideDashboard();
        else if (s.ShowDashboard && getChildFrames().size()==0)
            showDashboard(false);
        for(ChildFrame f : getChildFrames())
            f.programSettingsUpdated();
    }
    private void addMainActions() {
        addAction("show-rotate", e->new SliderToolWindow().show());
    }

    private Preferences getRecentFilesNode() { return program.prefs().user().node("RecentFiles"); }
    private void loadRecentFiles() {
        recentFiles.load(getRecentFilesNode());
        updateRecentFilesUI();
    }
    public void addRecentFile(String path, FileType type) {
        recentFiles.add(path, type);
        recentFiles.save(getRecentFilesNode());
        updateRecentFilesUI();
    }
    private void updateRecentFilesUI() {
        dashboard.setRecent(recentFiles);
        recentFilesMenu.setRecent(recentFiles);
        menuMerger.clearCache(recentFilesMenu.getComponent()); // clear so any previously-cached list is erased.
        menuMerger.sort(recentFilesMenu.getComponent());
    }

    public void cleanRecent(final boolean removeAll) {
        if (removeAll)
            recentFiles.clear();
        else {
            for (RecentFileList.RecentFile f : recentFiles.getFiles().toArray(new RecentFileList.RecentFile[recentFiles.size()]))
                if (!(new File(f.path).exists()))
                    recentFiles.remove(f);
        }
        recentFiles.save(getRecentFilesNode());
        updateRecentFilesUI();
    }


    /**
     * Create the gray back panel that fills the main frame.
     *
     * @return The back panel.
     */
    private JDesktopPane buildDesktop() {
        JDesktopPane pane = new JDesktopPane();
        pane.setBackground(Color.GRAY);
        pane.removeAll();
        return pane;
    }

    private JPanel buildStatusBar(JLabel lbl) {
        JPanel p = new JPanel();
        p.setLayout(new BoxLayout(p, BoxLayout.X_AXIS));
        p.setBorder(new BevelBorder(BevelBorder.LOWERED));
        p.add(lbl, BorderLayout.SOUTH);

        //Dimension d = lbl.getPreferredSize();
        //Insets ins = p.getBorder().getBorderInsets(p);
        //d.setSize(d.getWidth() + ins.left + ins.right, d.getHeight() + ins.top + ins.bottom);
        //p.setPreferredSize(d);

        this.add(p, BorderLayout.SOUTH);
        return p;
    }

    public static final int FileMenuSection = -100;
    public static final int EditMenuSection = -80;
    public static final int ViewMenuSection = -40;
    public static final int FormatMenuSection = -20;
    public static final int SettingsMenuSection = 900;
    public static final int WindowsMenuSection = 910;
    public static final int HelpMenuSection = 1000;

    private void createMainMenus() {
        FileMenu file = new FileMenu();
        SettingsMenu settings = new SettingsMenu();
        WindowsMenu windows = new WindowsMenu(desktop, this::windowMenuFilter);
        MergeMenu edit = new MergeMenu("&Edit"),
                  help = new MergeMenu("&Help");
        addMainMenu(file, FileMenuSection);
        addMainMenu(edit, EditMenuSection); //Add empty Edit menu as placeholder.
        addMainMenu(settings, SettingsMenuSection);
        addMainMenu(windows, WindowsMenuSection);
        addMainMenu(help, HelpMenuSection); //Add help menu
        if (asBool(System.getProperty("debug")))
            file.add(new UiAction("Test Icons", e->{ if (testIconsTimer.isRunning()) testIconsTimer.stop(); else testIconsTimer.start(); }).setKeyStroke("*I"));

        help.add(AppActions.SHOW_LOCAL_HELP);
        help.add(AppActions.SHOW_ONLINE_HELP);
        help.add(AppActions.OPEN_WEBSITE);
        help.addSeparator();
        help.add(AppActions.SEND_FEEDBACK);
        help.addSeparator();
        help.add(AppActions.SHOW_ABOUT_WINDOW);
        recentFilesMenu = new RecentFileMenuManager(file.recentFilesMenu());
    }
    /** Determines which windows should be listed in the "Windows" menu */
    private boolean windowMenuFilter(final JInternalFrame frame) {
        return frame != dashboard;
    }
    private void addMainMenu(final MergeMenu menu, final int mergePos) {
        menuBar.add(menu);
        menu.setMergePos(mergePos);
    }

    /**
     * Create the label at the bottom of the main frame.
     * <br><br>
     * This informative label holds helpful messages about modules or actions
     * that can be taken by the user.
     */
    private JLabel buildInfoLabel() {
        final JLabel label = new JLabel(Program.TITLE);
        label.setFont(new Font(label.getFont().getFontName(), 0, 16));
        label.setHorizontalAlignment(SwingConstants.LEFT);
        label.setSize(label.getPreferredSize());
        label.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(final MouseEvent e) {
                if (0 != (e.getModifiers() & InputEvent.SHIFT_MASK)) {
                    StringSelection sel = new StringSelection(label.getText());
                    Toolkit.getDefaultToolkit().getSystemClipboard().setContents(sel, sel);
                }
            }
        });
        return label;
    }

    private Timer testIconsTimer = new Timer(1000, e->testIcons());
    private void testIcons() {
        //if (toolBar.getComponentCount() == 0) {
            toolBar.removeAll();
            File dir = new File("E:\\home\\jobs\\RNA\\RNAstructure\\alignment-editor\\src\\ur_rna\\StructureEditor\\resources\\images");
            for (File f : dir.listFiles()) {
                String name = f.getName();
                if (!name.endsWith(".png")) continue;
                try {
                    ImageIcon icon = new ImageIcon(ImageIO.read(f));
                    String title = name + fmt(" %sx%s", icon.getIconWidth(), icon.getIconHeight());
                    System.out.println(title);
                    if (icon.getIconHeight() > 32 || icon.getIconWidth() > 32)
                        continue;
                    JButton b = new JButton(icon);
                    b.setToolTipText(title);
                    toolBar.add(b);
                } catch (Exception ex) {
                    System.err.println("Error with " + name);
                    ex.printStackTrace();
                }
            }
            toolBar.add(new JButton("T" + System.currentTimeMillis()));
//        } else {
//            toolBar.removeAll();
//        }
        toolBar.revalidate();
        toolBar.repaint();
    }

    /**
     * Build the main RNAstructure toolbar.
     * <br><br>
     * Note that adding the same action listeners to the buttons as their
     * corresponding menu items uses some hard-coded numbers, so if the order of
     * the menu items or tool bar items needs to be changed, those numbers
     * (indexes) might need to be changed too.
     */
    private JToolBar buildToolBar() {
        // Create the toolbar.
        JToolBar bar = new JToolBar();
        bar.setFloatable(true);
        return bar;
    }
    private void createMainToolButtons() {
        // Add Main buttons (additional buttons added later by child windows)
        addToolButton(AppActions.SHOW_OPEN_FILE, FileMenuSection); // actionButton(AppActions.SHOW_OPEN_STRUCTURE, false));
        addToolButton(AppActions.SHOW_DASHBOARD, FileMenuSection);
        //toolBar.add(new MergeButton.Separator(FileMenuSection));

//        addToolButton(AppActions.SHOW_ABOUT_WINDOW, HelpMenuSection);
//        addToolButton(AppActions.EXIT_PROGRAM, HelpMenuSection);

        //addToolButton(new UiAction("Tool Window", e->desktop.add(new SliderToolWindow())), EditMenuSection);
    }
    private MergeButton addToolButton(Action a, int mergePos) {
        //return (MergeButton )toolBar.add(new MergeButton(a, mergePos));
        MergeButton b = new MergeButton(a, mergePos);
        b.setMnemonic(0);
        return (MergeButton )toolBar.add(b);
    }
//    /**
//     * Copy the ActionCommand and ActionListener from an existing menu item.
//     *
//     * @param menuItemName The name or title of the existing menu item.
//     * @param destination  The AbstractButton (e.g. toolbar button) to which the existing action will be copied.
//     */
//    private void copyAction(String menuItemName, AbstractButton destination) {
//        MergeItem m = mainMenus.findByName(menuItemName);
//        if (m == null)
//            AppLog.getDefault().error("There is no existing menu item with the name or title '" + menuItemName + "'.");
//        else
//            copyAction(m.getMenuItem(), destination);
//    }

    /**
     * Copy the ActionCommand and ActionListener from the source to the destination AbstractButton
     *
     * @param source      The AbstractButton (e.g. menu item) from which the existing action will be copied.
     * @param destination The AbstractButton (e.g. toolbar button) to which the existing action will be copied.
     */
    private void copyAction(AbstractButton source, AbstractButton destination) {
        destination.setActionCommand(source.getActionCommand());
        for (ActionListener l : source.getActionListeners())
            destination.addActionListener(l);
    }

    protected JDesktopPane getDesktop() { return desktop; }

    /**
     * Get the toolbar.
     *
     * @return The toolbar.
     */
    public JToolBar getToolBar() { return toolBar; }

    /**
     * Reset the info bar at the bottom of the main window to its default.
     */
    public void resetInfoLabel() {
        setInfoLabel("For help, press F1.");
    }
//    public Consumer<String> InfoLabelConsumer = s -> {
//        if (s == null)
//            resetInfoLabel();
//        else
//            setInfoLabel(s);
//    };

    /**
     * Set the info bar at the bottom of the main window to contain the
     * specified text.
     *
     * @param text The text to set on the info bar.
     */
    public void setInfoLabel(String text) {
        infoLabel.setText(text);
    }

    public ChildFrame getActiveChild() {
        JInternalFrame frame = super.getActiveChild();
        return frame instanceof ChildFrame ? (ChildFrame)frame : null;
    }
    public DrawWindow getActiveDrawing() {
        JInternalFrame frame = super.getActiveChild();
        return frame instanceof DrawWindow ? (DrawWindow)frame : null;
    }

    public void updateActiveFrameInfo() {
        JInternalFrame activeFrame = desktop.getSelectedFrame();
        String title;

        Collection<? extends Component> menus = null;
        Collection<? extends Component> tools = null;

        if (activeFrame == null) {
            title = "";
        } else {
            title = ObjTools.toStr(activeFrame.getTitle(), "").trim();
            if (activeFrame instanceof ChildFrame) {
                menus = ((ChildFrame) activeFrame).getMenus();
                tools = ((ChildFrame) activeFrame).getToolbarButtons();
            }
        }

        if (title.length() == 0)
            setTitle(Program.TITLE);
        else
            setTitle(Program.TITLE + " - " + title);

        rebuildMainMenu(menus);
        rebuildToolBar(tools);
        menuBar.revalidate();
        toolBar.revalidate();

        if (ObjTools.contains(desktop.getAllFrames(), f->f instanceof ChildFrame))
            hideDashboard();
        else
            showDashboard(false);

        this.repaint();
    }

//    private JButton createToolButton(final String text, final String tip, final String icon, final String action) {
//        JButton b = createToolButton(text, tip, icon, (Action)null);
//        b.setActionCommand(action);
//        return b;
//    }
//    private JButton createToolButton(final String text, final String tip, final String icon, final Action action) {
//        JButton btn = action == null ? new JButton() : new JButton(action);
//        int mnemonic = KeyMnemonic.getMnemonic(text);
//        if (mnemonic != 0) btn.setMnemonic(mnemonic);
//        btn.setText(KeyMnemonic.stripMnemonics(text));
//        if (tip != null)
//            btn.setToolTipText(tip);
//        if (icon != null)
//            btn.setIcon(getIcon(icon));
//        //if (action != null)
//        //    setAction(btn, action); //copy the ActionCommand and ActionListener from the existing menu item.
//        return btn;
//    }

//    private void setAction(final AbstractButton button, final Action action) {
//        //button.setActionCommand(action.name());
//        //button.addActionListener(program);
//        button.setAction(action);
//    }

    private void test() {
        JToolBar t = new JToolBar();
        t.setFloatable(true);
        t.add(new JButton("Hello there"));
        this.add(t, BorderLayout.EAST);

        //t.setOrientation(SwingConstants.HORIZONTAL);
        t.addAncestorListener(new AncestorListener() {
            @Override
            public void ancestorAdded(AncestorEvent event) {
                if (SwingUtilities.getWindowAncestor(t).equals(MainFrame.this)) {
                    System.out.println("...In Main Frame");
                } else {
                    System.out.println("...Maybe floating");
                }
            }

            @Override
            public void ancestorRemoved(AncestorEvent event) {
                if (SwingUtilities.getWindowAncestor(t).equals(MainFrame.this)) {
                    System.out.println("...In Main Frame");
                } else {
                    System.out.println("...Maybe floating");
                }
            }

            @Override
            public void ancestorMoved(AncestorEvent event) {
            }
        });
    }

    private Component[] originalToolbar;
    private void rebuildToolBar(final Collection<? extends Component> tools) {
        JToolBar bar = this.toolBar;
        if (originalToolbar == null)
            originalToolbar = bar.getComponents();

//        Component[] newTools = tools == null ? MergeManager.EMPTY_COMPONENT_ARRAY : tools.toArray(new Component[tools.size()]);
//        Component[] merged = new Component[originalToolbar.length + newTools.length];
//        System.arraycopy(originalToolbar, 0, merged, 0, originalToolbar.length);
//        System.arraycopy(newTools, 0, merged, originalToolbar.length, newTools.length);
        bar.removeAll();
        Components.addAll(bar, originalToolbar);
        if (tools != null) {
            removeMnemonics(tools);
            Components.addAll(bar, tools);
        }
        menuMerger.sort(bar);


        for(Component c : toolBar.getComponents())
            if (c instanceof AbstractButton) {
            AbstractButton b = (AbstractButton)c;
                String tip = b.getToolTipText();
                if (b.getIcon() != null && tip != null)
                    b.setHideActionText(true);
                if (tip != null && tip.indexOf('(') == -1) {
                    Action a = b.getAction();
                    if (a != null) {
                        KeyStroke ks = UiAction.getKeyStroke(a);
                        if (ks != null)
                            b.setToolTipText(tip + " (" + AcceleratorKey.toString(ks) + ")");
                    }
                }
            }
        // Add the open sequence button to the toolbar.

        //newSequence.setIcon( ImageGrabber.tryGetImageIcon( "NewSequence.gif" ) );
//        bar.add(createToolButton("Open Structure", "Load an RNA structure from a CT file.", "file-open", AppActions.SHOW_OPEN_STRUCTURE));
//        bar.add(createToolButton("Rotate", "Rotate selection.", "rotate2", "show-rotate"));
//        bar.add(createToolButton("Rotate", "Rotate selection.", "rotate3", "show-rotate"));
//        bar.add(createToolButton("Rotate", "Rotate selection.", "rotate-left", "show-rotate"));
//        bar.add(createToolButton("Rotate", "Rotate selection.", "rotate", "show-rotate"));


//        bar.add(actionButton(AppActions.SHOW_OPEN_STRUCTURE, false));
//        //bar.add(AppActions.SHOW_OPEN_STRUCTURE).setHideActionText(false);
//        bar.add(AppActions.OPEN_DEMO);
//        bar.add(AppActions.SHOW_ABOUT_WINDOW);
//        bar.add(AppActions.SHOW_OPEN_SEQ);
//        bar.add(AppActions.EXIT_PROGRAM);
//        bar.add(new UiAction("Banana", e->desktop.add(new SliderToolWindow())));

        // bar.setPreferredSize(new Dimension(200, 500));

//        JButton newSequence = new JButton();
//        newSequence.setToolTipText( "New Sequence" );
//        newSequence.setActionCommand( newSequence.getToolTipText() );
//        newSequence.setIcon( ImageGrabber.tryGetImageIcon( "NewSequence.gif" ) );
//        copyAction("File->New Sequence", newSequence); //copy the ActionCommand and ActionListener from the existing menu item.
//        bar.add(newSequence);
//
//        // Add the open sequence button to the toolbar.
//        JButton openSequence = new JButton();
//        openSequence.setToolTipText( "Open Sequence" );
//        openSequence.setActionCommand( openSequence.getToolTipText() );
//        openSequence.setIcon( ImageGrabber.tryGetImageIcon( "OpenSequence.gif" ) );
//        copyAction("File->Open Sequence", openSequence); //copy the ActionCommand and ActionListener from the existing menu item.
//        bar.add(openSequence);
//
//        // Add a separator.
//        bar.addSeparator();
//
//        // Add the save sequence (as) button to the toolbar.
//        // Note that this button initially doesn't have an action associated
//        // with it, and is also disabled.
//        // Only the SequenceDisplayWindow can activate or deactivate it.
//        JButton saveSequence = new JButton();
//        saveSequence.setToolTipText( "Save" );
//        saveSequence.setName( "Save" );
//        saveSequence.setActionCommand( "Save Sequence" );
//        saveSequence.setIcon( ImageGrabber.tryGetImageIcon( "Save.gif" ) );
//        saveSequence.setSelection( false );
//        saveSequence.setFocusable( false );
//        bar.add(saveSequence);
//
//        // Add a separator.
//        bar.addSeparator();
//
//        // Add the draw button to the toolbar.
//        JButton draw = new JButton();
//        draw.setToolTipText( "Draw" );
//        draw.setActionCommand( draw.getToolTipText() );
//        draw.setIcon( ImageGrabber.tryGetImageIcon( "Draw.gif" ) );
//        copyAction("File->Draw", draw); //copy the ActionCommand and ActionListener from the existing menu item.
////		draw.addActionListener(
////			menu.getMenu( 0 ).getItem( 5 ).getActionListeners()[0] );
//        bar.add(draw);
//
//        // Add the fold RNA single strand button to the toolbar.
//        JButton fold = new JButton();
//        fold.setToolTipText( "Fold RNA Single Strand" );
//        fold.setActionCommand( fold.getToolTipText() );
//        fold.setIcon( ImageGrabber.tryGetImageIcon( "FoldRNASingle.gif" ) );
//        copyAction("RNA->FoldSingleStrand", fold); //copy the ActionCommand and ActionListener from the existing menu item.
////		fold.addActionListener(
////			menu.getMenu( 1 ).getItem( 0 ).getActionListeners()[0] );
//        bar.add(fold);
//
//        // Add the RNA OligoWalk button to the toolbar.
//        JButton oligo = new JButton();
//        oligo.setToolTipText( "RNA OligoWalk" );
//        oligo.setActionCommand( oligo.getToolTipText() );
//        oligo.setIcon( ImageGrabber.tryGetImageIcon( "OligoWalk.gif" ) );
//        copyAction("RNA->OligoWalk", oligo);
////		oligo.addActionListener(
////			menu.getMenu( 1 ).getItem( 17 )
////			.getActionListeners()[0] );
//        bar.add(oligo);
//
//        // Add the RNA Dynalign button to the toolbar.
//        JButton dynalign = new JButton();
//        dynalign.setToolTipText( "RNA Dynalign" );
//        dynalign.setActionCommand( dynalign.getToolTipText() );
//        dynalign.setIcon( ImageGrabber.tryGetImageIcon( "Dynalign.gif" ) );
//        copyAction("RNA->Dynalign", dynalign);
////		dynalign.addActionListener(
////			menu.getMenu( 1 ).getItem( 12 )
////			.getActionListeners()[0] );
//        bar.add(dynalign);
    }
    private void removeMnemonics(final Collection<? extends Component> components) {
        for(Component c : components)
            if (c instanceof AbstractButton)
                ((AbstractButton) c).setMnemonic(0);
    }

//    private void convertTooltipToRolloverRecursive(MergeMenu mi) {
//        convertTooltipToRollover(mi);
//        for (MergeItem sub : mi.getSubItems())
//            if (sub instanceof MergeMenu)
//                convertTooltipToRolloverRecursive((MergeMenu)sub);
//            else
//                convertTooltipToRollover(sub);
//    }
//    private void convertTooltipToRollover(MergeItem mi) {
//        final JMenuItem item = mi.uiItem;
//        final String tip = item.getToolTipText();
//        if (!ObjTools.isEmpty(tip)) {
//            RolloverListener existing = ObjTools.firstOfType(item.getMouseListeners(), RolloverListener.class);
//            if (existing != null) return;
//            item.setToolTipText(null);
//            item.addMouseListener(new RolloverListener(this.InfoLabelConsumer, tip));
//        }
//    }

    private void rebuildMainMenu(Collection<? extends Component> custom) {
        menuMerger.reset();
        if (custom != null)
            menuMerger.merge(menuBar, custom.toArray(new Component[custom.size()]));
        menuMerger.sort(menuBar);

//        FileMenu f = (FileMenu)menuBar.getComponent(0);
//        f.setMnemonic('8');

        // MenuList m = custom == null ? mainMenus : Menus.merge(mainMenus, custom);
        // convertTooltipToRolloverRecursive(m);
        if (OSInfo.isMac())
            removeMenuMnemonics(menuBar.getComponents());
        else
            verifyMenuMnemonics(menuBar.getComponents());
        //m.buildMenu(menuBar);
    }

    private void removeMenuMnemonics(final Component[] list) {
        for (Component mi : list) {
            if (mi instanceof AbstractButton)
                ((AbstractButton)mi).setMnemonic(0);
            if (mi instanceof Container)
                removeMenuMnemonics(((Container)mi).getComponents());
        }
    }
    private void verifyMenuMnemonics(final Component[] list) {
        for (Component mi : list) {
            String s;
            if (mi instanceof AbstractButton) {
                AbstractButton b = (AbstractButton) mi;
                if (b.getMnemonic() == 0 && !Strings.isEmpty(s = b.getText()))
                    b.setMnemonic(s.charAt(0));
            }
            if (mi instanceof Container)
                verifyMenuMnemonics(((Container)mi).getComponents());
        }
    }

//    // public JMenuBar getMainMenu() {        return menuBar;    }
//    public MenuList getMainMenus() { return mainMenus; }

    /**
     * Invoked when a child frame has been closed.
     */
    @Override
    public void internalFrameClosed(final InternalFrameEvent e) {
        updateActiveFrameInfo();
        super.internalFrameClosed(e);
    }
    /**
     * Invoked when a child frame is activated.
     */
    @Override
    public void internalFrameActivated(final InternalFrameEvent e) {
        updateActiveFrameInfo();
        super.internalFrameActivated(e);
    }
    public JButton getToolBarItem(final String name) {
        if (name == null) return null;
        for (Component c : toolBar.getComponents()) {
            if (c instanceof JButton && name.equals(c.getName()))
                return (JButton) c;
        }
        for (Component c : toolBar.getComponents()) {
            if (c instanceof JButton && name.equals(((JButton) c).getText()))
                return (JButton) c;
        }
        for (Component c : toolBar.getComponents()) {
            if (c instanceof JButton && name.equals(((JButton) c).getActionCommand()))
                return (JButton) c;
        }
        if (log.isDebugEnabled())
            log.warn("Unknown Toolbar Item: " + name);
        return null;
    }

    @Override
    protected void addChild(final MdiChildFrame child) {
        WindowsMenu.enableTabNext(child); // disables Ctrl+TAB traversal key
        super.addChild(child);
    }
    public void show(final ChildFrame child) { show(child, null); }
    public void show(final ChildFrame child, final Point location) {
        if (ObjTools.indexOf(desktop.getAllFrames(), child) == -1) {
            //child.setVisible(false);
            addChild(child);
        }
        if (location != null)
            child.setLocation(location);

        boolean noMaximize = !child.isMaximizable();
        if (!noMaximize)
            for(JInternalFrame f : getDesktop().getAllFrames())
                if (f != child && f != dashboard && f.isVisible() && !f.isMaximum()) {
                    noMaximize = true;
                    break;
                }

        if (!noMaximize)
            try {
                child.setMaximum(true);
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        child.setVisible(true);
//        if (!child.hasFocus())
//            child.requestFocus();
        if (getDesktop().getSelectedFrame() != child)
            activateChild(child);
    }
    public void showCentered(final ChildFrame child) {
        Dimension desktopSize = desktop.getSize(), size = child.getSize();
        hideDashboard();
        show(child, new Point((desktopSize.width - size.width) / 2, (desktopSize.height - size.height) / 2));
    }
    public DashboardFrame getDashboard() {
        return dashboard;
    }
    public void showDashboard(boolean force) {
        if (!force&&!program.settings().ShowDashboard) return;
        if (!ObjTools.contains(desktop.getAllFrames(), dashboard))
            desktop.add(dashboard);
        dashboard.pack();
        dashboard.updateLayout();
        dashboard.show();
    }
    public void hideDashboard() {
        if (dashboard!=null) {
            dashboard.hide();
            desktop.remove(dashboard);
        }
    }

    public List<ChildFrame> getChildFrames() {
        JInternalFrame[] all = getDesktop().getAllFrames();
        ArrayList<ChildFrame> list = new ArrayList<>(all.length);
        for (JInternalFrame f : all)
            if (f instanceof ChildFrame)
                list.add((ChildFrame)f);
        return list;
    }
}
