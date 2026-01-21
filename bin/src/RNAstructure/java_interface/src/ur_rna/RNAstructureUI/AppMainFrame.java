/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI;

import ur_rna.RNAstructureUI.menus.*;
import ur_rna.RNAstructureUI.menus.MenuList;
import ur_rna.RNAstructureUI.utilities.ImageGrabber;
import ur_rna.RNAstructureUI.windows.AboutWindow;
import ur_rna.RNAstructureUI.windows.ChildWindow;
import ur_rna.RNAstructureUI.windows.InternalWindow;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.OSInfo;
import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.Strings;
import ur_rna.Utilities.annotation.NotNull;
import ur_rna.Utilities.swing.Menus;
import ur_rna.Utilities.swing.MergeManager;
import ur_rna.Utilities.swing.MergeMenu;

import javax.swing.*;
import javax.swing.border.BevelBorder;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import javax.swing.event.InternalFrameListener;
import java.awt.*;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.function.Consumer;

/**
 * A class that creates the application frame of the RNAstructure GUI.
 *
 * @author Jessica S. Reuter
 * @author Richard M. Watson
 */
public class AppMainFrame
        extends JFrame {
    private static final long serialVersionUID = 20160128;
    private final AppLog log = AppLog.getDefault();
    private static boolean useSimpleFileChooser;

    private final JMenuBar menuBar;
    private final JToolBar toolBar;
    private final JPanel statusBar;
    private final JLabel infoLabel;
    public final JDesktopPane desktop;
    private static AppMainFrame activeRoot;
    private MenuList mainMenus;
	private MergeManager menuMerger = new MergeManager(); // used to sort and merge the main menus and toolbar buttons with those from child windows.

    /**
     * Constructor.
     * Create the frame, set the icon and layout, add the toolbar, back panel,
     * and info label, add a focus listener that shows the default menus, then
     * set the size and location before setting the frame to be visible.
     */
    public AppMainFrame() {
        super("RNAstructure");
        AppMainFrame.activeRoot = this;
        String iconString = "Icon.gif";
        setIconImage(ImageGrabber.tryGetImage(iconString));
        setLayout(new BorderLayout());
		setJMenuBar(menuBar = new JMenuBar());

        mainMenus = buildMainMenu();
        toolBar = buildToolBar();
        desktop = buildDesktop();


        add(toolBar, BorderLayout.PAGE_START);
        add(desktop, BorderLayout.CENTER);
        statusBar = buildStatusBar();
        infoLabel = buildInfoLabel(statusBar);
        resetInfoLabel();
        updateActiveFrameInfo();
        setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        setSize(1024, 730);
        setLocationRelativeTo(null);
        //System.out.println("Before vis");
        setVisible(true);

        showStartupInfo();
        //System.out.println("After vis");
    }

    private void showStartupInfo() {
        String val = RNAstructure.getPref("initial_run_" + RNAstructure.VERSION);
        if (!ObjTools.asBool(val)) {
            RNAstructure.setPref("initial_run_" + RNAstructure.VERSION, "1");
            showAboutWindow();
            RNAstructure.createDocsFolder();
        }
    }

    public void showAboutWindow() {
        AboutWindow f = new AboutWindow();
        f.setVisible(false);
        desktop.add(f);
        Dimension desktopSize = desktop.getSize(), size = f.getSize();
        f.setLocation((desktopSize.width - size.width) / 2, (desktopSize.height - size.height) / 2);
        f.setVisible(true);
        f.requestFocus();
        f.initFocus();
    }

    //	@Override
//	protected void processWindowEvent(WindowEvent e) {
//		try {
//			System.out.println("Window Event "+ e.getID() + " - " + e.toString());
//			super.processWindowEvent(e);
//		} catch (Throwable t) {
//			System.out.println("error " + t.getMessage());
//			t.printStackTrace();
//		}
//	}

	/**
	 * Create the gray back panel that fills the main frame.
	 *
	 * @return   The back panel.
	 */
	private JDesktopPane buildDesktop() {
		JDesktopPane pane = new JDesktopPane();
		pane.setBackground( Color.GRAY );
		pane.removeAll();
		return pane;
	}

	private JPanel buildStatusBar() {
		JPanel p = new JPanel();
		p.setPreferredSize(new Dimension(100, 16));
		p.setLayout(new BoxLayout(p, BoxLayout.X_AXIS));
		p.setBorder(new BevelBorder(BevelBorder.LOWERED));
		this.add(p, BorderLayout.SOUTH);
		return p;
	}

	public static final int FileMenuSection = -100;
	public static final int EditMenuSection = -80;
	public static final int PredictMenuSection = -60;
	public static final int NucleicAcidMenuSection = -40;
	public static final int HelpMenuSection = 1000;
	private MenuList buildMainMenu() {
		MenuList list = new MenuList();
		// Create the File menu.
		list.add(new FileMenu(), FileMenuSection);
		list.add(new MergeMenu("&Edit"), EditMenuSection); //Add empty Edit menu as placeholder.
		list.add(new PredictMenu(), PredictMenuSection);
		// Create the RNA menu.
		list.add(NucleicAcidMenu.createRNA(), NucleicAcidMenuSection);
		// Create the DNA menu.
		list.add(NucleicAcidMenu.createDNA(), NucleicAcidMenuSection);
		// Create the Help menu.
		list.add(new HelpMenu(), HelpMenuSection);
		for(JMenu m : list)
			menuBar.add(m);
		return list;
	}
	private void addMainMenu(final MergeMenu menu, final int mergePos) {
		menu.setMergePos(mergePos);
		menuBar.add(menu);
	}

	/**
	 * Create the label at the bottom of the main frame.
	 * <br><br>
	 * This informative label holds helpful messages about modules or actions
	 * that can be taken by the user.
	 */
	private JLabel buildInfoLabel(JComponent parent) {
		final JLabel label = new JLabel("RNAstructure");
		label.setFont(new Font(label.getFont().getFontName(), 0, 16));
		label.setHorizontalAlignment(SwingConstants.LEFT);
        label.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(final MouseEvent e) {
                if (0 != (e.getModifiers() & InputEvent.SHIFT_MASK)) {
                    StringSelection sel = new StringSelection(label.getText());
                    Toolkit.getDefaultToolkit().getSystemClipboard().setContents(sel, sel);
                }
            }
        });
		add(label, BorderLayout.SOUTH);
		return label;
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
		createToolbar(bar);
		return bar;
	}

	/**
	 * Copy the ActionCommand and ActionListener from an existing menu item.
	 * @param menuItemName The name or title of the existing menu item.
	 * @param destination  The AbstractButton (e.g. toolbar button) to which the existing action will be copied.
	 */
	private void copyAction(String  menuItemName, AbstractButton destination)  {
		JMenuItem m = Menus.findByName(mainMenus, menuItemName);
		if (m == null)
			AppLog.getDefault().error("There is no existing menu item with the name or title '" + menuItemName + "'.");
		else
			copyAction(m, destination);
	}

	/**
	 * Copy the ActionCommand and ActionListener from the source to the destination AbstractButton
	 * @param source The AbstractButton (e.g. menu item) from which the existing action will be copied.
	 * @param destination  The AbstractButton (e.g. toolbar button) to which the existing action will be copied.
	 */
	private void copyAction(AbstractButton source, AbstractButton destination)  {
		destination.setActionCommand(source.getActionCommand());
		for (ActionListener l : source.getActionListeners())
			destination.addActionListener(l);
	}

	public static AppMainFrame getFrame() { //TODO rename alias
		return activeRoot;
	}
    /**
     * Get the top frame of the application.
     *
     * @return The top frame.
     */
	@NotNull public static AppMainFrame getAppRoot() {
	// RMW: previously this function would return the first Frame in the list of Windows.
		//  If any other dialog windows were open prior to the creation of the AppMainFrame,
		//  this would result in a ClassCastException.  Better to return a frame only
		//  only if it is a true AppMainFrame.
		//if (activeRoot != null)
		return activeRoot;
//
//		for (java.awt.Frame f : AppMainFrame.getFrames())
//			if (f instanceof AppMainFrame)
//				return (AppMainFrame)f;
//		return null;
    }

    /**
     * Get the most recent internal frame in the application.
     *
     * @return The top frame, or null if no frame exists.
     */
    public InternalWindow getMostRecentFrame() {
        return (InternalWindow) desktop.getSelectedFrame();
    }

	/**
	 * Get the toolbar.
	 *
	 * @return   The toolbar.
	 */
	public JToolBar getToolBar() { return toolBar; }

	/**
	 * Reset the info bar at the bottom of the main window to its default.
	 */
	public void resetInfoLabel() {
		setInfoLabel( "For help, press F1." );
	}
    public Consumer<String> InfoLabelConsumer = new Consumer<String>() {
        @Override
        public void accept(final String s) {
            if (s == null)
                resetInfoLabel();
            else
                setInfoLabel(s);
        }
    };
    /**
     * Set the info bar at the bottom of the main window to contain the
     * specified text.
     *
     * @param text The text to set on the info bar.
     */
    public void setInfoLabel(String text) {
        //RMW: The previous implementation caused bugs in several cases,
        //   related to fact that the index of the LJabel can change
        //   within the container upon certain user events, such as moving
        //   the toolbar etc
        //  BAD CODE: ((JLabel)getContentPane().getComponent( 2 )).setText(" " + text);
        infoLabel.setText(" " + text);
    }
    public static boolean getUseSimpleFileChooser() {
        return useSimpleFileChooser;
    }
    public static void setUseSimpleFileChooser(boolean value) {
        useSimpleFileChooser = value;
    }

	public void updateActiveFrameInfo() {
		JInternalFrame activeFrame = desktop.getSelectedFrame();
		String title;

		MenuList menus = null;
		if (activeFrame == null) {
			title = "";
		} else {
			title = ObjTools.toStr(activeFrame.getTitle(), "").trim();
			if (activeFrame instanceof ChildWindow)
				menus =  ((ChildWindow) activeFrame).getCustomMenus();
			else if (activeFrame instanceof InternalWindow)
				menus =  ((InternalWindow) activeFrame).getCustomMenus();
		}

		if (title.length() == 0)
			setTitle( "RNAstructure" );
		else
			setTitle( "RNAstructure - " +  title);

		rebuildMainMenu(menus);
		rebuildToolBar();
		this.repaint();
	}

	private Component[] keepToolBarItems = null;
	private void rebuildToolBar() {
		JToolBar bar = this.toolBar;
		if (keepToolBarItems != null) {
			for (Component c : bar.getComponents())
				if (!ObjTools.contains(keepToolBarItems, c))
					bar.remove(c);
		}
	}
	private void createToolbar(JToolBar bar) {
		// Add the new sequence button to the toolbar.
		JButton newSequence = new JButton();
		newSequence.setToolTipText( "New Sequence" );
		newSequence.setActionCommand( newSequence.getToolTipText() );
		newSequence.setIcon( ImageGrabber.tryGetImageIcon( "NewSequence.gif" ) );
		copyAction("File->New Sequence", newSequence); //copy the ActionCommand and ActionListener from the existing menu item.
		bar.add(newSequence);

		// Add the open sequence button to the toolbar.
		JButton openSequence = new JButton();
		openSequence.setToolTipText( "Open Sequence" );
		openSequence.setActionCommand( openSequence.getToolTipText() );
		openSequence.setIcon( ImageGrabber.tryGetImageIcon( "OpenSequence.gif" ) );
		copyAction("File->Open Sequence", openSequence); //copy the ActionCommand and ActionListener from the existing menu item.
		bar.add(openSequence);

		bar.addSeparator();

		// Add the save sequence (as) button to the toolbar.
		// Note that this button initially doesn't have an action associated
		// with it, and is also disabled.
		// Only the SequenceDisplayWindow can activate or deactivate it.
		JButton saveSequence = new JButton();
		saveSequence.setToolTipText( "Save" );
		saveSequence.setName( "Save" );
		saveSequence.setActionCommand( "Save Sequence" );
		saveSequence.setIcon( ImageGrabber.tryGetImageIcon( "Save.gif" ) );
		saveSequence.setEnabled( false );
		saveSequence.setFocusable( false );
		bar.add(saveSequence);

		// Add a separator.
		bar.addSeparator();

		// Add the draw button to the toolbar.
		JButton draw = new JButton();
		draw.setToolTipText( "Draw" );
		draw.setActionCommand( draw.getToolTipText() );
		draw.setIcon( ImageGrabber.tryGetImageIcon( "Draw.gif" ) );
		copyAction("File->Draw", draw); //copy the ActionCommand and ActionListener from the existing menu item.
//		draw.addActionListener(
//			menu.getMenu( 0 ).getItem( 5 ).getActionListeners()[0] );
		bar.add(draw);

		// Add the fold RNA single strand button to the toolbar.
		JButton fold = new JButton();
		fold.setToolTipText( "Fold RNA Single Strand" );
		fold.setActionCommand( fold.getToolTipText() );
		fold.setIcon( ImageGrabber.tryGetImageIcon( "FoldRNASingle.gif" ) );
		copyAction("RNA->FoldSingleStrand", fold); //copy the ActionCommand and ActionListener from the existing menu item.
//		fold.addActionListener(
//			menu.getMenu( 1 ).getItem( 0 ).getActionListeners()[0] );
		bar.add(fold);

		// Add the RNA OligoWalk button to the toolbar.
		JButton oligo = new JButton();
		oligo.setToolTipText( "RNA OligoWalk" );
		oligo.setActionCommand( oligo.getToolTipText() );
		oligo.setIcon( ImageGrabber.tryGetImageIcon( "OligoWalk.gif" ) );
		copyAction("RNA->OligoWalk", oligo);
//		oligo.addActionListener(
//			menu.getMenu( 1 ).getItem( 17 )
//			.getActionListeners()[0] );
		bar.add(oligo);

		// Add the RNA Dynalign button to the toolbar.
		JButton dynalign = new JButton();
		dynalign.setToolTipText( "RNA Dynalign" );
		dynalign.setActionCommand( dynalign.getToolTipText() );
		dynalign.setIcon( ImageGrabber.tryGetImageIcon( "Dynalign.gif" ) );
		copyAction("RNA->Dynalign", dynalign);
//		dynalign.addActionListener(
//			menu.getMenu( 1 ).getItem( 12 )
//			.getActionListeners()[0] );
		bar.add(dynalign);

		keepToolBarItems = bar.getComponents();
	}

	private void convertTooltipsToRollovers(JMenuItem item) {
		final String tip = item.getToolTipText();
		if (!Strings.isEmpty(tip)) {
			RolloverListener existing = ObjTools.firstOfType(item.getMouseListeners(), RolloverListener.class);
			if (existing != null) return;
			item.setToolTipText(null);
			item.addMouseListener(new RolloverListener(this.InfoLabelConsumer, tip));
		}
		if (item instanceof JMenu)
			for (JMenuItem sub : Menus.getMenuItems(item))
				convertTooltipsToRollovers(sub);
	}

	private void rebuildMainMenu(MenuList custom) {
		menuMerger.reset();
		if (custom != null)
			menuMerger.merge(menuBar, custom.toArray(new Component[custom.size()]));
		menuMerger.sort(menuBar);

		for (JMenu m : Menus.getMenus(menuBar))
			convertTooltipsToRollovers(m);

		if (OSInfo.isMac())
			removeMenuMnemonics(menuBar.getComponents());
		else
			verifyMenuMnemonics(menuBar.getComponents());
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

	public JMenuBar getMainMenu() {
		return menuBar;
	}
	public MenuList getMenus() {
		return mainMenus;
	}
	public void addChild(final JInternalFrame window) {
		desktop.add(window);
		window.addInternalFrameListener(frameListener);
	}
//	protected final InternalFrameListener frameListener = new InternalFrameAdapter() {
//		@Override
//		public void internalFrameClosing(final InternalFrameEvent e) { frameClosing(e); }
//		@Override
//		public void internalFrameClosed(final InternalFrameEvent e) { frameClosed(e); }
//	};
	protected InternalFrameListener frameListener = new InternalFrameAdapter() {
		@Override
		public void internalFrameClosed(final InternalFrameEvent e) { updateActiveFrameInfo();	}
		@Override
		public void internalFrameActivated(final InternalFrameEvent e) {
			updateActiveFrameInfo();
		}
	};

	public JButton getToolBarItem(final String name) {
		if (name == null) return null;
		for (Component c : toolBar.getComponents()) {
			if (c instanceof JButton && name.equals(c.getName()))
				return (JButton)c;
		}
		for (Component c : toolBar.getComponents()) {
			if (c instanceof JButton && name.equals(((JButton) c).getText()))
				return (JButton)c;
		}
		for (Component c : toolBar.getComponents()) {
			if (c instanceof JButton && name.equals(((JButton) c).getActionCommand()))
				return (JButton)c;
		}
		if (log.isDebugEnabled())
			AppLog.getDefault().warn("Unknown Toolbar Item: " + name);
		return null;
	}
}

