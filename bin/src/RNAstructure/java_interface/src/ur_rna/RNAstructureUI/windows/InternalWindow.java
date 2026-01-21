/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.AppMainFrame;
import ur_rna.RNAstructureUI.RNAstructureBackendCalculator;
import ur_rna.RNAstructureUI.menus.MenuList;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.OSInfo;
import ur_rna.Utilities.swing.MergeMenu;

import javax.swing.*;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import javax.swing.event.InternalFrameListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyVetoException;
import java.beans.VetoableChangeListener;
import java.util.HashSet;

/**
 * The base class for all RNAstructure GUI dialogs and windows.
 * These are "children" of the JDesktopPane.
 * {@code AppMainFrame-->JDesktopPane-->[InternalWindow(s)...] }
 *
 * @author Jessica S. Reuter
 * @author Richard Watson
 */
public abstract class InternalWindow
	extends JInternalFrame implements ActionListener{
	private static final long serialVersionUID = 20120802;
	private final AppMainFrame mainFrame = AppMainFrame.getFrame();
	private final JDesktopPane desktop = mainFrame.desktop;
	private boolean disableClose;
	protected AppLog log = AppLog.getDefault();

	/**
	 * The back end calculator that's linked to this window.
	 */
	protected final RNAstructureBackendCalculator backend;

	/**
	 * Any custom menus that are associated to this window are added to the customMenus MenuList.
	 */
	private MenuList customMenus;

	protected HashSet<String> warnings = new HashSet<>();

	protected final InternalFrameListener frameListener = new InternalFrameAdapter() {
		@Override
		public void internalFrameClosing(final InternalFrameEvent e) {
			disableClose = !frameClosing(e);
		}
		@Override
		public void internalFrameClosed(final InternalFrameEvent e) { frameClosed(e); }
		@Override
		public void internalFrameActivated(final InternalFrameEvent e) { frameActivated(true, e);}
		@Override
		public void internalFrameDeactivated(final InternalFrameEvent e) { frameActivated(false, e);}
	};

	protected final VetoableChangeListener changeListener = new VetoableChangeListener() {
		/**
		 * This method gets called when a constrained property is changed.
		 *
		 * @param evt a <code>PropertyChangeEvent</code> object describing the
		 *            event source and the property that has changed.
		 * @throws PropertyVetoException if the recipient wishes the property
		 *                               change to be rolled back.
		 */
		@Override
		public void vetoableChange(final PropertyChangeEvent evt) throws PropertyVetoException {
			propertyChanged(evt.getPropertyName(), evt);
		}
	};

	protected void propertyChanged(final String name, final PropertyChangeEvent evt)
			throws PropertyVetoException  {
		if (name.equals("closed") && disableClose)
			throw new PropertyVetoException("Close operation cancelled.", evt);
	}
	/**
	 * Indicates that the InternalFrame is closing. Return true to allow the window to close, or return false to prevent it.
	 * @param e
	 * @return
	 */
	protected boolean frameClosing(final InternalFrameEvent e) { return true; }
	protected void frameClosed(final InternalFrameEvent e) { }
	protected void frameActivated(boolean active, final InternalFrameEvent e) { }

	/**
	 * Constructor.
	 */
	protected InternalWindow() {
		// Initialize the back end calculator and the dialog handler.
		backend = new RNAstructureBackendCalculator();

		// Set the frame defaults.
		setDefaultCloseOperation( JInternalFrame.DISPOSE_ON_CLOSE );
		setClosable( true );
		setFocusable( true );
		setIconifiable( true );
		setMaximizable( true );
		setResizable( true );

		// Add the window to the base frame.
		mainFrame.addChild( this );

		// Set the frame icon, if the look and feel allows for it.
		// If the OS is anything but a Mac, a frame icon can be set.
		if( !OSInfo.isMac() )
		    setFrameIcon( new ImageIcon( mainFrame.getIconImage() ) );

		// Add a focus listener so the window sets this window's associated
		// title and menu bar in the main frame.
		final JInternalFrame internal = this;
		addFocusListener( new FocusAdapter() { public void focusGained( FocusEvent e ) { onFocusGained(e); } });

		// Add a window closing listener that resets the main frame to its
		// original state when no more windows are open.
		addInternalFrameListener(frameListener);
//		addInternalFrameListener( new InternalFrameAdapter() {
//			public void internalFrameClosed( InternalFrameEvent e ) {
//				if( mainFrame.getMostRecentFrame() == null ) {
//					mainFrame.resetFrame();
//				}
//			}
//		});
	}

	/**
	 * An immutable object that represents a command invoked by some part of the program, e.g. in response to a user action.
	 */
	protected static class CommandInfo {
		/**
		 * The name or title of the command.
		 */
		public final String command;
		/**
		 * The source of the command -- e.g. a Component or other object class that invoked or caused the action/command.
		 */
		public final Object source;
		/**
		 * A reference to the original ActionEvent that invoked this command, if any.
		 */
		public final ActionEvent event;

		public CommandInfo(final String command) {
			this(command, null, null);
		}
		public CommandInfo(final ActionEvent event) { this(event.getActionCommand(), event.getSource(), event);  }
		public CommandInfo(final String command, final Object source) { this(command, source, null); }
		public CommandInfo(final String command, final Object source, final ActionEvent originalEvent) {
			this.command = command;
			this.source = source;
			this.event = originalEvent;
		}

		public String getCommand() { return command; }
		public int getId() { return event == null ? 0 : event.getID(); }
		public Object getSource() { return source; }
	}
	/**
	 * Perform an action specific to this window.
	 * @param ci        Information about the command.
	 */
	protected void processCommand(final CommandInfo ci) {
		// The command was not handled by the derived class, so show a warning.
		if (log.isDebugEnabled()) {
			if (!warnings.contains(ci.command)) {
				AppLog.getDefault().warn("Unprocessed UI command: '%s' (in %s; from %s)", ci.getCommand(),  describeSource(InternalWindow.this), describeSource(ci.source));
				warnings.add(ci.command);
			}
		}
	}

	protected void invokeCommand(String command) {
		processCommand(new CommandInfo(command, this));
	}
	protected void invokeCommand(String command, ActionEvent event) {
		processCommand(new CommandInfo(command, this, event));
	}
	protected void invokeCommand(String command, Object source) {
		processCommand(new CommandInfo(command, source));
	}
//	protected CommandInfo processCommand(String command) {
//		CommandInfo ci = new CommandInfo(command);
//		processCommand(command, ci);
//		return ci;
//	}

	protected void onFocusGained(FocusEvent e) {
//		AppMainFrame frame =
//				(AppMainFrame)getDesktopPane()
//						.getTopLevelAncestor();
//		setTitles( getTitle() );
//		frame.setJMenuBar( menuBar );
//		menuBar.revalidate();
//		menuBar.repaint();
//		frame.repaint();
//		desktop.getDesktopManager().activateFrame( internal );
	}

	@Override
	public void actionPerformed(final ActionEvent e) {
		processCommand(new CommandInfo(e));
	}

	private String describeSource(Object o) {
		if (o == null) return "NULL";
		String type = o.getClass().getSimpleName(),
				title = null;
		if (o instanceof Component)
			title = ((Component) o).getName();
		if (o instanceof JInternalFrame)
			title = ((JInternalFrame) o).getTitle();
		else if (o instanceof Dialog)
			title = ((Dialog) o).getTitle();
		else if (o instanceof AbstractButton)
			title = ((AbstractButton) o).getText();
		if (title == null)
			return type;
		return String.format("%s '%s'", type, title);
	}
//	public void actionPerformed( ActionEvent e ) {
//		Object source = e.getSource();
//		Component c = source instanceof Component ? (Component)source : null;
//		doActions( e.getActionCommand(), c, e );
//	}

	public void setCustomMenus(MenuList menus) { customMenus = menus; }
	public MenuList getCustomMenus() {
		if (customMenus == null) {
			customMenus = new MenuList();
			MergeMenu[] list = createCustomMenus();
			if (list!=null)
				customMenus.add(list);
		}
		return customMenus;
	}

	/**
	 * Set the caption text on this window (which is also shown on the main RNAstructure window caption bar).
	 *
	 * @param title   The new title.
	 */
	protected void setCaption( String title ) {
		super.setTitle(title);
		notifyUpdated();
		//AppMainFrame.getFrame().setTitle("RNAstructure -- " + title);
	}
	/**
	 * Call this if changes are made to the window that could require
	 * the JDesktopPane or the AppMainFrame to update.
	 * These could include changes to menus, title, etc.
	 */
	protected void notifyUpdated() {
        if (this.isActiveFrame())
		    mainFrame.updateActiveFrameInfo();
	}

	/**
	 * Override to create one or more custom menus that should be
	 * added to the main frame whenever this internal frame is active.
	 * @return an array of MergeMenu menus, or null if this
	 * window does not require any custom menus.
     */
	protected MergeMenu[] createCustomMenus() { return null; }

	protected boolean isActiveFrame() {
		return this == desktop.getSelectedFrame();
	}
	/**
	 * View the window.
	 */
	public void showWindow() {
		pack();
		setLocation(0, 0);
		setVisible(true);
		//if (!isActiveFrame())
		//	desktop.setSelectedFrame(this);
		//this.grabFocus();
		toFront();
		grabFocus();

		//desktop.setSelectedFrame(this);
//		menuBar.revalidate();
//		menuBar.repaint();
	}
}
