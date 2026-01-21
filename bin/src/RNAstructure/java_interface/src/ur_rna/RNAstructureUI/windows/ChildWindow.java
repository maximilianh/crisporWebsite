/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.AppMainFrame;
import ur_rna.RNAstructureUI.menus.MenuList;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.OSInfo;

import javax.swing.*;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import javax.swing.event.InternalFrameListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyVetoException;
import java.beans.VetoableChangeListener;
import java.util.HashSet;

/**
 * The base class for all RNAstructure GUI dialogs and windows.
 * These are "children" of the JDesktopPane.
 * {@code AppMainFrame-->JDesktopPane-->[ChildWindow(s)...] }
 *
 * @author Jessica S. Reuter
 * @author Richard Watson (major revisions)
 */
public abstract class ChildWindow
	extends JInternalFrame implements FocusListener {
	private static final long serialVersionUID = 20120802;
	private boolean disableClose;
	protected final AppMainFrame mainFrame;
	protected MenuList customMenus;
	protected AppLog log = AppLog.getDefault();

	protected HashSet<String> warnings = new HashSet<>();

	protected final ActionListener commandListener = new ActionListener() {
		@Override
		public void actionPerformed(final ActionEvent e) {
			CommandInfo ci = new CommandInfo(e.getActionCommand(), e.getID(), e.getSource());
			processCommand(ci.getCommand(), ci);
			if (!ci.isHandled() && log.isDebugEnabled()) {
				if (!warnings.contains(ci.getCommand())) {
					AppLog.getDefault().warn("Unprocessed UI command: '%s' (in %s; from %s)", ci.getCommand(),  describeSource(ChildWindow.this), describeSource(ci.source));
					warnings.add(ci.getCommand());
				}
			}
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
	};
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
	protected ChildWindow() {
		// Set the frame defaults.
		super("ChildWindow", true, true, true, true);
		setDefaultCloseOperation( JInternalFrame.DISPOSE_ON_CLOSE );
//		setClosable( true );
//		setFocusable( true );
//		setIconifiable( false );
//		setMaximizable( false );
//		setResizable( false );

		addVetoableChangeListener(changeListener);
		mainFrame = AppMainFrame.getAppRoot();
		mainFrame.addChild(this);

		// Set the frame icon, if the look and feel allows for it.
		// If the OS is anything but a Mac, a frame icon can be set.
		if( !OSInfo.isMac() )
		    setFrameIcon( new ImageIcon( mainFrame.getIconImage() ) );

		addFocusListener(this);
		/**
		 * Listen for window events (see the internalFrameXXXX methods)
 		 */
		addInternalFrameListener(frameListener);
	}

	protected static class CommandInfo {
		private String command;
		private int id;
		private Object source;

		private boolean handled;
		public ActionEvent event;

		public CommandInfo(final String command, final int id, final Object source) {
			this.command = command;
			this.id = id;
			this.source = source;
		}
		public CommandInfo(final String command, final Object source) {
			this.command = command;
			this.source = source;
		}
		public CommandInfo(final String command) {
			this.command = command;
		}
		public String getCommand() { return command; }
		public int getId() { return id; }
		public Object getSource() { return source; }
		public boolean isHandled() { return handled; }
		public void markHandled() { this.handled = true; }
	}
	/**
	 * Do any actions specific to this window.
	 * @param ci        Information about the command.
	 */
	protected abstract void processCommand(String command, final CommandInfo ci);
	protected CommandInfo processCommand(String command) {
		CommandInfo ci = new CommandInfo(command);
		processCommand(command, ci);
		return ci;
	}

	/**
	 * Invoked when a component gains the keyboard focus.
	 *
	 * @param e
	 */
	@Override
	public void focusGained(final FocusEvent e) {
		//mainFrame.desktop.getDesktopManager().activateFrame( this );
		//updateMainFrame();
	}

	/**
	 * Call this if changes are made to the window that could require
	 * the JDesktopPane or the AppMainFrame to update.
	 * These could include changes to menus, title, etc.
	 */
	protected void updateMainFrame() {
		mainFrame.updateActiveFrameInfo();
	}

	public AppMainFrame getMainWindow() {
		return mainFrame;
	}

	public MenuList getCustomMenus() { return customMenus; }

	/**
	 * Invoked when a component loses the keyboard focus.
	 *
	 * @param e
	 */
	@Override
	public void focusLost(final FocusEvent e) { }
	/**
	 * View the window.
	 */
	public void showWindow() {
		pack();
		setLocation(0, 0);
		setVisible(true);
	}
}
