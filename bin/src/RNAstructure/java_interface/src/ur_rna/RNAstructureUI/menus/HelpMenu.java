/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.menus;

import ur_rna.RNAstructureUI.AppMainFrame;
import ur_rna.RNAstructureUI.RNAstructure;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.PathTools;

import javax.accessibility.AccessibleContext;
import javax.swing.*;
import javax.swing.text.JTextComponent;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.net.URI;
import java.util.HashMap;
import java.util.Objects;

import static ur_rna.Utilities.Strings.asBool;
import static ur_rna.Utilities.Strings.isEmpty;

/**
 * A class that creates a "Help" menu.
 *
 * @author Jessica S. Reuter
 */
public class HelpMenu
	extends MainMenu {
	private static final long serialVersionUID = 20120802;
	private static final String WEB_URL = "http://rna.urmc.rochester.edu/";
	private static final String HELP_URL = WEB_URL + "GUI/html/Contents.html";
	private static final String LOCAL_HELP_PAGE = "Contents.html";

	/**
	 * Constructor.
	 */
	public HelpMenu() {
		super( "Help" );
		addItem( "Program Help", "View the RNAstructure Manual.", "F1");
		addItem( "Online Help", "View the online RNAstructure help.", "*F1");

		if (PathTools.findLocalPath("examples") != null) {
			addSeparator();
			addItem("Extract Examples", "Extract example files to a folder of your choice.");
		}
		addSeparator();

		if (asBool(System.getProperty("show-gui-info-menu"))) {
			addItem("Show GUI Info", "Add mouse-over information for GUI Components.");
			addItem("Copy GUI Info Tree", "Copy the GUI Component hierarchy.");
			addSeparator();
		}

		addItem( "About RNAstructure...",
			"Display program information, version number, and copyright." );
	}

	@Override
	protected void onMenuAction(String command, final ActionEvent ev) {
		// If the command is to show the about screen, do so.
		if( command.startsWith( "About" ) ) {
			try {
				//Dialogs.showMessage( ImageGrabber.getImageLabel( "logo.png" ) );
				AppMainFrame.getFrame().showAboutWindow();
			} catch( Exception e ) {
				Dialogs.showError( "error showing about screen.\n" + AppLog.getErrorInfo(e));
			}
		}

		// Otherwise, show the help in a browser window.
		else if ("Program Help".equals(command)){
			String uri = HELP_URL;
			File page = PathTools.findLocalPath("manual/html/" + LOCAL_HELP_PAGE);
			if (page == null) page = PathTools.findLocalPath("manual/GUI/html/" + LOCAL_HELP_PAGE);
			if (page != null)
				uri = page.toURI().toString();
			browse(uri);
		}
		// Otherwise, show the help in a browser window.
		else if ("Online Help".equals(command)){
			browseHelp();
		}
		else if ("Extract Examples".equals(command)) {
			RNAstructure.ExtractExamples();
		}
		else if ("Show GUI Info".equals(command)) {
			for (Window w: Window.getWindows())
				addGuiInfo(w);
		}
		else if ("Copy GUI Info Tree".equals(command)) {
			String info = getGuiTree(KeyboardFocusManager.getCurrentKeyboardFocusManager().getActiveWindow());
			Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new StringSelection(info), null);
		}
	}
	private void addGuiInfo(final Component c) {
		if (c instanceof Window) {
			Window w = (Window)c;
			for (Window ow : w.getOwnedWindows())
				addGuiInfo(ow);
		}
		c.removeMouseListener(guiInfoListener);
		c.addMouseListener(guiInfoListener);
		if (c instanceof Container) {
			for (Component co : ((Container) c).getComponents())
				addGuiInfo(co);
		}
	}
	private final MouseAdapter guiInfoListener = new MouseAdapter() {
		@Override
		public void mouseClicked(final MouseEvent e) {
			showInfo(e);
		}
		@Override
		public void mouseEntered(final MouseEvent e) {
			showInfo(e);
		}
		void showInfo(MouseEvent e) {
			AppMainFrame.getFrame().setInfoLabel(getGuiInfo(e.getSource()));
		}
	};

	private String getGuiTree(Component c) {
		HashMap<Component,Integer> found = new HashMap<>();
		StringBuilder sb = new StringBuilder();
		appendGuiTree("", c, found, sb);
		return sb.toString();
	}
	private void appendGuiTree(String level, Component c, HashMap<Component,Integer> found, StringBuilder sb) {
		Integer id = found.get(c);
		if (id != null) {
			sb.append(level).append(id).append("\t(ref) ").append(c.getClass().getSimpleName());
			return;
		}

		id = found.size()+1;
		found.put(c, id);
		sb.append(level).append(id).append('\t').append(getGuiInfo(c)).append('\n');

		level = level + "..";
		if (c instanceof Window) {
			for (Window ow : ((Window)c).getOwnedWindows())
				appendGuiTree(level, ow, found, sb);
		}
		if (c instanceof Container) {
			for (Component co : ((Container) c).getComponents())
				appendGuiTree(level, co, found, sb);
		}
	}

	private String getGuiInfo(Object source) {
		if (source == null)
			return "Source: NULL";
		if (!(source instanceof Component))
			return "Source Type: " + source.getClass().getSimpleName();
		return getComponentInfo((Component)source);
	}

	public static String getComponentInfo(Component c) {
		StringBuilder sb = new StringBuilder();

		sb.append(c.getClass().getSimpleName());

		appendInfo(sb, "Name",c.getName());

		if (c instanceof JComponent) appendInfo(sb, "ToolTip",((JComponent) c).getToolTipText());

		if (c instanceof JTextComponent) appendInfo(sb, "InputText",((JTextComponent) c).getText());
		else if (c instanceof TextComponent) appendInfo(sb, "InputText",((TextComponent) c).getText());

		if (c instanceof AbstractButton ) appendInfo(sb, "Caption",((AbstractButton) c).getText()); // includes JButton, JMenuItem
		else if (c instanceof Button) appendInfo(sb, "Caption",((Button) c).getLabel());
		else if (c instanceof JLabel) appendInfo(sb, "Caption",((JLabel) c).getText());
		else if (c instanceof Label) appendInfo(sb, "Caption",((Label) c).getText());
		else if (c instanceof Frame) appendInfo(sb, "Caption",((Frame) c).getTitle());
		else if (c instanceof Dialog) appendInfo(sb, "Caption",((Dialog) c).getTitle());

		if (c instanceof JList) appendInfo(sb, "SelItem",(""+((JList) c).getSelectedValue()));
		else if (c instanceof JComboBox) appendInfo(sb, "SelItem",""+(((JComboBox) c).getSelectedItem()));

		AccessibleContext ac = c.getAccessibleContext();
		if (ac != null) {
			if (!Objects.equals(c.getName(),ac.getAccessibleName()))
				appendInfo(sb, "AccessName",ac.getAccessibleName());
			appendInfo(sb, "AccessDesc",ac.getAccessibleDescription());
		}

		return sb.toString();
	}

	public static void appendInfo(StringBuilder sb, String name, String value) {
		if (value != null)
			sb.append(' ').append(name).append(": \"").append(value).append('\"');
	}

	public static void browseHelp() {
		browse(HELP_URL);
	}
	public static void browseWebsite() {
		browse(WEB_URL);
	}
	public static void browse(String uri) {
		try {
			if (isEmpty(uri)) uri = WEB_URL;
			Desktop.getDesktop().browse(new URI(uri));
		} catch( Exception e ) {
			Dialogs.showError( "error showing webpage in browser: " + uri);
		}
	}
}
