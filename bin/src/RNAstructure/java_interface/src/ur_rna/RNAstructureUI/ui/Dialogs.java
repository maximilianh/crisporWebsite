/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */
package ur_rna.RNAstructureUI.ui;

import ur_rna.RNAstructureUI.utilities.RequestFocusListener;

import javax.swing.*;
import java.awt.*;
import java.io.Serializable;

import static ur_rna.Utilities.ObjTools.toStr;
import static ur_rna.Utilities.Strings.escapeHTML;

/**
 * A class that handles simple dialogs used in the RNAstructure GUI.
 *
 * @author Jessica S. Reuter
 * @author Richard M. Watson
 *
 */
public abstract class Dialogs
	implements Serializable {
	private static final long serialVersionUID = 20120802;
	/**
	 * Show an error dialog.
	 *
	 * @param message   The message on the dialog.
	 */
	public static void showError(String message) {
		showError(message, "RNAstructure error");
	}
	public static void showError(String message, String title) {
		JOptionPane.showMessageDialog(null, fmtMsg(message), title, JOptionPane.ERROR_MESSAGE);
	}
	/**
	 * Show an input dialog that gets a string.
	 *
	 * @param title   The title (caption) on the dialog window.
	 * @return        The new text.
	 */
	public static String getInput(String title, String message, Object defaultValue) {
		Object result = JOptionPane.showInputDialog(null, message, title,
				JOptionPane.QUESTION_MESSAGE, null,null, defaultValue);
		return toStr(result, null);
	}
	/**
	 * Show an input dialog that gets a double within a certain range.
	 *
	 * @param field   The double input field.
	 * @return        The double selected from this field.
	 */
	public static Double getInput(NumberField.DoubleField field) {
		JPanel panel = new JPanel( new BorderLayout() );
		panel.add( new JLabel( field.getName() ), BorderLayout.NORTH );
		panel.add( field );
		field.setName("input");
		field.addAncestorListener(new RequestFocusListener());
		JOptionPane.showOptionDialog(
			null, panel, "RNAstructure Input", JOptionPane.OK_CANCEL_OPTION,
			JOptionPane.INFORMATION_MESSAGE, null, null, null );
		return field.getValue();
	}
	/**
	 * Show an input dialog that gets an integer within a certain range.
	 *
	 * @param field   The integer input field.
	 * @return        The integer selected from this field.
	 */
	public static Integer getInput(NumberField.IntegerField field) {
		JPanel panel = new JPanel( new BorderLayout() );
		panel.add( new JLabel( field.getName() ), BorderLayout.NORTH );
		panel.add( field );

		field.setName("input");
		field.addAncestorListener(new RequestFocusListener());
		JOptionPane.showOptionDialog(
			null, panel, "RNAstructure Input", JOptionPane.OK_CANCEL_OPTION,
			JOptionPane.INFORMATION_MESSAGE, null, null, null );
		return field.getValue();
	}

	/**
	 * Show a message dialog using String.format to create the message.
	 * @param msg The message on the dialog.
	 */
	public static void showMessage(String msg) {
		JOptionPane.showMessageDialog(null, fmtMsg(msg), "RNAstructure",  JOptionPane.INFORMATION_MESSAGE);
	}
	/**
	 * Show a message dialog using String.format to create the message.
	 * @param format The message on the dialog.
	 * @param args Arguments to String.format
	 */
	public static void showMessage(String format, Object... args) {
		showMessage( String.format(format,args));
	}
//	/**
//	 * Show a message dialog using String.format to create the message.
//	 * @param messageContainsHtml Whether or not the message intentionally includes HTML for formatting.
//	 * @param format The message on the dialog.
//	 * @param args Arguments to String.format
//	 */
//	public static void showMessage(boolean messageContainsHtml, String format, Object... args) {
//		showMessage(messageContainsHtml,  String.format(format,args));
//	}
	/**
	 * Show a dialog that gives the user one of three choices to make.
	 *
	 * @param message   The message on the dialog.
	 * @return          The response: "YES", "NO", or "CANCEL".
	 */
	public static String showYesNoCancel(Object message) {
		int choice = JOptionPane.showConfirmDialog(
			null, message, "RNAstructure",
			JOptionPane.YES_NO_CANCEL_OPTION );
		if( choice == JOptionPane.YES_OPTION ) { return "YES"; }
		else if( choice == JOptionPane.NO_OPTION ) { return "NO"; }
		else { return "CANCEL"; }
	}
	/**
	 * Show a dialog that gives the user one of three choices to make.
	 *
	 * @param msg   The message on the dialog.
	 * @return      The response: "YES", "NO", or "CANCEL".
	 */
	public static String showYesNoCancel(String msg) {
		StringBuilder builder = new StringBuilder( fmtMsg(msg) );
		return showYesNoCancel( builder );
	}
	/**
	 * Show a dialog that gives the user one of two choices to make.
	 *
	 * @param message   The message on the dialog.
	 * @return          The response: "OK" or "CANCEL".
	 */
	public static boolean showConfirm(Object message) {
		return JOptionPane.OK_OPTION == JOptionPane.showConfirmDialog(
				null, message, "RNAstructure", JOptionPane.OK_CANCEL_OPTION );
	}
	/**
	 * Show a dialog that gives the user one of two choices to make.
	 *
	 * @param msg   The message on the dialog.
	 * @return      The response: "OK" or "CANCEL".
	 */
	public static boolean showConfirm(String msg) {
		return JOptionPane.OK_OPTION == JOptionPane.showConfirmDialog(
				null, fmtMsg(msg), "RNAstructure", JOptionPane.OK_CANCEL_OPTION );
	}

	/**
	 * Show a warning dialog.
	 *
	 * @param message   The message on the dialog.
	 */
	public static void showWarning(Component message) {
		JOptionPane.showMessageDialog(null, message, "RNAstructure", JOptionPane.WARNING_MESSAGE );
	}

	/**
	 * Show a warning dialog.
	 *
	 * @param message   The message on the dialog.
	 */
	public static void showWarning(String message) {
		JOptionPane.showMessageDialog(null, fmtMsg(message), "RNAstructure", JOptionPane.WARNING_MESSAGE );
	}

	/**
	 * Show a warning dialog.
	 *
	 * @param message   The message on the dialog.
	 */
	public static void showWarning(String message, boolean scrollable) {
		if (scrollable)
			JOptionPane.showMessageDialog(null, createTextPanel(message), "RNAstructure", JOptionPane.WARNING_MESSAGE );
		else
			JOptionPane.showMessageDialog(null, fmtMsg(message), "RNAstructure", JOptionPane.WARNING_MESSAGE );
	}

	private static String fmtMsg(String message) {
		return message; // return htmlCenter(message)
	}

	/**
	 * Converts a simple text message to an html-formatted message
	 * by prepending {@code <html><center>} and escaping all html
	 * special characters within the message. However if the message
	 * already starts with {@code <html>}, it is returned as-is.
	 * @param msg The message to convert to HTML
     */
	public static String htmlCenter(String msg) {
		if (msg.substring(0,6).equals("<html>")) return msg;
		if (msg.substring(0,9).equals("<no-html>")) return msg.substring(9);
		return "<html><center>" + nl2br(escapeHTML(msg)) + "</center></html>";
	}

	public static String nohtml(String msg) {
		if (msg.substring(0,4).equals("<html>")) return msg;
		return "<no-html>" + msg;
	}

	public static String nl2br(final String s) {
		return s.replace("\n", "<br>").replace("\t", "&nbsp;&nbsp;&nbsp;&nbsp;");
	}

	public static JPanel createTextPanel(String content) {
		JPanel p = new JPanel(new BorderLayout());
		JTextArea text = new JTextArea();
		text.setBackground(p.getBackground());
		text.setForeground(p.getForeground());
		text.setFont(p.getFont());
		text.setLineWrap(true);
		text.setWrapStyleWord(true);
		JScrollPane scroll = new JScrollPane(text);
		scroll.setBorder(BorderFactory.createEmptyBorder(0,0,0,0));
		p.add(scroll, BorderLayout.CENTER);
		p.setPreferredSize(new Dimension(560, 300));
		text.setText(content);
		return p;
	}
}