/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.ui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.util.ArrayList;

/**
 * A class that creates a panel containing a group of radio buttons.
 *
 * @author Jessica S. Reuter
 */
public class RadioButtonPanel
	extends JPanel {
	private static final long serialVersionUID = 20120802;

	/**
	 * The button group connected to this panel.
	 */
	private ButtonGroup group;

	/**
	 * The button group, as a list, for ease of access to each.
	 */
	private ArrayList<JRadioButton> group2;

	/**
	 * The names of the buttons in the group.
	 */
	private String[] names;

	/**
	 * Constructor.
	 *
	 * @param title      The title of the panel.
	 * @param layout     The layout of the panel, horizontal or vertical.
	 * @param listener   The action listener connected to the buttons.
	 * @param names      The names of the buttons.
	 */
	private RadioButtonPanel(
		String title, String layout, ActionListener listener,
		String... names ) {

		// Determine the number of rows and columns in the panel.
		int rows = ( layout.equals( "Horizontal" ) ) ? 1 : names.length;
		int cols = ( layout.equals( "Vertical" ) ) ? 1 : names.length;

		// Set the panel's layout.
		setLayout( new GridLayout( rows, cols ) );

		// Create its titled border, if necessary
		if( title != null ) {
			new BorderBuilder().makeTitledBorder( title, this );
		}

		// Create the button group and determine how many buttons should be
		// in it. Also, create the button list and save the button names.
		group = new ButtonGroup();
		int numButtons = names.length;
		group2 = new ArrayList<JRadioButton>( numButtons );
		this.names = names;

		// Create a radio button for each name, and add each button to the
		// panel and to the button group.
		for( int i = 1; i <= numButtons; i++ ) {
			JRadioButton button = new JRadioButton( names[i-1] );
			button.setActionCommand( names[i-1] );
			if( listener != null ) { button.addActionListener( listener ); }
			if( i == 1 ) { button.setSelected( true ); }
			group.add( button );
			group2.add( button );
			add( button );
		}
	}

	/**
	 * Get the index of the selected button from this panel.
	 *
	 * @return   The button index, one-indexed.
	 */
	public int getSelectedIndex() {
		for( int i = 1; i <= names.length; i++ ) {
			if( names[i-1].equals( getSelectedName() ) ) { return i; }
		}
		return 0;
	}

	/**
	 * Get the command of the selected button from this panel.
	 *
	 * @return   The button name.
	 */
	public String getSelectedName() {
		return group.getSelection().getActionCommand();
	}

	/**
	 * Create a button panel with a horizontal layout.
	 *
	 * @param names   The names of the buttons.
	 */
	public static RadioButtonPanel makeHorizontal( String[] names ) {
		return new RadioButtonPanel( null, "Horizontal", null, names );
	}

	/**
	 * Create a button panel with a vertical layout.
	 *
	 * @param names   The names of the buttons.
	 */
	public static RadioButtonPanel makeVertical( String[] names ) {
		return new RadioButtonPanel( null, "Vertical", null, names );
	}

	/**
	 * Create a button panel with a title and horizontal layout.
	 *
	 * @param title   The title of the panel.
	 * @param names   The names of the buttons.
	 */
	public static RadioButtonPanel makeHorizontal(
		String title, String... names ) {
		return new RadioButtonPanel( title, "Horizontal", null, names );
	}

	/**
	 * Create a button panel with a title and vertical layout.
	 *
	 * @param title   The title of the panel.
	 * @param names   The names of the buttons.
	 */
	public static RadioButtonPanel makeVertical(
		String title, String... names ) {
		return new RadioButtonPanel( title, "Vertical", null, names );
	}

	/**
	 * Create a button panel with a title, horizontal layout, and an action
	 * listener associated with the buttons.
	 *
	 * @param title      The title of the panel.
	 * @param listener   The listener that controls the actions done by buttons.
	 * @param names      The names of the buttons.
	 */
	public static RadioButtonPanel makeHorizontal(
		String title, ActionListener listener, String... names ) {
		return new RadioButtonPanel( title, "Horizontal", listener, names );
	}

	/**
	 * Create a button panel with a title, vertical layout, and an action
	 * listener associated with the buttons.
	 *
	 * @param title      The title of the panel.
	 * @param listener   The listener that controls the actions done by buttons.
	 * @param names      The names of the buttons.
	 */
	public static RadioButtonPanel makeVertical(
		String title, ActionListener listener, String... names ) {
		return new RadioButtonPanel( title, "Vertical", listener, names );
	}

	/**
	 * Set a particular button in this group to be selected.
	 *
	 * @param name   The name of the button to be selected.
	 */
	public void setSelectedButton( String name ) {

		// Get the button to select.
		JRadioButton selectedButton = null;
		for( JRadioButton button: group2 ) {
			if( button.getText().equals( name ) ) {
				selectedButton = button;
				break;
			}
		}

		// Set the button selected in both groups.
		group.setSelected( selectedButton.getModel(), true );
		for( JRadioButton button: group2 ) {
			boolean doSelect =  button.getText().equals( name );
			button.setSelected( doSelect );
		}
	}
}
