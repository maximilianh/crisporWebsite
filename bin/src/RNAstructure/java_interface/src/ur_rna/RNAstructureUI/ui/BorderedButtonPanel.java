/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.ui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;

/**
 * A class that produces a panel which contains a button with a border
 * around it.
 *
 * @author Jessica S. Reuter
 */
public class BorderedButtonPanel
	extends JPanel {
	private static final long serialVersionUID = 20120802;

	/**
	 * Text Constructor.
	 * <br><br>
	 * By default, this constructor specifies a 5 pixel border all around the
	 * button.
	 *
	 * @param listener   The action listener connected to the button.
	 * @param text       The button text.
	 */
	public BorderedButtonPanel( ActionListener listener, String text ) {
		this( listener, text, 5 );
	}

	/**
	 * Border Constructor.
	 *
	 * @param listener   The action listener connected to the button.
	 * @param text       The button text.
	 * @param border     The border size around the button.
	 */
	public BorderedButtonPanel(
		ActionListener listener, String text, int border ) {

		// Create the button.
		JButton button = new JButton( text );
		button.setPreferredSize( new Dimension( 125, 30 ) );
		button.addActionListener( listener );

		// Add the button to the panel.
		add( button );

		// Add the border.
		BorderBuilder borderBuilder = new BorderBuilder();
		borderBuilder.makeEqualBorder( border, this );
	}
}
