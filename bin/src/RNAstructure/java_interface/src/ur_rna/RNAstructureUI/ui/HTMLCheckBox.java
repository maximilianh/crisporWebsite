/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.ui;

import javax.swing.*;

/**
 * A class that creates a check box whose text is HTML.
 *
 * @author Jessica S. Reuter
 */
public class HTMLCheckBox
	extends JCheckBox {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param text       The text of the box.
	 * @param selected   The flag that says if this check box is selected.
	 */
	private HTMLCheckBox( String text, boolean selected ) {
		super( "<html>" + text.replace( " ", "&nbsp;" ) );
		setName(text);
		// Add an empty border of 10 pixels around the entire check box.
		new BorderBuilder().makeEqualBorder( 10, this );

		// Set the check box selected or deselected.
		setSelected( selected );
	}

	/**
	 * Create an empty check box.
	 *
	 * @param text   The text of the box.
	 * @return       The empty check box.
	 */
	public static HTMLCheckBox createEmptyBox( String text ) {
		return new HTMLCheckBox( text, false );
	}

	/**
	 * Create a selected check box.
	 *
	 * @param text   The text of the box.
	 * @return       The selected check box.
	 */
	public static HTMLCheckBox createSelectedBox( String text ) {
		return new HTMLCheckBox( text, true );
	}
}
