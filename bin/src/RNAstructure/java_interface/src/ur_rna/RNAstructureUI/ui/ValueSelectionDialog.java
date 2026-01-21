/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */
package ur_rna.RNAstructureUI.ui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
/**
 * A class which creates a dialog that selects values from labeled text fields.
 *
 * @author Jessica S. Reuter
 */
public abstract class ValueSelectionDialog
	extends JDialog {
	private static final long serialVersionUID = 20120802;
	/**
	 * The array of number fields this dialog holds.
	 */
	protected NumberField[] fields;
	/**
	 * Constructor.
	 *
	 * @param title   The title of the dialog.
	 */
	protected ValueSelectionDialog( String title ) {
		setTitle( title );
		setLayout( new BorderLayout() );
	}
	/**
	 * Build the value selection dialog.
	 *
	 * @param selectText   The text on the selection button.
	 * @param fields       The list of number fields shown on this panel.
	 */
	public void buildDialog( String selectText, NumberField... fields ) {
		// Save the fields.
		this.fields = fields;
		// Create the fields panel.
		int numFields = fields.length;
		JPanel fieldPanel = new JPanel( new GridLayout( 1, numFields ) );
		for( int i = 1; i <= numFields; i++ ) {
			// Create the field label.
			JLabel label = new JLabel( fields[i-1].getName() );
			label.setHorizontalAlignment( JLabel.CENTER );
			new BorderBuilder().makeTopBorder( 5, label );
			// Set the panel size and put it in a panel.
			int height = fields[i-1].getPreferredSize().height;
			Dimension size = new Dimension( 50, height );
			fields[i-1].setPreferredSize( size );
			JPanel fieldWrapper = new JPanel();
			fieldWrapper.add( fields[i-1] );
			// Put the label and field together into a panel, then add that
			// panel to the fields panel.
			JPanel mini = new JPanel( new BorderLayout() );
			mini.add( label, BorderLayout.NORTH );
			mini.add( fieldWrapper );
			fieldPanel.add( mini );
		}
		// Create the selection button.
		JPanel selectPanel = new JPanel();
		JButton select = new JButton( selectText );
		select.setPreferredSize( new Dimension( 125, 30 ) );
		select.setActionCommand( selectText + " " + getTitle() );
		select.addActionListener( createSelectionAction() );
		selectPanel.add( select );
		// Create the "Cancel" button.
		JPanel cancelPanel = new JPanel();
		JButton cancel = new JButton( "Cancel" );
		cancel.setPreferredSize( new Dimension( 125, 30 ) );
		cancel.addActionListener( new ActionListener() {
			public void actionPerformed( ActionEvent e ) {
				dispose();
			}
		});
		cancelPanel.add( cancel );
		// Add borders to buttons.
		BorderBuilder borderBuilder = new BorderBuilder();
		borderBuilder.makeEqualBorder( 5, selectPanel );
		borderBuilder.makeEqualBorder( 5, cancelPanel );
		// Add the buttons to a panel.
		JPanel buttons = new JPanel( new GridLayout( 1, 2 ) );
		buttons.add( cancelPanel );
		buttons.add( selectPanel );
		// Create any extra objects in a panel, if more things need to be
		// shown on the dialog.
		JPanel additionalPieces = createAdditionalControls();
		// Add the component panels to the dialog and show it.
		add( fieldPanel, BorderLayout.NORTH );
		if( additionalPieces != null ) {
			add( additionalPieces, BorderLayout.CENTER );
		}
		add( buttons, BorderLayout.SOUTH );
		pack();
		setLocationRelativeTo( null );
		setVisible( true );
	}
	/**
	 * Create any additional controls that need to show up on this dialog.
	 * <br><br>
	 * This method does not need to be overridden; it's only to provide an
	 * ability to add more to this dialog when the situation calls for it.
	 */
	public JPanel createAdditionalControls() { return null; }
	/**
	 * Create the action that occurs once values are selected.
	 */
	public abstract ActionListener createSelectionAction();
}
