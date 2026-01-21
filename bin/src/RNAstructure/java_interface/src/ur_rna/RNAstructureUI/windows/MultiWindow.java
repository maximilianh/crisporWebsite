/**
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.menus.ConstraintsMenu;
import ur_rna.RNAstructureUI.ui.*;
import ur_rna.RNAstructureUI.ui.FieldPanel.FilePanel;
import ur_rna.RNAstructureUI.ui.NumberField.IntegerField;
import ur_rna.RNAstructureUI.utilities.*;
import ur_rna.Utilities.swing.MergeMenu;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;

/**
 * A class responsible for creating and displaying an input window that takes
 * multiple arbitrary sequences and CT files as input.
 *
 * @author Jessica S. Reuter
 */
public abstract class MultiWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * The width of each column of the window (input controls in column 1,
	 * sequence set in column 2).
	 */
	protected final int columnWidth = 350;

	/**
	 * The list of drawing windows shown after the calculation is done.
	 */
	protected ArrayList<DrawingWindow> imageDialogs;

	/**
	 * The box that holds a particular window's unique options.
	 */
	protected Box options;

	/**
	 * The scroll pane holding the sequence set.
	 */
	private ScrollerPane scroll;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 * @param type   The type of window.
	 */
	protected MultiWindow( String acid, String type ) {
		super( acid, acid + " " + type );
		scroll.setSizes( columnWidth, 390 );
		imageDialogs = new ArrayList<DrawingWindow>();
	}

	/**
	 * Activate data input for this window.
	 */
	protected abstract void activate();

	/**
	 * Add a tuple to the calculation back end.
	 *
	 * @param seq   The sequence file.
	 * @param ct    The ct file.
	 */
	protected abstract void addBackendTuple( String seq, String ct );

	/**
	 * Build an input options box specific to this window.
	 */
	protected abstract void buildSpecificOptionsPanel();

	/**
	 * Delete a tuple from the calculation back end.
	 *
	 * @param index   The index of the tuple to delete.
	 */
	protected abstract void deleteBackendTuple( int index );

	@Override
	protected void processCommand(final CommandInfo ci) {
		// Get any input controls from the window that handle actions.
		FilePanel input = (FilePanel)getInputControl( 1 );
		JPanel seqSetPanel = (JPanel)getInputControl( 5 );
		FieldPanel delete =
			(FieldPanel) getInputControl( seqSetPanel, 3 )
			.getComponent(1);
		ScrollerPane scroll = (ScrollerPane)getInputControl( seqSetPanel, 2 );
		JTextArea area = (JTextArea)(scroll.getViewport().getView());

		// If the action comes from the "Add" button, add a sequence
		// to the sequence set.
		if( ci.command.equals( "ADD -->" ) ) {
			input.saveRecent();

			// Get the sequence and CT file names. If either of those files
			// were empty, return without adding them to the sequence set.
			String inFile = input.getFile( 1 );
			String outFile = input.getFile( 2 );
			if( inFile.equals( "" ) || outFile.equals( "" ) ) { return; }

			// Add the seq/ct tuple to the calculation back end.
			addBackendTuple( inFile, outFile );

			// Repaint the sequence set pane.
			invokeCommand( "Repaint",ci.event);

			// Clear the sequence and CT file names from the input fields.
			input.setFile( 1, "" );
			input.setFile( 2, "" );

			// If at least two sequences have been added, enable the menus.
			if( getNumSequences() >= 2 ) { getCustomMenus().enableMenus(); }
		}

		// If the action comes from the "CT File" button, try to select a CT
		// file, and if one was selected set its name.
		else if( ci.command.equals( "CT File" ) ) {
			int index = 2;
			String ct = StandardFileChooser.getSaveName(FileFilters.CT, input.getFile(index));
			if( !ct.equals( "" ) ) { input.setFile( index, ct ); }
		}

		// If the action comes from the delete button, delete the sequence/ct
		// tuple specified.
		else if( ci.command.equals( "Delete Sequence" ) ) {
			Integer index = ((IntegerField)delete.getField( 1 )).getValue();
			if( ( index < 1 ) || ( index > getNumSequences() ) ) {
				String message = "Cannot delete the selected sequence/ct pair.";
				Dialogs.showError( message );
			} else { deleteBackendTuple( index ); }
			invokeCommand( "Repaint", ci.event);
		}

		// If the command is to repaint the sequence set, do so.
		else if( ci.command.equals( "Repaint" ) ) {
			area.setText( getSequenceSetAsString() );
			((IntegerField)delete.getField( 1 ))
				.resetField( getNumSequences(), 0, getNumSequences() );
			area.repaint();
		}

		// If the action comes from the "Sequence File" button, get a sequence
		// file and create a default output file name.
		else if( ci.command.equals( "Sequence File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String seq = StandardFileChooser.getOpenName(FileFilters.Sequence);
			if( seq.equals( "" ) ) { return; }

			// Set the sequence file name and the default output file name in
			// the input panel.
			String defaultOut = getOutputFile( seq, "ct" );
			input.setFile( 1, seq );
			input.setFile( 2, defaultOut );
		} else
			super.processCommand(ci);
	}

	/**
	 * Get the number of sequences this window handles.
	 *
	 * @return   The number of sequences this window handles.
	 */
	protected abstract int getNumSequences();

	/**
	 * Get a representation of the sequence set as a string.
	 *
	 * @return   The sequence set as a string.
	 */
	protected abstract String getSequenceSetAsString();

	@Override
	protected void makeInputControls() {

		// Create a border builder for padding the input controls.
		BorderBuilder borderBuilder = new BorderBuilder();

		// Create the sequence set input panel.
		FileField seq = FileField.createDisabled( "Sequence File" ).inputFile(FileFilters.Sequence);
		FileField ct = FileField.createEnabled( "CT File" ).outputFile(FileFilters.CT);
		FilePanel input = new FilePanel( this, seq, ct );
		input.setPanelWidth( columnWidth );
		input.makePanel();

		// Create the "ADD" button.
		JButton addButton = new JButton( "ADD -->" );
		addButton.addActionListener( this );

		// Create the options panel that is unique to a specific window.
		options = Box.createVerticalBox();
		Dimension optionPanelSize = new Dimension( columnWidth, 305 );
		options.setPreferredSize( optionPanelSize );
		options.setMinimumSize( optionPanelSize );
		options.setMaximumSize( optionPanelSize );
		options.setSize( optionPanelSize );
		buildSpecificOptionsPanel();

		// Create the sequence set panel. This section is only in a block for
		// organizational purposes, since multiple components are placed on the
		// panel, but only the large panel is explicitly placed in the window.
		JPanel scrollPanel = new JPanel( new BorderLayout() );
		{

			// Create the label for the sequence set panel.
			JLabel label = new JLabel( "Sequence Set:" );
			label.setHorizontalAlignment( JLabel.CENTER );
			borderBuilder.makeTopBorder( 10, label );

			// Create the text area where the seq/ct pairs are displayed.
			JTextArea area = new JTextArea();
			area.setEditable( false );
			area.setDisabledTextColor( Color.black );
			area.setFont( new Font( Font.MONOSPACED, 0, 12 ) );
			area.setTabSize( 4 );

			// Create the scroll pane that holds the text area.
			scroll = new ScrollerPane( area, 0, 0 );
			scroll.setBarPolicies( ScrollerPane.ALWAYS, ScrollerPane.ALWAYS );
			borderBuilder.makeEqualBorder( 10, scroll );

			// Create the delete sequence panel.
			JButton deleteButton = new JButton( "Delete Sequence" );
			deleteButton.addActionListener( this );
			IntegerField del = new IntegerField( "", 0, 0, 0 );
			FieldPanel delete = new FieldPanel( del );
			delete.setPanelWidth( 100 );
			delete.makePanel();
			JPanel deletePanel = new JPanel( new BorderLayout() );
			borderBuilder.makeEqualBorder( 10, deletePanel );
			deletePanel.add( deleteButton, BorderLayout.WEST );
			deletePanel.add( delete );

			// Build the main scroll panel.
			scrollPanel.add( label, BorderLayout.NORTH );
			scrollPanel.add( scroll, BorderLayout.CENTER );
			scrollPanel.add( deletePanel, BorderLayout.SOUTH );
		}

		// Add the components in their proper places.
		setGrid( 1, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, input );
		setFillCenter();
		setPad( 40, 10 );
		placeComponent( 0, 1, addButton );
		setPad( 0, 0 );
		setFillHorizontal();
		placeComponent( 0, 2, options );
		makeStartButton( 0, 3 );
		setInsets( 0, 0, 0, 0 );
		setPad( 0, 0 );
		setGrid( 1, 4 );
		placeComponent( 1, 0, scrollPanel );
	}

	@Override
	protected boolean runMainCalculation() {

		// Check to make sure at least one sequence has been inputted. If no
		// sequences have been inputted, show an error and don't do any
		// calculations without sufficient data.
		if( getNumSequences() <= 0 ) {
			String message = "No file names given in sequence set.";
			Dialogs.showError( message );
			return false;
		}

		// Do the main calculation for multiple sequences.
		return runMainCalculationMultipleSequences();
	}

	/**
	 * Run the multiple sequence dependent module calculation in the back end.
	 */
	protected abstract boolean runMainCalculationMultipleSequences();

	@Override
	protected MergeMenu[] createCustomMenus() {
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();
		return new ConstraintsMenu[]{ temperature };
	}

	/**
	 * Show the image dialogs, if possible.
	 */
	protected void showImageDialogs() {

		// Determine if any dialogs ran into an error while being built.
		// If an error did occur, show an error dialog and return without
		// showing any dialogs.
		for( DrawingWindow window: imageDialogs ) {
			if( window.isError() ) {
				String message = "error creating image results window.";
				Dialogs.showError( message );
				return;
			}
		}

		// Show all image dialogs.
		for( DrawingWindow dialog: imageDialogs ) { dialog.showWindow(); }
	}
}
