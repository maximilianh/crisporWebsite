/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.menus.ConstraintsMenu;
import ur_rna.RNAstructureUI.ui.FieldPanel;
import ur_rna.RNAstructureUI.ui.FieldPanel.FilePanel;
import ur_rna.RNAstructureUI.ui.FileField;
import ur_rna.RNAstructureUI.ui.HTMLCheckBox;
import ur_rna.RNAstructureUI.ui.NumberField.IntegerField;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.*;
import ur_rna.Utilities.swing.MergeMenu;


/**
 * A class responsible for initializing and running the Fold module, which
 * folds a strand of nucleic acids into their likely lowest free energy
 * structure.
 *
 * @author Jessica S. Reuter
 */
public class FoldSingleWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public FoldSingleWindow( String acid ) {
		super( acid, acid + " Single Strand Fold" );
	}

	@Override
	protected void processCommand(final CommandInfo ci) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );
		FieldPanel params = (FieldPanel)getInputControl( 3 );

		// If the action comes from the "CT File" button, try to select a CT
		// file, and if one was selected set its name.
		if( ci.command.equals( "CT File" ) ) {
			int index = 2;
			String file = StandardFileChooser.getSaveName(FileFilters.CT, files.getFile(index));
			if( file != null ) { files.setFile( index, file ); }
		}

		// If the action comes from the "Sequence File" button, get a sequence
		// file, initialize a data structure, and create a default output file
		// name.
		else if( ci.command.equals( "Sequence File" ) ) {
			// Attempt to select the file.
			// If no file was selected, return.

			String file = StandardFileChooser.getOpenName(FileFilters.Sequence);
//			if (event instanceof FileSelectedEvent)
//				file = ((FileSelectedEvent) event).getFile().toString();
//			else
//				file

			if( file == null ) { return; }
			// Fill in the window data.
			setSequenceFile(file);
		}

		// If the command is to fill in the window's data automatically,
		// initialize a data structure, and create a default output file name.
		else if( ci.command.startsWith( "Sequence File Auto" ) ) {
			// Get the file whose data is being filled in.
			setSequenceFile(ci.command.split( ";" )[1]);
		} else
			super.processCommand(ci);
	}

	public void setSequenceFile(String file) {
		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );
		FieldPanel params = (FieldPanel)getInputControl( 3 );

		// Create a data structure.
		// If an error occurred creating the data structure, show an error
		// and return.
		String result = backend.buildFoldDataStructure( file, isRNA );
		if( !verifyBackendResult(result, "File: " + file + "\nRNA: " + isRNA) ) { return; }

		// Reset the window size text field based on the sequence length.
		int size = backend.getFoldWindowSize();
		((IntegerField)params.getField( 3 )).resetField( size, 0, Integer.MAX_VALUE );

		// Set the sequence file name and the default output file name in
		// the input panel.
		// Then, enable the menus.
		String defaultOut = getOutputFile( file, "ct" );
		files.setFile( 1, file );
		files.setFile( 2, defaultOut );
		getCustomMenus().enableMenus();
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField seq = FileField.createDisabled( "Sequence File" ).inputFile(FileFilters.Sequence);
		FileField ct = FileField.createEnabled( "CT File" ).outputFile(FileFilters.CT);
		FilePanel files = new FilePanel( this, seq, ct );
		files.setPanelWidth( 600 );
		files.makePanel();

		// Create the box that allows the user to create a save file.
		HTMLCheckBox box =
			HTMLCheckBox.createEmptyBox( "Generate Save File" );

		// Create the parameter panel.
		IntegerField energy =
			new IntegerField( "Max % Energy Difference", 10, 0 );
		IntegerField structures =
			new IntegerField( "Max Number of Structures", 20, 1 );
		IntegerField window = new IntegerField( "Window Size", 0, 0 );
		FieldPanel params = new FieldPanel( energy, structures, window );
		params.setPanelWidth( 300 );
		params.makePanel();

		// Add the components in their proper places.
		setGrid( 2, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, files );
		setFillCenter(); 
		placeComponent( 0, 1, box );
		setFillHorizontal();
		setGrid( 1, 1 );
		placeComponent( 0, 2, params );
		makeStartButton( 1, 2 );
	}

	@Override
	protected boolean runMainCalculation() {

		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );
		HTMLCheckBox box = (HTMLCheckBox)getInputControl( 2 );
		FieldPanel params = (FieldPanel)getInputControl( 3 );
		files.saveRecent();

		// Get the data from the file input panel.
		// If an error occurred while retrieving data, return.
		files.getFile( 1 );
		String ct = files.getFile( 2 );
		if( files.isError() ) { return false; }

		// Get the data from the parameters panel.
		Integer percent = ((IntegerField)params.getField( 1 )).getValue();
		Integer structs = ((IntegerField)params.getField( 2 )).getValue();
		Integer window = ((IntegerField)params.getField( 3 )).getValue();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress(true, true);
		String result = backend.runFold( ct, percent, structs, window, box.isSelected() );
		if( !displayCalcError( result ))  return false;
		// If the user wants to draw structures, draw them.
		return drawStructures(ct);
	}
	/**
	 * Set the window's data automatically.
	 *
	 * @param file   The input file.
	 */
	public void setDataAutomatically( String file ) {
		invokeCommand("Sequence File Auto;" + file, this);
	}

	@Override
	protected MergeMenu[] createCustomMenus() {

		// Create the temperature menu.
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();

		// Create the menu that handles forced constraints.
		ConstraintsMenu forced = new ConstraintsMenu( backend );
		forced.addGeneralSection();
		forced.addMaxPairingDistanceSection();
		forced.addSHAPESection();
		forced.addShowResetSection();
		forced.addSaveRestoreSection();

		// Create the maximum loop menu.
		ConstraintsMenu loop = new ConstraintsMenu( backend );
		loop.buildMaxLoopMenu();

		// Return the array of variable menus.
		return new ConstraintsMenu[]{ temperature, forced, loop };
	}
}
