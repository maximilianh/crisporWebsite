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
import ur_rna.RNAstructureUI.ui.NumberField.IntegerField;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.*;
import ur_rna.Utilities.swing.MergeMenu;

/**
 * A class responsible for initializing and running the refold module, which
 * refolds a previous folding calculation from a save file.
 *
 * @author Jessica S. Reuter
 */
public class RefoldWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 * <br><br>
	 * Note that since the folding save file knows its nucleic acid type, the
	 * exact nucleic acid type this window uses doesn't matter, and the
	 * nucleic acid type is specified as an empty string.
	 */
	public RefoldWindow() {
		super( "", "Refold From Save File" );
	}

	@Override
	protected void processCommand(final CommandInfo ci) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );
		FieldPanel params = (FieldPanel)getInputControl( 2 );

		// If the action comes from the "CT File" button, try to select a CT
		// file, and if one was selected set its name.
		if( ci.command.equals( "CT File" ) ) {
			int index = 2;
			String file = StandardFileChooser.getSaveName(FileFilters.CT, files.getFile(index));
			if( file != null ) { files.setFile( index, file ); }
		}

		// If the action comes from the "Save File" button, get a SAV file,
		// initialize a data structure, and create a default output file name.
		else if( ci.command.equals( "Save File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String file = StandardFileChooser.getOpenName(FileFilters.FoldingSav);
			if( file == null ) { return; }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildRefoldSingleDataStructure( file );
			if( !verifyBackendResult(result, "File: %s", file) ) { return; }

			// Set the folding save file name and the default output file name
			// in the input panel.
			// Then, enable the menus.
			String defaultOut = getOutputFile( file, "ct" );
			files.setFile( 1, file );
			files.setFile( 2, defaultOut );
			getCustomMenus().enableMenus();

			// Reset the window size text field based on the sequence length.
			int size = backend.getRefoldWindowSize();
			((IntegerField)params.getField( 3 ))
				.resetField( size, 0, Integer.MAX_VALUE );
		} else
			super.processCommand(ci);
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField save = FileField.createDisabled( "Save File" ).inputFile(FileFilters.FoldingSav);
		FileField ct = FileField.createEnabled( "CT File" ).outputFile(FileFilters.CT);
		FilePanel files = new FilePanel( this, save, ct );
		files.setPanelWidth( 600 );
		files.makePanel();

		// Create the parameter panel.
		IntegerField energy =
			new IntegerField( "Max % Energy Difference", 10, 1 );
		IntegerField structures =
			new IntegerField( "Max Number of Structures", 20, 1 );
		IntegerField window = new IntegerField( "Window Size", 0, 0 );
		FieldPanel params = new FieldPanel( energy, structures, window );
		params.setPanelWidth( 250 );
		params.makePanel();

		// Add the components in their proper places.
		setGrid( 2, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, files );
		setGrid( 1, 1 );
		placeComponent( 0, 1, params );
		makeStartButton( 1, 1 );
	}

	@Override
	protected boolean runMainCalculation() {

		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );
		FieldPanel params = (FieldPanel)getInputControl( 2 );
		files.saveRecent();

		// Get the data from the file input panel.
		// If an error occurred while retrieving data, return.
		files.getFile( 1 );
		String ct = files.getFile( 2 );
		if( files.isError() ) { return false; }

		// Get the data from the parameters panel.
		Integer percent = ((IntegerField)params.getField( 1 )).getValue();
		Integer structures = ((IntegerField)params.getField( 2 )).getValue();
		Integer window = ((IntegerField)params.getField( 3 )).getValue();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress();
		String result =
			backend.runRefold( ct, percent, structures, window );
		displayCalcError( result );
		if( !result.equals( "" ) ) { return false; }

		// If the user wants to draw structures, draw them.
		return drawStructures(ct);
	}

	@Override
	protected MergeMenu[] createCustomMenus() {

		// Create the temperature menu.
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();

		// Create the menu that handles forced constraints.
		ConstraintsMenu forced = new ConstraintsMenu( backend );
		forced.addShowResetSection();

		// Create the maximum loop menu.
		ConstraintsMenu loop = new ConstraintsMenu( backend );
		loop.buildMaxLoopMenu();

		// Return the array of variable menus.
		return new ConstraintsMenu[]{ temperature, forced, loop };
	}
}
