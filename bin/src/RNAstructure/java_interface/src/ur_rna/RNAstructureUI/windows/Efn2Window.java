/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.menus.ConstraintsMenu;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.ui.FieldPanel.FilePanel;
import ur_rna.RNAstructureUI.ui.FileField;
import ur_rna.RNAstructureUI.ui.HTMLCheckBox;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.*;
import ur_rna.Utilities.swing.MergeMenu;

/**
 * A class responsible for initializing and running the efn2 (Energy Function
 * 2) module, which calculates the folding free energy of structures.
 *
 * @author Jessica S. Reuter
 */
public class Efn2Window
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public Efn2Window( String acid ) {
		super( acid, acid + " Efn2" );
	}

	@Override
	protected void processCommand(final CommandInfo ci) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "CT File" button, get a CT file,
		// initialize a data structure, and create a default output file name.
		if( ci.command.equals( "CT File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String file = StandardFileChooser.getOpenName(FileFilters.CT);
			if( file == null ) { return; }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildEfn2DataStructure( file, isRNA );
			if( !verifyBackendResult(result, "File: %s\nRNA: %s", file, isRNA) ) { return; }

			// Set the CT file name and the default output file name in the
			// input panel.
			// Then, enable the menus.
			String defaultOut = getOutputFile( file, "out" );
			files.setFile( 1, file );
			files.setFile( 2, defaultOut );
			getCustomMenus().enableMenus();
		}

		// If the action comes from the "Output File" button, try to select an
		// output file, and if one was selected set its name.
		else if( ci.command.equals( "Output File" ) ) {
			int index = 2;
			String file = StandardFileChooser.getSaveName(FileFilters.OUT, files.getFile(index));
			if( file != null ) { files.setFile( index, file ); }
		} else
			super.processCommand(ci);
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField ct = FileField.createDisabled( "CT File" ).inputFile(FileFilters.CT);
		FileField out = FileField.createEnabled( "Output File" ).setFilters(FileFilters.OUT);
		FilePanel files = new FilePanel( this, ct, out );
		files.setPanelWidth( 600 );
		files.makePanel();

		// Make the thermodynamic details box.
		HTMLCheckBox box =
			HTMLCheckBox.createEmptyBox( "Write Thermodynamic Details File" );

		// Add the components in their proper places.
		setGrid( 1, 1 );
		setFillHorizontal();
		placeComponent( 0, 0, files );
		setGrid( 2, 1 );
		setFillCenter();
		placeComponent( 0, 1, box );
		setAnchorNorth();
		setGrid( 1, 1 );
		makeStartButton( 2, 0 );
	}

	@Override
	protected boolean runMainCalculation() {

		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );
		HTMLCheckBox box = (HTMLCheckBox)getInputControl( 2 );
		files.saveRecent();

		// Get all data from the window.
		// If an error occurred while retrieving data, return.
		files.getFile( 1 );
		String out = files.getFile( 2 );
		if( files.isError() ) { return false; }
		boolean writeDetails = box.isSelected();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress(false);
		String result =
			backend.runEfn2( out, writeDetails );
		displayCalcError( result );
		if( !result.equals( "" ) ) { return false; }

		// Show a message saying that the calculation has finished.
		Dialogs.showMessage( "Efn2 Complete." );
		return true;
	}

	@Override
	protected MergeMenu[] createCustomMenus() {
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();
		return new ConstraintsMenu[]{ temperature };
	}
}
