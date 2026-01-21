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
import ur_rna.RNAstructureUI.ui.RadioButtonPanel;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.*;
import ur_rna.Utilities.swing.MergeMenu;

import java.awt.*;

/**
 * A class responsible for initializing and running the OligoScreen module,
 * which calculates thermodynamic properties for a list of oligonucleotides.
 *
 * @author Jessica S. Reuter
 */
public class OligoScreenWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * The oligomer chemistry.
	 */
	private static String oligomerChemistry;

	/**
	 * Constructor.
	 * <br><br>
	 * This module has DNA hard-coded as its nucleic acid type to start with,
	 * but the nucleic acid type is dynamic and can be changed in the window.
	 */
	public OligoScreenWindow() {
		super( "", "OligoScreen" );
		isRNA = false;
	}

	@Override
	protected void processCommand(final CommandInfo ci) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "Oligomer List" button, get a list
		// file, initialize a data structure, and create a default output
		// file name.
		if( ci.command.equals( "Oligomer List" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String file = StandardFileChooser.getOpenName(FileFilters.List);
			if( file == null ) { return; }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildOligoScreenDataStructure();
			if( !verifyBackendResult(result, null) ) { return; }

			// Set the list file name and the default output file name in the
			// input panel.
			// Then, enable the menus.
			String defaultOut = getOutputFile( file, "rep" );
			files.setFile( 1, file );
			files.setFile( 2, defaultOut );
			getCustomMenus().enableMenus();
		}

		// If the action comes from the "Report File" button, try to
		// select an output file, and if one was selected set its name.
		else if( ci.command.equals( "Report File" ) ) {
			int index = 2;
			String rep = StandardFileChooser.getSaveName(FileFilters.Report, files.getFile(index));;
			if( !rep.equals( "" ) ) { files.setFile( index, rep ); }
		}

		// If the action comes from the oligomer chemistry group, set the
		// current oligomer chemistry.
		else if( ci.command.equals( "DNA" ) || ci.command.equals( "RNA" ) ) {
			oligomerChemistry = ci.command;
		} else
			super.processCommand(ci);
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField list = FileField.createDisabled( "Oligomer List" ).inputFile(FileFilters.List);
		FileField report = FileField.createEnabled( "Report File" ).outputFile(FileFilters.Report);
		FilePanel files = new FilePanel( this, list, report );
		files.setPanelWidth( 600 );
		files.makePanel();

		// Create the oligo chemistry panel.
		RadioButtonPanel chem = RadioButtonPanel.makeHorizontal(
			"Oligomer Chemistry", this, "DNA", "RNA" );
		int height = chem.getPreferredSize().height;
		chem.setPreferredSize( new Dimension( 200, height ) );
		if( oligomerChemistry != null ) {
			chem.setSelectedButton( oligomerChemistry );
		}

		// Add the components in their proper places.
		setGrid( 2, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, files );
		setGrid( 1, 1 );
		setAnchorCenter();
		setInsets( 0, 10, 0, 0 );
		placeComponent( 0, 1, chem );
		makeStartButton( 1, 1 );
	}

	@Override
	protected boolean runMainCalculation() {

		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );
		RadioButtonPanel chem = (RadioButtonPanel)getInputControl( 2 );
		files.saveRecent();

		// Get the data from the file input panel.
		// If an error occurred while retrieving data, return.
		String list = files.getFile( 1 );
		String report = files.getFile( 2 );
		if( files.isError() ) { return false; }

		// Determine the nucleic acid chemistry type.
		isRNA = chem.getSelectedName().equals( "RNA" );

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress(false);
		String result =
			backend.runOligoScreen( list, report, isRNA );
		displayCalcError( result );
		if( !result.equals( "" ) ) { return false; }

		// Show a message saying that the calculation has finished.
		Dialogs.showMessage( "OligoScreen Complete." );
		return true;
	}

	@Override
	protected MergeMenu[] createCustomMenus() {
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();
		return new ConstraintsMenu[]{ temperature };
	}
}
