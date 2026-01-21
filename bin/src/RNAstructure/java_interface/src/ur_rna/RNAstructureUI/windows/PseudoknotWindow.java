/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.menus.ConstraintsMenu;
import ur_rna.RNAstructureUI.ui.FieldPanel.FilePanel;
import ur_rna.RNAstructureUI.ui.FileField;
import ur_rna.RNAstructureUI.ui.RadioButtonPanel;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.*;
import ur_rna.Utilities.swing.MergeMenu;

import java.awt.*;

/**
 * A class responsible for initializing and running the RemovePseudoknots
 * module, which removes pseudoknots from a nucleic acid sequence. This module
 * is only used with RNA.
 *
 * @author Jessica S. Reuter
 */
public class PseudoknotWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public PseudoknotWindow( String acid ) {
		super( acid, "Break " + acid + " Pseudoknots" );
	}

	@Override
	protected void processCommand(final CommandInfo ci) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "Input CT File" button, get a CT file,
		// initialize a data structure, and create a default output file name.
		if( ci.command.equals( "Input CT File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String file = StandardFileChooser.getOpenName(FileFilters.CT);
			if( file == null ) { return; }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildRemovePseudoknotsDataStructure( file, isRNA );
			if( !verifyBackendResult(result, "File: %s\nRNA: %s", file, isRNA) ) { return; }


			// Set the CT file name and the default output file name in the
			// input panel.
			// Then, enable the menus.
			String root = file.substring( 0, file.lastIndexOf( "." ) );
			String defaultOut = root + "_no_pseudo.ct";
			files.setFile( 1, file );
			files.setFile( 2, defaultOut );
			getCustomMenus().enableMenus();
		}

		// If the action comes from the "Output CT File" button, try to select
		// a CT file, and if one was selected set its name.
		else if( ci.command.equals( "Output CT File" ) ) {
			int index = 2;
			String file = StandardFileChooser.getSaveName(FileFilters.CT, files.getFile(index));
			if( file != null ) { files.setFile( index, file ); }
		} else
			super.processCommand(ci);
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField ctIn = FileField.createDisabled( "Input CT File" ).inputFile(FileFilters.CT);
		FileField ctOut = FileField.createEnabled( "Output CT File" ).outputFile(FileFilters.CT);
		FilePanel files = new FilePanel( this, ctIn, ctOut );
		files.setPanelWidth( 600 );
		files.makePanel();

		// Create the pseudoknot options panel.
		String[] buttons = {
			"Minimize Folding Free Energy",
			"Maximize Number of Pairs"
		};
		RadioButtonPanel options = RadioButtonPanel.makeVertical(
			"Pseudoknot Breakage Options", buttons );
		int height = options.getPreferredSize().height;
		options.setPreferredSize( new Dimension( 250, height ) );

		// Put the components in their proper places.
		setGrid( 2, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, files );
		setGrid( 1, 1 );
		placeComponent( 0, 1, options );
		makeStartButton( 1, 1 );
	}

	@Override
	protected boolean runMainCalculation() {

		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );
		RadioButtonPanel options = (RadioButtonPanel)getInputControl( 2 );
		files.saveRecent();

		// Get the data from the window.
		// If an error occurred while retrieving data, return.
		files.getFile( 1 );
		String ct = files.getFile( 2 );
		if( files.isError() ) { return false; }

		// Determine if free energy should be minimized.
		boolean minimize = ( options.getSelectedIndex() == 1 );

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress(false);
		String result =
			backend.runPseudoknotRemoval( ct, minimize );
		displayCalcError( result );
		if( !result.equals( "" ) ) { return false; }

		// If the user wants to draw structures, draw them.
		return drawStructures(ct);
	}

	@Override
	protected MergeMenu[] createCustomMenus() {
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();
		return new ConstraintsMenu[]{ temperature };
	}
}
