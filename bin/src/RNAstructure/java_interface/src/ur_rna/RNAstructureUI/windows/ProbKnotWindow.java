/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.ui.FieldPanel;
import ur_rna.RNAstructureUI.ui.FieldPanel.FilePanel;
import ur_rna.RNAstructureUI.ui.FileField;
import ur_rna.RNAstructureUI.ui.NumberField.IntegerField;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.*;
import ur_rna.Utilities.swing.MergeMenu;

/**
 * A class responsible for initializing and running the ProbKnot module, which
 * predicts pseudoknots in a strand of nucleic acids.
 *
 * @author Jessica S. Reuter
 */
public class ProbKnotWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public ProbKnotWindow( String acid ) {
		super( acid, "Identify Pseudoknots" );
	}

	@Override
	protected void processCommand(final CommandInfo ci) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "CT File" button, try to select a CT
		// file, and if one was selected set its name.
		if( ci.command.equals( "CT File" ) ) {
			int index = 2;
			String file = StandardFileChooser.getSaveName(FileFilters.CT, files.getFile(index));
			if( file != null ) { files.setFile( index, file ); }
		}

		// If the action comes from the "Partition Function Save File" button,
		// get a PFS file, initialize a data structure, and create a default
		// output file name.
		else if( ci.command.equals( "Partition Function Save File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String file = StandardFileChooser.getOpenName(FileFilters.PartitionSav);
			if( file == null ) { return; }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildProbKnotDataStructure( file );
			if( !verifyBackendResult(result, "File: %s", file) ) { return; }

			// Set the partition function save file name and the default
			// output file name in the input panel.
			// Then, enable the menus.
			String defaultOut = getOutputFile( file, "ct" );
			files.setFile( 1, file );
			files.setFile( 2, defaultOut );
			getCustomMenus().enableMenus();
		} else
			super.processCommand(ci);
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField pfs =
			FileField.createDisabled( "Partition Function Save File" ).inputFile(FileFilters.PartitionSav);
		FileField ct = FileField.createEnabled( "CT File" ).outputFile(FileFilters.CT);
		FilePanel files = new FilePanel( this, pfs, ct );
		files.setPanelWidth( 600 );
		files.makePanel();

		// Create the parameter panel.
		IntegerField iterations = new IntegerField( "Iterations", 1, 1 );
		IntegerField helix = new IntegerField( "Minimum Helix Length", 3, 1 );
		FieldPanel params = new FieldPanel( iterations, helix );
		params.setPanelWidth( 350 );
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
		Integer iterations = ((IntegerField)params.getField( 1 )).getValue();
		Integer minHelix = ((IntegerField)params.getField( 2 )).getValue();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress(false);
		String result =
			backend.runPseudoknotPrediction( ct, iterations, minHelix );
		displayCalcError( result );
		if( !result.equals( "" ) ) { return false; }

		// If the user wants to draw structures, draw them.
		if (promptDraw()) {
			DrawingWindow drawing = new DrawingWindow( ct );
			if(drawing.isError()) return false;
			drawing.showWindow();
		}
		return true;
	}

	@Override
	protected MergeMenu[] createCustomMenus() { return null; }
}
