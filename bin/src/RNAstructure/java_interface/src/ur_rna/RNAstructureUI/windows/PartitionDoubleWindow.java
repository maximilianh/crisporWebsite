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
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.FileFilters;
import ur_rna.Utilities.swing.MergeMenu;

/**
 * A class responsible for initializing and running the bipartition module,
 * which calculates the partition function for two strands of nucleic acids.
 *
 * @author Jessica S. Reuter
 */
public class PartitionDoubleWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public PartitionDoubleWindow( String acid ) {
		super( acid, acid + " Bimolecular Partition Function" );
	}

	@Override
	protected void processCommand(final CommandInfo ci) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "Save File" button, try to select a
		// partition function save file, and if one was selected set its name.
		if( ci.command.equals( "Save File" ) ) {
			int index = 3;
			String file = StandardFileChooser.getSaveName(FileFilters.PartitionSav, files.getFile(index));;
			if( file != null ) { files.setFile( index, file ); }
		}

		// If the action comes from the "Sequence File 1" button, try to
		// select a sequence file, and if one was selected set its name.
		else if( ci.command.equals( "Sequence File 1" ) ) {
			String file = StandardFileChooser.getOpenName(FileFilters.Sequence);
			if( file != null ) { files.setFile( 1, file ); }
		}

		// If the action comes from the "Sequence File 2" button, get a second
		// sequence file, initialize a data structure, and create a default
		// output file name.
		else if( ci.command.equals( "Sequence File 2" ) ) {

			// Try to select the second sequence file, and if one was selected
			// set its name. If a sequence file wasn't selected, return.
			String file = StandardFileChooser.getOpenName(FileFilters.Sequence);
			if( file != null ) { files.setFile( 2, file ); }
			else { return; }

			// Get the two input sequence files from the panel.
			String seq1 = files.getFile( 1 );
			String seq2 = files.getFile( 2 );

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildBipartitionDataStructure( seq1, seq2, isRNA );
			if( !verifyBackendResult(result, "Seq1: %s\nSeq2: %s\nRNA: %s", seq1, seq2, isRNA) ) { return; }

			// Build the default output file name and set it in the panel.
			String defaultOut = combineFileNames( seq1, seq2, "pfs" );
			files.setFile( 3, defaultOut );

			// Enable the menus.
			getCustomMenus().enableMenus();
		} else
			super.processCommand(ci);
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField seq1 = FileField.createDisabled( "Sequence File 1" ).inputFile(FileFilters.Sequence);
		FileField seq2 = FileField.createDisabled( "Sequence File 2" ).inputFile(FileFilters.Sequence);
		FileField pfs = FileField.createEnabled( "Save File" ).outputFile(FileFilters.PartitionSav);
		FilePanel files = new FilePanel( this, seq1, seq2, pfs );
		files.setPanelWidth( 600 );
		files.makePanel();

		// Put the components in their proper places.
		setGrid( 2, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, files );
		makeStartButton( 0, 1 );
	}

	@Override
	protected boolean runMainCalculation() {
		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );
		files.saveRecent();

		// Get the data from the file input panel.
		// If an error occurred while retrieving data, return.
		files.getFile( 1 );
		files.getFile( 2 );
		String pfs = files.getFile( 3 );
		if( files.isError() ) { return false; }

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress();
		String result = backend.runBipartition( pfs );
		displayCalcError( result );
		if( !result.equals( "" ) ) { return false; }

		// Show the dot plot.
		if (promptDraw()) {
			showDrawProgress();
			DrawingWindow drawing = new DrawingWindow(pfs);
			if (drawing.isError()) return false;
			drawing.showWindow();
		}
		return true;
	}

	@Override
	protected MergeMenu[] createCustomMenus() {
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();
		return new ConstraintsMenu[]{ temperature };
	}
}
