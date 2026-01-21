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
import ur_rna.RNAstructureUI.ui.NumberField.DoubleField;
import ur_rna.RNAstructureUI.ui.NumberField.FloatField;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.*;
import ur_rna.Utilities.swing.MergeMenu;

/**
 * A class responsible for initializing and running the AllSub module, which
 * calculates the free energy of a varied group of suboptimal structures.
 *
 * @author Jessica S. Reuter
 */
public class FoldSuboptimalWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public FoldSuboptimalWindow( String acid ) {
		super( acid, "Generate All Suboptimal " + acid + " Structures" );
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

		// If the action comes from the "Sequence File" button, get a sequence
		// file, initialize a data structure, and create a default output file
		// name.
		else if( ci.command.equals( "Sequence File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String file = StandardFileChooser.getOpenName(FileFilters.Sequence);
			if( file == null ) { return; }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildAllSubDataStructure( file, isRNA );
			if( !verifyBackendResult(result, "File: %s\nRNA: %s", file, isRNA) ) { return; }

			// Reset the percent energy difference text field based on the
			// sequence length.
			float size1 = backend.getSuboptimalPercentDiff();
			((FloatField)params.getField( 1 ))
				.resetField( size1, 0, Float.MAX_VALUE );

			// Reset the absolute energy difference text field based on the
			// sequence length.
			double size2 = backend.getSuboptimalAbsoluteDiff();
			((DoubleField)params.getField( 2 ))
				.resetField( size2, 0, Double.MAX_VALUE );

			// Set the sequence file name and the default output file name in
			// the input panel.
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
		FileField seq = FileField.createDisabled( "Sequence File" ).inputFile(FileFilters.Sequence);
		FileField ct = FileField.createEnabled( "CT File" ).outputFile(FileFilters.CT);
		FilePanel files = new FilePanel( this, seq, ct );
		files.setPanelWidth( 600 );
		files.makePanel();

		// Create the suboptimal structures parameter panel.
		FloatField percent =
			new FloatField( "Max % Energy Difference", 0, 0 );
		DoubleField absolute =
			new DoubleField( "Max Absolute Energy Difference", 0, 0 );
		FieldPanel params = new FieldPanel( percent, absolute );
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
		Float percent = ((FloatField)params.getField( 1 )).getValue();
		Double abs = ((DoubleField)params.getField( 2 )).getValue();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress(true, true);
		String result = backend.runAllSub( ct, percent, abs );
		displayCalcError( result );
		if( !result.equals( "" ) ) { return false; }

		// If the user wants to draw structures, draw them.
		if (promptDraw(null, UserOptions.DRAW_PROMPT)) {  //always prompt because it takes a long time.
			showDrawProgress("Drawing Structures (this may take a long time)...");
			DrawingWindow drawing = new DrawingWindow( ct );
			if(drawing.isError())
				return false;
			drawing.showWindow();
		}
		return true;
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

		// Return the array of variable menus.
		return new ConstraintsMenu[]{ temperature, forced };
	}
}
