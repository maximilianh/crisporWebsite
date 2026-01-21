/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.menus.ConstraintsMenu;
import ur_rna.RNAstructureUI.ui.*;
import ur_rna.RNAstructureUI.ui.FieldPanel.FilePanel;
import ur_rna.RNAstructureUI.ui.NumberField.DoubleField;
import ur_rna.RNAstructureUI.ui.NumberField.FloatField;
import ur_rna.RNAstructureUI.ui.NumberField.IntegerField;
import ur_rna.RNAstructureUI.utilities.FileFilters;
import ur_rna.Utilities.swing.MergeMenu;

import javax.swing.text.JTextComponent;

/**
 * A class responsible for initializing and running the AccessFold module, which
 * folds two strands of nucleic acids into their likely hybrid structure.
 *
 * @author Richard M. Watson
 */
public class AccessFoldWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20160714;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public AccessFoldWindow( String acid ) {
		super( acid, acid + " AccessFold" );
	}

	@Override
	protected void processCommand(final CommandInfo ci) {
		JTextComponent field;
		String file;
		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "CT File" button, try to select a CT
		// file, and if one was selected set its name.
		switch (ci.command) {
			case "CT File":
				field = files.getField(ci.command);
				file = StandardFileChooser.getSaveName(FileFilters.CT, field.getText());
				if (file != null) { field.setText(file); }
				break;

			// If the action comes from the "Sequence File 1" button, try to
			// select a sequence file, and if one was selected set its name.
			case "Sequence File 1":
				field = files.getField(1);
				file = StandardFileChooser.getOpenName(FileFilters.Sequence);
				if (file != null) { field.setText(file); }
				break;

			// If the action comes from the "Sequence File 2" button, get a second
			// sequence file, initialize a data structure, and create a default
			// output file name.
			// Note that all this can only be done if the first sequence file has
			// already been selected.
			case "Sequence File 2": {
				// If the first sequence file hasn't been selected yet, show an
				// error and return.
				if (files.getFile(1).equals("")) {
					String message = "Sequence File 1 hasn't been selected yet.";
					Dialogs.showError(message);
					return;
				}

				// Try to select the second sequence file, and if one was selected
				// set its name. If a sequence file wasn't selected, return.
				field = files.getField(2);
				file = StandardFileChooser.getOpenName(FileFilters.Sequence);
				if (file != null) { field.setText(file); } else { return; }

				// Get the two input sequence files from the panel.
				String seq1 = files.getFile(1);
				String seq2 = files.getFile(2);

				// Create a data structure.
				// If an error occurred creating the data structure, show an error
				// and return.
				String result = backend.buildAccessFoldDataStructure(seq1, seq2, isRNA);
				if (!verifyBackendResult(result, "Seq1: %s\nSeq2: %s\nRNA: %s", seq1, seq2, isRNA)) { return; }

				// Build the default output file name and set it in the panel.
				String defaultOut = combineFileNames(seq1, seq2, "ct");
				files.setFile(3, defaultOut);

				// Enable the menus.
				getCustomMenus().enableMenus();
				break;
			}
			default:
				super.processCommand(ci);
				break;
		}
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField seq1 = FileField.createDisabled( "Sequence File 1" ).inputFile(FileFilters.Sequence);
		FileField seq2 = FileField.createDisabled( "Sequence File 2" ).inputFile(FileFilters.Sequence);
		FileField ct = FileField.createEnabled( "CT File" ).outputFile(FileFilters.CT);
		FilePanel files = new FilePanel( this, seq1, seq2, ct );
		files.setPanelWidth( 600 );
		files.makePanel();

		// Create the box that allows the user to create a save file.
		HTMLCheckBox box =
			HTMLCheckBox.createEmptyBox( "Generate Save File" );

		// Create the parameter panel.
		DoubleField gamma        = new DoubleField( "Accessibility factor (gamma)", 0.4, 0 );
		FloatField energy       = new FloatField( "Max % Energy Difference", 50, 0 );
		IntegerField structures = new IntegerField( "Max Number of Structures", 20, 1 );
		IntegerField window     = new IntegerField( "Window Size", 0, 0 );
		FieldPanel params = new FieldPanel( gamma, energy, structures, window );
		params.setPanelWidth( 300 );
		params.makePanel();

		// Add the components in their proper places.
		setGrid( 2, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, files );
		setFillCenter(); 
		placeComponent( 0, 1, box );
		setGrid( 1, 1 );
		setFillHorizontal(); 
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
		files.getFile( 2 );
		String ct = files.getFile( 3 );
		if( files.isError() ) { return false; }

		// Get the data from the parameters panel.
		Double gamma = ((DoubleField)params.getField( 1 )).getValue();
		Float percent = ((FloatField)params.getField( 2 )).getValue();
		Integer structures = ((IntegerField)params.getField( 3 )).getValue();
		Integer window = ((IntegerField)params.getField( 4 )).getValue();

		// Get whether a save file should be written, and whether
		// intramolecular pairs should be forbidden.
		boolean save = box.isSelected();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress();
		String result = backend.runAccessFold( ct, percent, structures, gamma, window, save);
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

		// Create the menu that determines if unimolecular pairs are allowed.
		// ConstraintsMenu pairs = new ConstraintsMenu( backend );
		// pairs.buildUnimolecularMenu();

		// Create the maximum loop menu.
		ConstraintsMenu loop = new ConstraintsMenu( backend );
		loop.buildMaxLoopMenu();

		// Return the array of variable menus.
		return new ConstraintsMenu[]{ temperature, loop };
	}
}
