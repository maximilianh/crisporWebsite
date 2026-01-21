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
 * A class responsible for initializing and running the Dynalign refold
 * module, which refolds a previous Dynalign calculation from a save file.
 *
 * @author Jessica S. Reuter
 */
public class DynalignRefoldWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 * <br><br>
	 * Note that since the Dynalign save file knows its nucleic acid type, the
	 * exact nucleic acid type this window uses doesn't matter, and the
	 * nucleic acid type is specified as an empty string.
	 */
	public DynalignRefoldWindow() {
		super( "", "Refold From Dynalign Save File" );
	}

	@Override
	protected void processCommand(final CommandInfo ci) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "Alignment File" button, try to select
		// an alignment file, and if one was selected set its name.
		if( ci.command.equals( "Alignment File" ) ) {
			String file = StandardFileChooser.getSaveName(FileFilters.Alignment);
			if( file != null ) { files.setFile( 4, file ); }
		}

		// If the action comes from one of the "CT File" buttons, try to
		// select a CT file, and if one was selected set its name, based on
		// its index.
		else if( ci.command.startsWith( "CT File" ) ) {
			String file = StandardFileChooser.getSaveName(FileFilters.CT);
			if( file != null ) {
				int index = ( ci.command.endsWith( "1" ) ) ? 2 : 3;
				files.setFile( index, file );
			}
		}

		// If the action comes from the "Save File" button, get a Dynalign
		// save file and initialize a data structure.
		else if( ci.command.equals( "Save File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			// Otherwise, set the file in the input panel.
			String file = StandardFileChooser.getOpenName(FileFilters.DynalignSav);
			if( file == null ) { return; }
			else { files.setFile( 1, file ); }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildRefoldDynalignDataStructure( file );
			if( !verifyBackendResult(result, "File: %s", file) ) { return; }

			// Enable the menus.
			getCustomMenus().enableMenus();
		} else
			super.processCommand(ci);
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField save = FileField.createDisabled( "Save File" ).inputFile(FileFilters.DynalignSav);
		FileField ct1 = FileField.createEnabled( "CT File 1" ).outputFile(FileFilters.CT);
		FileField ct2 = FileField.createEnabled( "CT File 2" ).outputFile(FileFilters.CT);
		FileField align = FileField.createEnabled( "Alignment File" ).outputFile(FileFilters.Alignment);
		FilePanel files = new FilePanel( this, save, ct1, ct2, align );
		files.setPanelWidth( 600 );
		files.makePanel();

		IntegerField energy =
			new IntegerField( "Max % Energy Difference", 20, 0 );
		IntegerField structures =
			new IntegerField( "Max Number of Structures", 750, 1 );
		IntegerField windowStruct =
			new IntegerField( "Structure Window Size", 0, 0 );
		IntegerField windowAlign =
			new IntegerField( "Alignment Window Size", 0, 0 );
		FieldPanel params =
			new FieldPanel( energy, structures, windowStruct, windowAlign );
		params.setPanelWidth( 300 );
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
		String ct1 = files.getFile( 2 );
		String ct2 = files.getFile( 3 );
		String align = files.getFile( 4 );
		if( files.isError() ) { return false; }

		// Get the data from the parameters panel.
		Integer percent = ((IntegerField)params.getField( 1 )).getValue();
		Integer structures = ((IntegerField)params.getField( 2 )).getValue();
		Integer windowStr = ((IntegerField)params.getField( 3 )).getValue();
		Integer windowAli = ((IntegerField)params.getField( 4 )).getValue();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress();
		String result =
			backend.runDynalignRefold(
				ct1, ct2, align, percent, structures, windowStr, windowAli );
		displayCalcError( result );
		if( !result.equals( "" ) ) { return false; }

		// If the user wants to draw structures, draw them.
		if (promptDraw()) {
			showDrawProgress();
			DrawingWindow one = new DrawingWindow( ct1 );
			DrawingWindow two = new DrawingWindow( ct2 );
			if(one.isError() || two.isError()) return false;
			one.showWindow();
			two.showWindow();
		}
		return true;
	}

	@Override
	protected MergeMenu[] createCustomMenus() { return null; }
}
