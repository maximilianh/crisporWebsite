/**
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.ui.*;
import ur_rna.RNAstructureUI.ui.FieldPanel.FilePanel;
import ur_rna.RNAstructureUI.ui.NumberField.DoubleField;
import ur_rna.RNAstructureUI.ui.NumberField.IntegerField;
import ur_rna.RNAstructureUI.utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * A class responsible for initializing and running the Multilign module,
 * which calculates common secondary structures of multiple sequences using
 * multiple iterations of Dynalign.
 *
 * @author Jessica S. Reuter
 */
public class MultilignWindow extends MultiWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public MultilignWindow( String acid ) { super( acid, "Multilign" ); }

	@Override
	protected void activate() { backend.activateMultilign(); }

	@Override
	protected void addBackendTuple( String seq, String ct ) {
		backend.addMultilignTuple( seq, ct );
	}

	@Override
	protected void buildSpecificOptionsPanel() {

		// Create a border builder for padding input controls.
		BorderBuilder borderBuilder = new BorderBuilder();

		// Create a slightly shrunken width for objects inside the tabbed pane.
		int tabWidth = columnWidth - 10;

		// Create the alignment file input panel.
		FileField align = FileField.createEnabled( "Alignment File" );
		FilePanel alignPanel = new FilePanel( this, align );
		alignPanel.setPanelWidth( columnWidth );
		alignPanel.makePanel();

		// Create the tabbed pane for parameters.
		JTabbedPane tabbedPane = new JTabbedPane();
		tabbedPane.setPreferredSize( new Dimension( columnWidth, 375 ) );
		borderBuilder.makeEqualBorder( 10, tabbedPane );

		// Create the Dynalign parameters tab and add it to the tabbed pane.
		// This is placed inside a block for organizational purposes.
		{
			// Create the suboptimal structure parameters panel.
			IntegerField energy =
				new IntegerField( "Max % Energy Difference", 20, 0 );
			IntegerField structures =
				new IntegerField( "Max Number of Structures", 20, 1 );
			IntegerField windowStruct =
				new IntegerField( "Structure Window Size", 2, 0 );
			IntegerField windowAlign =
				new IntegerField( "Alignment Window Size", 1, 0 );
			DoubleField gap =
				new DoubleField( "Gap Penalty", 0.4 );
			FieldPanel suboptimal = new FieldPanel(
				energy, structures, windowStruct, windowAlign, gap );
			suboptimal.setPanelWidth( tabWidth );
			suboptimal.makePanel();

			// Create the base pairs insertion check box.
			HTMLCheckBox bpBox =
				HTMLCheckBox.createSelectedBox( "Single BP Inserts Allowed" );
			borderBuilder.makeLeftBorder( 10, bpBox );

			// Create the save/alignment files generation check box.
			HTMLCheckBox saveBox = HTMLCheckBox.createEmptyBox(
				"Generate Save Files and Alignment Files" );
			borderBuilder.makeLeftBorder( 10, saveBox );

			// Create a smaller panel to hold the two text boxes.
			JPanel boxPanel = new JPanel( new GridLayout( 0, 1 ) );
			boxPanel.setPreferredSize( new Dimension( tabWidth, 30 ) );
			boxPanel.add( bpBox );
			boxPanel.add( saveBox );

			// Put all the Dynalign parameters together into one box, and add
			// add that box to the tabbed pane.
			Box dynalign = Box.createVerticalBox();
			dynalign.add( suboptimal );
			dynalign.add( boxPanel );
			dynalign.add( Box.createVerticalGlue() );
			tabbedPane.addTab( "Dynalign Parameters", dynalign );
		}

		// Create the Multilign specific options panel and add it to the
		// tabbed pane.
		IntegerField iterations = new IntegerField( "Iterations", 2, 1 );
		IntegerField maxPairs = new IntegerField( "MaxPairs", 0, 0 );
		DoubleField dsv = new DoubleField( "maxdsvchange", 1, 0, 99 );
		FieldPanel multilign = new FieldPanel( iterations, maxPairs, dsv );
		multilign.setPanelWidth( tabWidth );
		multilign.makePanel();
		Box multilignBox = Box.createVerticalBox();
		multilignBox.add( multilign );
		multilignBox.add( Box.createVerticalGlue() );
		tabbedPane.addTab( "Multilign Parameters", multilignBox );

		// Add the alignment panel and tabbed pane to the options box.
		options.add( alignPanel );
		options.add( tabbedPane );
	}

	@Override
	protected void deleteBackendTuple( int index ) {
		backend.deleteMultilignTuple( index );
	}

	@Override
	protected void processCommand(final CommandInfo ci) {
		// Get any input controls from the window that handle actions.
		Box opts = (Box)getInputControl( 3 );
		FilePanel alignPanel = (FilePanel)getInputControl( opts, 1 );
		JTabbedPane tabs = (JTabbedPane)getInputControl( opts, 2 );

		Box tabBox = (Box)getInputControl( tabs, 2 );
		FieldPanel multilignPanel = (FieldPanel)getInputControl( tabBox, 1 );

		// If the command comes from the "ADD" button or the "Delete Sequence"
		// button, set the new max pairs value. Then, put the default file
		// name in for the alignment file, if it hasn't been done already.
		if( ci.command.equals( "ADD -->" ) ) {
			// Let the MultiWindow parent process the event first.
			super.processCommand(ci);
			// Set the max pairs.
			int maxPairs = backend.getMultilignMaxPairs();
			((IntegerField)multilignPanel.getField( 2 )).resetField( maxPairs );

			// Check to see if the alignment file should be set.
			boolean setAlignFile =
				( alignPanel.getFile( 1 ).equals( "" ) ) &&
				( getNumSequences() == 2 );

			// Set the alignment file if it should be set.
			if( setAlignFile ) {
				String align =
					System.getProperty( "user.home" ) +
					System.getProperty( "file.separator" ) + "multi.ali";
				alignPanel.setFile( 1, align );
			}
		}

		// If the action comes from the "Alignment File" button, try to select
		// an alignment file, and if one was selected set its name.
		else if( ci.command.equals( "Alignment File" ) ) {
			int index = 1;
			String file = StandardFileChooser.getSaveName(FileFilters.Alignment, alignPanel.getFile(index));;
			if( file != null ) { alignPanel.setFile( index, file ); }
		}

		// If the command comes from the "Delete Sequence" button, set the new
		// max pairs value.
		else if( ci.command.equals( "Delete Sequence" ) ) {
			// Let the MultiWindow parent process the event first.
			super.processCommand(ci);
			int maxPairs = backend.getMultilignMaxPairs();
			((IntegerField)multilignPanel.getField( 2 )).resetField( maxPairs );
		} else
			// Let the MultiWindow parent process the event also.
			super.processCommand(ci);


	}

	@Override
	protected int getNumSequences() {
		return backend.getNumMultilignSequences();
	}

	@Override
	protected String getSequenceSetAsString() {
		return backend.getMultilignSequenceSetData();
	}

	@Override
	protected boolean runMainCalculationMultipleSequences() {
		// Show a progress bar.
		showProgress();

		// Get the specific input options panels for reference.
		Box opts = (Box)getInputControl( 3 );
		JTabbedPane tabs = (JTabbedPane)getInputControl( opts, 2 );

		// Get the alignment file string. If an error occurred getting that
		// string, show an error and do not continue.
		FilePanel align = (FilePanel)getInputControl( opts, 1 );
		align.saveRecent();
		String alignString = align.getFile( 1 );
		if( align.isError() ) { return false; }

		// Get the suboptimal structure parameters.
		Box dynalignTab = (Box)getInputControl( tabs, 1 );
		FieldPanel suboptimal =
			(FieldPanel)getInputControl( dynalignTab, 1 );
		Integer percent = ((IntegerField)suboptimal.getField( 1 )).getValue();
		Integer structs = ((IntegerField)suboptimal.getField( 2 )).getValue();
		Integer bpWin = ((IntegerField)suboptimal.getField( 3 )).getValue();
		Integer alignWin = ((IntegerField)suboptimal.getField( 4 )).getValue();
		Double gap = ((DoubleField)suboptimal.getField( 5 )).getValue();

		// Set whether base pair inserts are allowed and whether save files
		// and alignment files should be generated, based on if the proper
		// boxes are checked.
		JPanel boxes = (JPanel)getInputControl( dynalignTab, 2 );
		boolean insert =
			((HTMLCheckBox)getInputControl( boxes, 1 )).isSelected();
		boolean saveFiles =
			((HTMLCheckBox)getInputControl( boxes, 2 )).isSelected();

		// Get data from the Multilign specific panel.
		Box multilignTab = (Box)getInputControl( tabs, 2 );
		FieldPanel multilign =
			(FieldPanel)getInputControl( multilignTab, 1 );
		Integer cycles = ((IntegerField)multilign.getField( 1 )).getValue();
		Integer maxPairs = ((IntegerField)multilign.getField( 2 )).getValue();
		Double maxDsv = ((DoubleField)multilign.getField( 3 )).getValue();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress();
		String result =
			backend.runMultilign(
				percent, structs, bpWin, alignWin, gap, insert, maxDsv,
				maxPairs, cycles, alignString, saveFiles, isRNA );
		displayCalcError( result );
		if( !result.equals( "" ) ) { return false; }

		// Ask if the user wants to draw structures, and return if not.
		if (!promptDraw()) return true;
		showDrawProgress();
		// If the user wants to draw structures, build the structure dialogs.
		// Show them if they were built correctly.
		for( int i = 1; i <= getNumSequences(); i++ ) {
			String ct = backend.getMultilignCT( i );
			imageDialogs.add( new DrawingWindow( ct ) );
		}
		showImageDialogs();
		return true;
	}
}
