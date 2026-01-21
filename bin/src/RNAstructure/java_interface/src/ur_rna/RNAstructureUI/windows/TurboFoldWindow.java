/**
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.ui.BorderBuilder;
import ur_rna.RNAstructureUI.ui.FieldPanel;
import ur_rna.RNAstructureUI.ui.NumberField.DoubleField;
import ur_rna.RNAstructureUI.ui.NumberField.IntegerField;
import ur_rna.RNAstructureUI.ui.RadioButtonPanel;

import javax.swing.*;
import java.awt.*;

/**
 * A class responsible for initializing and running the TurboFold module,
 * which calculates common secondary structures of multiple sequences using
 * the partition function and probabalistic alignments. This module is only
 * used with RNA.
 *
 * @author Jessica S. Reuter
 */
public class TurboFoldWindow
	extends MultiWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * The layout that holds the dynamic options applicable to a mode.
	 */
	private CardLayout cards;

	/**
	 * An array that holds the names of the TurboFold modes.
	 */
	private static final String[] modeNames = {
		"Maximum Expected Accuracy",
		"ProbKnot/TurboKnot",
		"Probability Threshold"
	};

	/**
	 * The panel that holds the dynamic options applicable to a mode.
	 */
	private JPanel specificPanel;

	/**
	 * Constructor.
	 * <br><br>
	 * Note that since this module is unique to RNA, the nucleic acid type is
	 * hard-coded as RNA.
	 */
	public TurboFoldWindow() { super( "RNA", "TurboFold" ); }

	@Override
	protected void activate() { backend.activateTurboFold(); }

	@Override
	protected void addBackendTuple( String seq, String ct ) {
		backend.addTurboFoldTuple( seq, ct );
	}

	@Override
	protected void buildSpecificOptionsPanel() {

		// Create a border builder for padding input controls.
		BorderBuilder borderBuilder = new BorderBuilder();

		// Create the mode selection panel.
		RadioButtonPanel modePanel =
			RadioButtonPanel.makeVertical( "Mode", this, modeNames );

		// Create the general TurboFold options panel.
		DoubleField turboGamma = new DoubleField( "Gamma", 0.3 );
		IntegerField turboIterations = new IntegerField( "Iterations", 3, 1 );
		FieldPanel general = new FieldPanel( turboGamma, turboIterations );
		general.setPanelWidth( columnWidth );
		general.makePanel();
		borderBuilder.makeTitledBorder( "General Options", general );

		// Create the maximum expected accuracy mode options box.
		DoubleField score =
			new DoubleField( "Max % Score Difference", 50.0, 0.0, 99.0 );
		IntegerField structures =
			new IntegerField( "Max Number of Structures", 1000, 1 );
		IntegerField window = new IntegerField( "Window Size", 5, 0 );
		DoubleField meaGamma = new DoubleField( "MEA Gamma", 1 );
		FieldPanel mea =
			new FieldPanel( score, structures, window, meaGamma );
		mea.setPanelWidth( columnWidth );
		mea.makePanel();
		Box meaBox = Box.createVerticalBox();
		meaBox.add( mea );
		meaBox.add( Box.createVerticalStrut( columnWidth ) );

		// Create the pseudoknots mode options panel.
		IntegerField pkIterations =
			new IntegerField( "Pseudoknot Iterations", 1, 1 );
		IntegerField helix = new IntegerField( "Minimum Helix Length", 3, 1 );
		FieldPanel knot = new FieldPanel( pkIterations, helix );
		knot.setPanelWidth( columnWidth );
		knot.makePanel();
		Box knotBox = Box.createVerticalBox();
		knotBox.add( knot );
		knotBox.add( Box.createVerticalStrut( columnWidth ) );

		// Create the threshold mode options panel.
		DoubleField thresh = new DoubleField( "Threshold", 0, 0, 99 );
		FieldPanel threshold = new FieldPanel( thresh );
		threshold.setPanelWidth( columnWidth );
		threshold.makePanel();
		Box thresholdBox = Box.createVerticalBox();
		thresholdBox.add( threshold );
		thresholdBox.add( Box.createVerticalStrut( columnWidth ) );

		// Create a panel that holds the mode-specific options.
		cards = new CardLayout();
		specificPanel = new JPanel( cards );
		specificPanel.add( meaBox, modeNames[0] );
		specificPanel.add( knotBox, modeNames[1] );
		specificPanel.add( thresholdBox, modeNames[2] );
		cards.show( specificPanel, modeNames[0] );

		// Add the panels to the options box.
		options.add( general );
		options.add( modePanel );
		options.add( specificPanel );
		borderBuilder.makeLeftBorder( 10, options );
	}

	@Override
	protected void deleteBackendTuple( int index ) {
		backend.deleteTurboFoldTuple( index );
	}

	@Override
	protected void processCommand(final CommandInfo ci) {
		// If the command is one of the mode names, show that mode's options.
		for( String name: modeNames ) {
			if( ci.command.equals( name ) ) {
				cards.show( specificPanel, name );
				return;
			}
		}
		// Call the superclass module action method.
		super.processCommand(ci);
	}

	@Override
	protected int getNumSequences() {
		return backend.getNumTurboFoldSequences();
	}

	@Override
	protected String getSequenceSetAsString() {
		return backend.getTurboFoldSequenceSetData();
	}

	@Override
	protected boolean runMainCalculationMultipleSequences() {

		// Show a progress bar.
		showProgress(false);

		// Get the specific input options panel for reference.
		JPanel opts = (JPanel)getInputControl( options, 3 );

		// Get the mode from the window, and initialize a string to hold the
		// result of a calculation.
		RadioButtonPanel modePanel =
			(RadioButtonPanel)getInputControl( options, 2 );
		String mode = modePanel.getSelectedName();
		String result = "";

		// Get the general TurboFold parameters.
		FieldPanel general = (FieldPanel)getInputControl( options, 1 );
		double turboGamma = ((DoubleField)general.getField( 1 )).getValue();
		int turboIterations = ((IntegerField)general.getField( 2 )).getValue();

		// If selected mode is maximum expected accuracy, run TurboFold in
		// maximum expected accuracy mode.
		if( mode.equals( "Maximum Expected Accuracy" ) ) {

			// Get the data values for maximum expected accuracy mode.
			Box optsMode = (Box)getInputControl( opts, 1 );
			FieldPanel mea = (FieldPanel)getInputControl( optsMode, 1 );
			double percent = ((DoubleField)mea.getField( 1 )).getValue();
			int structures = ((IntegerField)mea.getField( 2 )).getValue();
			int window = ((IntegerField)mea.getField( 3 )).getValue();
			double meaGamma = ((DoubleField)mea.getField( 4 )).getValue();

			// If all validation of input went fine, run TurboFold in maximum
			// expected accuracy mode.
			result = backend.runTurboFoldMaximumExpectedAccuracy(
				turboGamma, turboIterations, percent, structures, window,
				meaGamma );
		}

		// If selected mode is pseudoknots, run TurboFold in pseudoknots mode.
		else if( mode.equals( "ProbKnot/TurboKnot" ) ) {

			// Get the data values for pseudoknots mode.
			Box optsMode = (Box)getInputControl( opts, 2 );
			FieldPanel knot = (FieldPanel)getInputControl( optsMode, 1 );
			int pkIterations = ((IntegerField)knot.getField( 1 )).getValue();
			int helix = ((IntegerField)knot.getField( 2 )).getValue();

			// If all validation of input went fine, run TurboFold in
			// pseudoknots mode.
			result = backend.runTurboFoldPseudoknot(
				turboGamma, turboIterations, pkIterations, helix );
		}

		// If selected mode is probability thresholds, run TurboFold in
		// threshold mode.
		else if( mode.equals( "Probability Threshold" ) ) {

			// Get the data values for threshold mode.
			Box optsMode = (Box)getInputControl( opts, 3 );
			FieldPanel thresh = (FieldPanel)getInputControl( optsMode, 1 );
			Double threshold = ((DoubleField)thresh.getField( 1 )).getValue();

			// If all validation of input went fine, run TurboFold in
			// threshold mode.
			result = backend.runTurboFoldThreshold(
				turboGamma, turboIterations, threshold );
		}

		// Finish the module by closing the window and checking for errors.
		// If an error occurred during calculation, return.
		displayCalcError( result );
		if( !result.equals( "" ) ) { return false; }

		// Ask if the user wants to draw structures and dot plots.
		// If not, return.
		if( !promptDraw("Do you want to draw\nstructures and dot plots?", null)) { return true; }

		showDrawProgress();
		// If the user wants to draw, build the structure and plot dialogs.
		// Show them if they were built correctly.
		for( int i = 1; i <= getNumSequences(); i++ ) {
			String ct = backend.getTurboFoldCT( i );
			imageDialogs.add( new DrawingWindow( ct ) );
			String pfs = backend.getTurboFoldSaveFile( i );
			imageDialogs.add( new DrawingWindow( pfs ) );
		}
		showImageDialogs();
		return true;
	}
}
