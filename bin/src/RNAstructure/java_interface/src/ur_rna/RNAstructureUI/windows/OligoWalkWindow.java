/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.AppMainFrame;
import ur_rna.RNAstructureUI.menus.ConstraintsMenu;
import ur_rna.RNAstructureUI.menus.MainMenu;
import ur_rna.RNAstructureUI.ui.*;
import ur_rna.RNAstructureUI.ui.FieldPanel.FilePanel;
import ur_rna.RNAstructureUI.ui.NumberField.IntegerField;
import ur_rna.RNAstructureUI.utilities.*;
import ur_rna.Utilities.swing.MergeMenu;

import javax.swing.*;
import javax.swing.JSpinner.DefaultEditor;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.*;
import java.io.File;

/**
 * A class that handles creation of windows for the OligoWalk module. The
 * OligoWalk module creates two linked windows, an input window and an output
 * (results) window. The two windows share a back end calculator. This module
 * is only used with RNA.
 * <br><br>
 * Note also that this class violates the usual 500 line class length coding
 * standard limit. This is because both the input window and the results window
 * contain objects unique to themselves, and the results window, which is more
 * complex than the input window, relies on the same data structure as the
 * input window.
 *
 * @author Jessica S. Reuter
 */
public class OligoWalkWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * A boolean that tells whether this window is an input window or not.
	 */
	private boolean isInputWindow;

	/*
	 * The panel that holds input controls.
	 */
	private JPanel panelInput;

	/**
	 * The box that holds the results.
	 */
	private JPanel panelResults;

	/**
	 * Constructor.
	 * <br><br>
	 * Note that since this module is unique to RNA, the nucleic acid type is
	 * hard-coded as RNA.
	 */
	public OligoWalkWindow() {
		super( "RNA", "RNA OligoWalk" );

		// Set the window as an input window and set the appropriate menus.
		isInputWindow = true;
	}

	/**
	 * Build the results window for this OligoWalk calculation.
	 *
	 * @param file     The report file written for this calculation.
	 * @param target   The length of the entire target.
	 * @param start    The start index where hybridization can occur.
	 * @param stop     The stop index where hybridization can occur.
	 */
	private void buildResultsWindow(
		String file, int target, int start, int stop ) {

		// Set window type to a results window and set its title.
		// Then, replace the input window menu with the results window menu.
		isInputWindow = false;
		setCaption( "OligoWalk Results: " + file );

		getCustomMenus().clear();
		getCustomMenus().add(new OligoWalkResultsMenu());


		// Add a key listener that allows the user to move between oligos.
		addKeyListener( new KeyAdapter() {

			@Override
			public void keyPressed( KeyEvent e ) {
				int keyCode = e.getKeyCode();
				Integer move = 0;
				if( keyCode == KeyEvent.VK_RIGHT ) { move = 1; }
				else if( keyCode == KeyEvent.VK_LEFT ) { move = -1; }
				invokeCommand( "Move Oligo " + Integer.toString( move ), this);
			}
		});

		// Get the pieces of the graph.
		JPanel graphPanel = (JPanel)panelResults.getComponent( 2 );
		ScrollerPane scroller = (ScrollerPane)graphPanel.getComponent( 1 );
		JPanel graphWrapper = (JPanel)scroller.getViewport().getView();
		JLabel seqLabel = (JLabel)graphWrapper.getComponent( 0 );
		Box graph = (Box)graphWrapper.getComponent( 2 );
		Box top = (Box)graph.getComponent( 0 );
		Box bottom = (Box)graph.getComponent( 1 );

		// Add the bars to the graph.
		// This loop starts at 0 and ends at target + 1 so the spacers which
		// frame the left and right sides of the graph can be constructed at
		// the same time as the graph bars.
		int numBarsForEachIndex = 3;
		Dimension barSize = new Dimension( 10, 10 );
		for( int i = 0; i <= target + 1; i++ ) {

			// Create the bars for each index.
			// If the index is in hybridization range, make buttons.
			// Otherwise, make spacers.
			JComponent[] components = new JComponent[numBarsForEachIndex];
			for( int j = 1; j <= numBarsForEachIndex; j++ ) {
				if( ( i >= start ) && ( i <= stop ) ) {
					BarButton button = new BarButton();
					button.setActionCommand(
						"Set Oligo " + Integer.toString( i ) );
					button.addActionListener( this );
					button.setBorder( new LineBorder( Color.BLACK, 1 ) );
					components[j-1] = button;
				} else { components[j-1] = new JLabel(); }
			}

			// Set the sizes for the components.
			for( Component component: components ) {
				component.setPreferredSize( barSize );
				component.setMinimumSize( barSize );
				component.setMaximumSize( barSize );
			}

			// Set the alignments for the components.
			components[0].setAlignmentY( JButton.BOTTOM_ALIGNMENT );
			components[1].setAlignmentY( JButton.TOP_ALIGNMENT );
			components[2].setAlignmentY( JButton.TOP_ALIGNMENT );

			// Create the stack of buttons for the bottom half of the graph.
			Box lowerStack = Box.createVerticalBox();
			lowerStack.setAlignmentY( JButton.TOP_ALIGNMENT );
			lowerStack.add( components[1] );
			lowerStack.add( components[2] );

			// Add the column of bars in its proper place.
			top.add( components[0] );
			bottom.add( lowerStack );
		}

		// Set the graph's proper size.
		graph.setPreferredSize( new Dimension( graph.getWidth(), 400 ) );

		// Set the sequence on the graph.
		seqLabel.setText( backend.getOligoTargetSequence() );

		// Set the input panel invisible and the results panel visible, then
		// resize the window accordingly.
		panelInput.setVisible( false );
		panelResults.setVisible( true );
		pack();

		// Set the results panel name and show the first oligo.
		String initType = getCustomMenus().findByTreePath("Oligo Graph->[0]").getText();
		String stopString = Integer.toString( stop );
		panelResults.setName( initType + ";0/" + stopString );
		invokeCommand( "Menu " + initType, this);
		setOligoData( backend.getGraphRegionBegin() );

		notifyUpdated();
	}

	/**
	 * Create the panel which holds input controls.
	 */
	private void createInputPanel() {

		// Create the file input panel.
		FileField ct = FileField.createDisabled( "CT File" ).inputFile(FileFilters.CT);
		FileField report = FileField.createEnabled( "Report File" ).outputFile(FileFilters.Report);
		FilePanel files = new FilePanel( this, ct, report );
		files.setPanelWidth( 600 );
		files.makePanel();

		// Create the mode selection panel.
		String[] modeOptions = {
			"Break Local Structure",
			"Refold Whole RNA for Each Oligomer (Slowest)",
			"Do Not Consider Target Structure (Fastest)"
		};
		RadioButtonPanel modes =
			RadioButtonPanel.makeVertical( "Mode", modeOptions );

		// Create the oligo chemistry panel.
		String[] chemOptions = { "DNA", "RNA" };
		RadioButtonPanel chem = RadioButtonPanel.makeHorizontal(
			"Oligomer Chemistry", chemOptions );
		int chemPanelHeight = chem.getPreferredSize().height;
		chem.setPreferredSize( new Dimension( 200, chemPanelHeight ) );

		// Create the check box that gives the suboptimal structures option.
		HTMLCheckBox suboptimal = HTMLCheckBox.createEmptyBox(
			"Include Target Suboptimal Structures in " +
			"Free Energy Calculation" );

		// Create the oligo length panel. Also add a change listener to the
		// spinner so it updates the other spinners when it changes.
		SpinnerPanel oligoLength = new SpinnerPanel(
			"Oligo Length:", 18, 1, Integer.MAX_VALUE );
		oligoLength.getSpinner().addChangeListener( new ChangeListener() {
			public void stateChanged( ChangeEvent e ) {
				invokeCommand( "Update Target Region Range", this);
			}
		});

		// Create the oligo concentration panel.
		IntegerField conc = new IntegerField( "Oligo Concentration", 1, 1 );
		FieldPanel concField = new FieldPanel( conc );
		concField.setPanelWidth( 200 );
		concField.makePanel();
		StringComboBox concBox = new StringComboBox( "M", "mM", "uM", "nM", "pM" );
		concBox.setSelectedItem( "uM" );
		JPanel concPanel = new JPanel();
		concPanel.add( concField );
		concPanel.add( concBox );

		// Create the target walk limits panel.
		BorderBuilder limitBuilder = new BorderBuilder();
		SpinnerPanel start = new SpinnerPanel( "Start:", 0, 0, 0 );
		limitBuilder.makeRightBorder( 5, start );
		SpinnerPanel stop = new SpinnerPanel( "Stop:", 0, 0, 0 );
		JPanel limitsPanel = new JPanel();
		limitBuilder.makeTitledBorder(
			"Target Structure Limits for Walk", limitsPanel );
		limitsPanel.add( start );
		limitsPanel.add( stop );

		// Create the start button.
		JButton startButton = new JButton( "START" );
		startButton.addActionListener( this );

		// Add the input controls to the window and make sure the input panel
		// is visible.
		panelInput = new JPanel( new GridBagLayout() );
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.fill = GridBagConstraints.CENTER;
		constraints.gridwidth = 2;
		constraints.gridheight = 1;
		constraints.gridx = 0;
		constraints.gridy = 0;
		panelInput.add( files, constraints );
		constraints.gridy = 1;
		panelInput.add( modes, constraints );
		constraints.insets = new Insets( 5, 0, 0, 0 );
		constraints.gridy = 2;
		panelInput.add( chem, constraints );
		constraints.insets = new Insets( 0, 0, 0, 0 );
		constraints.gridy = 3;
		panelInput.add( suboptimal, constraints );
		constraints.gridy = 4;
		panelInput.add( oligoLength, constraints );
		constraints.gridy = 5;
		panelInput.add( concPanel, constraints );
		constraints.gridwidth = 1;
		constraints.insets = new Insets( 0, 10, 10, 0 );
		constraints.gridy = 6;
		panelInput.add( limitsPanel, constraints );
		constraints.insets = new Insets( 0, 0, 10, 0 );
		constraints.gridx = 1;
		constraints.ipadx = 15;
		constraints.ipady = 15;
		panelInput.add( startButton, constraints );
		panelInput.setVisible( true );
	}

	/**
	 * Create the panel which holds results data.
	 * <br><br>
	 * The panel created by this method only builds the general layout of the
	 * results window; it doesn't include real results.
	 */
	private void createResultsPanel() {

		// Create the labels panel with empty labels to begin with.
		JPanel allLabels = new JPanel( new GridLayout( 3, 3 ) );
		for( int i = 1; i <= 9; i++ ) {
			JLabel label = new JLabel();
			label.setFont( new Font( Font.MONOSPACED, 1, 12 ) );
			label.setHorizontalAlignment( JLabel.CENTER );
			allLabels.add( label );
		}

		// Create the panel that moves the oligo bar graph by looping over
		// the possible movement increments.
		JPanel buttonPanel = new JPanel( new GridLayout( 1, 5 ) );
		String[] buttonTexts = { "<<", "<", "Go...", ">", ">>" };
		Integer[] increments = { -10, -1, 0, 1, 10 };
		for( int i = 1; i <= buttonTexts.length; i++ ) {

			// Get the text and oligo increment for the next button.
			String text = buttonTexts[i-1];
			String increment = Integer.toString( increments[i-1] );

			// Create the next button.
			JButton button = new JButton( text );
			button.setActionCommand( "Move Oligo " + increment );
			button.addActionListener( this );

			// Create a padding panel and add the button to it to create a
			// single complete button panel, then add that panel to the main
			// button panel.
			JPanel mini = new JPanel();
			new BorderBuilder().makeEqualBorder( 5, mini );
			mini.add( button );
			buttonPanel.add( mini );
		}
		buttonPanel.setPreferredSize(
			new Dimension( buttonPanel.getPreferredSize().width, 100 ) );

		// Create the graph panel.
		Box graph = Box.createVerticalBox();
		Box top = Box.createHorizontalBox();
		top.setAlignmentX( Box.LEFT_ALIGNMENT );
		Box bottom = Box.createHorizontalBox();
		bottom.setAlignmentX( Box.LEFT_ALIGNMENT );
		graph.add( top );
		graph.add( bottom );

		// Create the panel that holds the graph legend.
		JPanel legend = new JPanel( new BorderLayout() );
		legend.add( new JLabel( "MAX" ), BorderLayout.NORTH );
		legend.add( new JLabel( "MIN" ), BorderLayout.SOUTH );
		Component[] legendLabels = legend.getComponents();
		for( Component label: legendLabels ) {
			((JLabel)label).setFont( new Font( Font.MONOSPACED, 1, 12 ) );
			((JLabel)label).setHorizontalAlignment( JLabel.RIGHT );
			new BorderBuilder().makeRightBorder( 5, (JLabel)label );
		}

		// Create the sequence label.
		JLabel sequence = new JLabel( "SEQUENCE" );
		sequence.setFont( new Font( Font.MONOSPACED, 1, 16 ) );
		sequence.setHorizontalAlignment( JLabel.LEFT );

		// Create the oligo label.
		final JLabel oligo = new JLabel( "3OLIGO5" );
		oligo.setFont( new Font( Font.MONOSPACED, 1, 16 ) );
		oligo.setHorizontalAlignment( JLabel.LEFT );
		oligo.addMouseListener( new MouseAdapter() {

			@Override
			public void mousePressed( MouseEvent e ) {

				// Get the current oligo data.
				String oligoSequence =
					oligo.getText().trim()
					.replace( "3", "" ).replace( "5", "" );
				Integer oligoIndex =
					Integer.parseInt( oligo.getName().split( " " )[1] );
				boolean oligoIsRNA =
					backend.getOligoLabelData( oligoIndex )
					.contains( "RNA" );

				// Determine if the oligo can be folded.
				// If the oligo cannot be folded, return.
				boolean oligo = backend.canFoldOligoOligo( oligoIndex );
				boolean self = backend.canFoldOligoSelf( oligoIndex );
				if( !oligo && !self ) { return; }

				// Create a temporary structure file name to use later.
				String tempFile =
					System.getProperty( "user.home" ) +
					System.getProperty( "file.separator" ) + "Oligo" +
					Integer.toString( oligoIndex ) + ".ct";

				// If only oligo-self data is present, prepare a drawing of a
				// unimolecular oligo.
				if( !oligo ) { //self is true
					backend.foldOligo(
						oligoSequence, oligoIndex, false,
						oligoIsRNA, tempFile );
				}

				// If only oligo-oligo data is present, prepare a drawing of a
				// bimolecular oligo.
				else if( !self ) { // oligo is true (due to 'if' condition)
					backend.foldOligo(
						oligoSequence, oligoIndex, true,
						oligoIsRNA, tempFile );
				}

				// If both oligo-self and oligo-oligo data are present, show a
				// dialog that allows the user to choose which one to draw.
				else  {  //oligo and self are both true.
					// Create the dialog.
					// If the user decided not to draw structures, return.
					int choice = JOptionPane.showOptionDialog(
						null, "Would you like to draw structures?",
						"RNAstructure", JOptionPane.YES_NO_CANCEL_OPTION,
						JOptionPane.QUESTION_MESSAGE, null,
						new String[]{ "Unimolecular", "Bimolecular", "Cancel" },
						"Cancel" );
					if( choice == JOptionPane.CANCEL_OPTION ) { return; }

					// Determine which type of folding the user wants.
					boolean bimolecular = ( choice == JOptionPane.NO_OPTION );

					// Prepare a structure.
					backend.foldOligo(
						oligoSequence, oligoIndex, bimolecular, oligoIsRNA,
						tempFile );
				}

				showDrawProgress();
				// Attempt to draw a structure.
				// If that isn't successful, show an error.
				DrawingWindow drawing = new DrawingWindow( tempFile );
				if(!drawing.isError()) { drawing.showWindow(); }

				// Remove the temporary structure file if it exists.
				new File( tempFile ).delete(); //delete does not raise an error if the file does not exist.
			}
		});

		// Put all the scrollable graph components into a container.
		JPanel scrollableGraph = new JPanel( new BorderLayout() );
		scrollableGraph.add( sequence, BorderLayout.NORTH );
		scrollableGraph.add( oligo, BorderLayout.CENTER );
		scrollableGraph.add( graph, BorderLayout.SOUTH );

		// Add the scrollable graph panel to a scroll pane.
		ScrollerPane scroll = new ScrollerPane( scrollableGraph, 875, 460 );
		scroll.setBarPolicies( ScrollerPane.ALWAYS, ScrollerPane.NEVER );
		scroll.setAutoscrolls( true );

		// Put the graph and its legend together into a panel.
		JPanel completeGraph = new JPanel( new BorderLayout() );
		completeGraph.add( legend, BorderLayout.WEST );
		completeGraph.add( scroll );
		new BorderBuilder().makeEqualBorder( 10, completeGraph );

		// Place all the major components where they should be, and make sure
		// the results panel is invisible.
		panelResults = new JPanel( new BorderLayout() );
		panelResults.add( allLabels, BorderLayout.NORTH );
		panelResults.add( buttonPanel, BorderLayout.CENTER );
		panelResults.add( completeGraph, BorderLayout.SOUTH );
		panelResults.setVisible( false );
	}

	@Override
	protected void processCommand(final CommandInfo ci) {

		// If the window is an input window, do any actions that handle input.
		if( isInputWindow ) {

			// Get any input controls from the window that handle actions.
			FilePanel files = (FilePanel)getInputControl( panelInput, 1 );
			Integer length =
				((SpinnerPanel)getInputControl( panelInput, 5 )).getValue();
			JPanel hybrid = (JPanel)getInputControl( panelInput, 7 );
			SpinnerPanel start = (SpinnerPanel)getInputControl( hybrid, 1 );
			SpinnerPanel stop = (SpinnerPanel)getInputControl( hybrid, 2 );

			// If the action comes from the "CT File" button, get a CT file,
			// initialize a data structure, and create an output file name.
			if( ci.command.equals( "CT File" ) ) {

				// Attempt to select the file.
				// If no file was selected, return.
				String file = StandardFileChooser.getOpenName(FileFilters.CT);
				if( file == null ) { return; }

				// Create a data structure.
				// If an error occurred creating the data structure, show an
				// error and return.
				String result =
					backend.buildOligoWalkDataStructure( file, isRNA );
				if( !verifyBackendResult(result, "File: %s\nRNA: %s", file, isRNA) ) { return; }

				// Set the target region range.
				invokeCommand( "Update Target Region Range", ci.event);

				// Set the CT file name and the default output file name in
				// the input panel.
				// Then, enable the menus.
				String defaultOut = getOutputFile( file, "rep" );
				files.setFile( 1, file );
				files.setFile( 2, defaultOut );
				getCustomMenus().enableMenus();
			}

			// If the action comes from the "Report File" button, try to
			// select a report file, and if one was selected set its name.
			else if( ci.command.equals( "Report File" ) ) {
				int index = 2;
				String file = StandardFileChooser.getSaveName(FileFilters.Report, files.getFile(index));
				if( file != null ) { files.setFile( index, file ); }
			}

			// If the command is to update the target region range, do it
			// based on the sequence length.
			else if( ci.command.equals( "Update Target Region Range" ) ) {
				int max = backend.determineOligoMaximum( length );
				start.resetModel( 1, 1, max );
				stop.resetModel( max, 1, max );
			} else
				super.processCommand(ci);
		}

		// Otherwise, do any actions for the results window.
		else {

			// Get the pieces of the graph.
			JPanel graphPanel = (JPanel)panelResults.getComponent( 2 );
			JPanel legend = (JPanel)graphPanel.getComponent( 0 );
			ScrollerPane scroller = (ScrollerPane)graphPanel.getComponent( 1 );
			JPanel graphWrapper = (JPanel)scroller.getViewport().getView();
			Box graph = (Box)graphWrapper.getComponent( 2 );
			Box top = (Box)graph.getComponent( 0 );
			Box bottom = (Box)graph.getComponent( 1 );

			// If the command comes from the results menu, set the name of the
			// results panel to the name of that menu item, then view the
			// proper bars.
			if( ci.command.startsWith( "Menu" ) ) {

				// Create variables to hold the bar colors and the string that
				// holds a particular type of oligo data.
				Color color1 = null;
				Color color2 = null;

				// Remove any colors from the label panel.
				JPanel labelPanel = (JPanel)panelResults.getComponent( 0 );
				int numLabels = labelPanel.getComponents().length;
				for( int i = 1; i <= numLabels; i++ ) {
					labelPanel.getComponent( i - 1 )
						.setForeground( Color.BLACK );
				}

				// Set the panel name to reflect the new panel being shown.
				String name = panelResults.getName();
				String oligoData = name.substring( name.indexOf( ";" ) + 1 );
				String panelType = ci.command.substring( 5 );
				panelResults.setName( panelType + ";" + oligoData );

				// Get all possible oligo data, and declare a string that will
				// be used for only the selected type of data.
				int graphHeight =
					((JPanel)scroller.getViewport().getView())
					.getComponent( 2 ).getPreferredSize().height;
				String[] allOligoData =
					backend.getAllOligoData( graphHeight ).split( ";" );
				String oligoTypeData = null;

				// Set colors and bounds depending on the results panel type.
				// Also set the appropriate label colors after they are
				// properly identified.
				// And identify the proper data string for each type.
				if( ci.command.endsWith( "Overall and Duplex" ) ) {
					color1 = Color.BLUE;
					color2 = new Color( 17, 167, 77 );
					oligoTypeData = allOligoData[0];
					labelPanel.getComponent( 1 ).setForeground( color1 );
					labelPanel.getComponent( 4 ).setForeground( color2 );
				} else if( ci.command.endsWith( "Overall" ) ) {
					color1 = Color.BLUE;
					oligoTypeData = allOligoData[1];
					labelPanel.getComponent( 1 ).setForeground( color1 );
				} else if( ci.command.endsWith( "Duplex" ) ) {
					color1 = new Color( 17, 167, 77 );
					oligoTypeData = allOligoData[2];
					labelPanel.getComponent( 4 ).setForeground( color1 );
				} else if( ci.command.contains( "Target" ) ) {
					color1 = new Color( 114, 67, 212 );
					oligoTypeData = allOligoData[3];
					labelPanel.getComponent( 2 ).setForeground( color1 );
				} else if( ci.command.endsWith( "Unimolecular" ) ) {
					color1 = Color.MAGENTA;
					oligoTypeData = allOligoData[4];
					labelPanel.getComponent( 5 ).setForeground( color1 );
				} else if( ci.command.endsWith( "Bimolecular" ) ) {
					color1 = new Color( 139, 25, 14 );
					oligoTypeData = allOligoData[5];
					labelPanel.getComponent( 8 ).setForeground( color1 );
				} else
					super.processCommand(ci);

				// Get the bar size data for the appropriate panel.
				// The first two indices in this array are the minimum and
				// maximum bounds of the panel, and the others are groups of
				// lengths for each column of bars in the graph.
				String[] sizeData = oligoTypeData.trim().split( " " );

				// Set the maximum and minimum bounds on the legend.
				((JLabel)legend.getComponent( 0 )).setText( sizeData[1] );
				((JLabel)legend.getComponent( 1 )).setText( sizeData[0] );

				// Set the sizes and colors of the bars.
				// At the same time, save the maximum height of any bars > 0.
				Double maxBarSize = 0.0;
				int start = backend.getGraphRegionBegin();
				int end = backend.getGraphRegionEnd();
				for( int i = start; i <= end; i++ ) {

					// Get the lengths of the upper, middle, and lower bars.
					int lengthIndex = ( i - start ) + 2;
					String[] lengths = sizeData[lengthIndex].split( "," );
					Double upper = Double.parseDouble( lengths[0] );
					Double lower1 = Double.parseDouble( lengths[1] );
					Double lower2 = Double.parseDouble( lengths[2] );

					// Save the maximum bar.
					maxBarSize = Math.max( maxBarSize, upper );

					// Get the bars for this index.
					BarButton upperBar = (BarButton)top.getComponent( i );
					Box lowerBox = (Box)bottom.getComponent( i );
					BarButton lower1Bar = (BarButton)lowerBox.getComponent( 0 );
					BarButton lower2Bar = (BarButton)lowerBox.getComponent( 1 );

					// Set the size of the upper bar.
					Dimension upperSize =
						new Dimension( 10, upper.intValue() );
					upperBar.setPreferredSize( upperSize );
					upperBar.setMinimumSize( upperBar.getPreferredSize() );
					upperBar.setMaximumSize( upperBar.getPreferredSize() );
					upperBar.setSize( upperBar.getPreferredSize() );

					// Set the size of the middle bar.
					Dimension lower1Size =
						new Dimension( 10, lower1.intValue() );
					lower1Bar.setPreferredSize( lower1Size );
					lower1Bar.setMinimumSize( lower1Bar.getPreferredSize() );
					lower1Bar.setMaximumSize( lower1Bar.getPreferredSize() );
					lower1Bar.setSize( lower1Bar.getPreferredSize() );

					// Set the size of the lower bar.
					Dimension lower2Size =
						new Dimension( 10, lower2.intValue() );
					lower2Bar.setPreferredSize( lower2Size );
					lower2Bar.setMinimumSize( lower2Bar.getPreferredSize() );
					lower2Bar.setMaximumSize( lower2Bar.getPreferredSize() );
					lower2Bar.setSize( lower2Bar.getPreferredSize() );

					// Set the color of the bars, assuming no stacked bars.
					upperBar.setBarColor( color1 );
					lower1Bar.setBarColor( color1 );
					lower2Bar.setBarColor( color1 );

					// Color any bars that are stacked, if applicable.
					if( color2 != null ) {
						if( lower2 != 0.0 ) { lower2Bar.setBarColor( color2 ); }
						else { lower1Bar.setBarColor( color2 ); }
					}
				}

				// Set the height of the upper graph panel.
				Dimension topSize = new Dimension(
					top.getPreferredSize().width, maxBarSize.intValue() );
				top.setPreferredSize( topSize );
				top.setMinimumSize( topSize );
				top.setMaximumSize( topSize );
				top.setSize( topSize );

				// Set the current oligo selected.
				String currentString = name.substring(
					name.indexOf( ";" ) + 1, name.indexOf( "/" ) );
				int currentIndex = Integer.parseInt( currentString );
				setOligoData( currentIndex );

				// Refresh the panel.
				panelResults.revalidate();
				panelResults.repaint();
			}

			// If the command comes from an oligo mover button, move to the
			// next proper oligo.
			else if( ci.command.startsWith( "Move Oligo" ) ) {
				String name = panelResults.getName();
				String currentString = name.substring(
					name.indexOf( ";" ) + 1, name.indexOf( "/" ) );
				int currentIndex = Integer.parseInt( currentString );
				int moveAmount = Integer.parseInt( ci.command.split( " " )[2] );
				if( moveAmount == 0 ) { new OligoSelectionDialog(); }
				else { setOligoData( currentIndex + moveAmount ); }
			}

			// If the command comes from a bar button, select it.
			else if( ci.command.startsWith( "Set Oligo" ) ) {
				int index = Integer.parseInt( ci.command.split( " " )[2] );
				setOligoData( index );
			} else
				super.processCommand(ci);
		}
	}

	/**
	 * {@inheritDoc}
	 * <br><br>
	 * Note that this version of the function does NOT use the placeComponent
	 * method, because it builds on a panel rather than the content pane.
	 * <br><br>
	 * Unlike other derived classes of ModuleWindow, it holds two panels, an
	 * input panel and a results panel, each of which are visible when needed.
	 */
	@Override
	protected void makeInputControls() {

		// Create the input panel and the results panel.
		// Note that the results panel will not be fully filled in with data
		// until a calculation is done.
		createInputPanel();
		createResultsPanel();

		// Add both the input panel and the results panel to the window.
		add( panelInput );
		add( panelResults );
	}

	@Override
	protected boolean runMainCalculation() {

		// Get data from the file input panel.
		// If an error occurred while retrieving data, return.
		FilePanel files = (FilePanel)getInputControl( panelInput, 1 );
		files.saveRecent();

		files.getFile( 1 );
		String reportFile = files.getFile( 2 );
		if( files.isError() ) { return false; }

		// Get the calculation mode.
		RadioButtonPanel modePanel =
			(RadioButtonPanel)getInputControl( panelInput, 2 );
		int modeNum = modePanel.getSelectedIndex();

		// Get the chemistry type, the suboptimal structures mode, and the
		// oligo length.
		String chem =
			((RadioButtonPanel)getInputControl( panelInput, 3 ))
			.getSelectedName();
		boolean suboptimal =
			((HTMLCheckBox)getInputControl( panelInput, 4 )).isSelected();
		Integer length =
			((SpinnerPanel)getInputControl( panelInput, 5 )).getValue();

		// Determine the oligo concentration.
		JPanel concPanel = (JPanel)getInputControl( panelInput, 6 );
		Integer concAmount =
			((IntegerField)((FieldPanel)getInputControl( concPanel, 1 ))
			.getField( 1 )).getValue();
		String concUnit =
			((StringComboBox)getInputControl( concPanel, 2 ))
				.getSelectedItem();

		// Get the start and stop indices of hybridization.
		JPanel hybrid = (JPanel)getInputControl( panelInput, 7 );
		Integer start =
			((SpinnerPanel)getInputControl( hybrid, 1 )).getValue();
		Integer stop =
			((SpinnerPanel)getInputControl( hybrid, 2 )).getValue();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgress();
		String result =
			backend.runOligoWalk(
				reportFile, modeNum, chem, suboptimal, length, concAmount,
				concUnit, start, stop );
		if (!displayCalcError(result)) return false;

		// Rebuild the current window as a results window.
		buildResultsWindow(
			reportFile, backend.getOligoTargetLength(), start, stop );
		return false; // always return false because this window should not be closed.
	}

	@Override
	protected MergeMenu[] createCustomMenus() {

		// If the window is an input window, create a temperature menu.
		// Otherwise, set a specialized graph switching menu.
		if( isInputWindow ) {
			ConstraintsMenu temperature = new ConstraintsMenu( backend );
			temperature.buildTemperatureMenu();
			return new ConstraintsMenu[]{ temperature };
		} else {
			OligoWalkResultsMenu resultsMenu = new OligoWalkResultsMenu();
			return new MergeMenu[]{ resultsMenu };
		}
	}

	/**
	 * Set the data for a particular oligo.
	 *
	 * @param index   The index of the oligo to set.
	 */
	private void setOligoData( int index ) {

		// Get the current name of the results panel, which holds useful info
		// about its current state.
		String name = panelResults.getName();

		// Make sure the index isn't out of bounds, and if it is, reset it
		// to be within bounds.
		// Then, convert the index into a string.
		int numOligos =
			Integer.parseInt( name.substring( name.indexOf( "/" ) + 1 ) );
		if( index < 1 ) { index = 1; }
		else if( index > numOligos ) { index = numOligos; }
		String indexString = Integer.toString( index );

		// Get data for the requested oligo.
		// Note that this data contains all text and all values necessary for
		// the labels, formatted as HTML.
		String[] dataArray = backend.getOligoLabelData( index ).split( ";" );
		int numLabels = dataArray.length;
		JPanel labelPanel = (JPanel)panelResults.getComponent( 0 );
		for( int i = 1; i <= numLabels; i++ ) {
			((JLabel)labelPanel.getComponent( i-1 )).setText( dataArray[i-1] );
		}

		// Set the oligo text and name.
		JPanel graphPanel = (JPanel)panelResults.getComponent( 2 );
		ScrollerPane scroller = (ScrollerPane)graphPanel.getComponent( 1 );
		JPanel graphWrapper = (JPanel)scroller.getViewport().getView();
		JLabel oligoLabel = (JLabel)graphWrapper.getComponent( 1 );
		oligoLabel.setText( backend.getDisplayedOligo( index ) );
		oligoLabel.setName( "Oligo " + Integer.toString( index ) );

		// Go through all the oligos.
		// If the oligo is the specified index, select it; if not, deselect it.
		Box graph = (Box)graphWrapper.getComponent( 2 );
		Box top = (Box)graph.getComponent( 0 );
		Box bottom = (Box)graph.getComponent( 1 );
		int start = backend.getGraphRegionBegin();
		int end = backend.getGraphRegionEnd();
		for( int i = start; i <= end; i++ ) {

			// Get the upper two possible bars for this oligo.
			// The lowest possible bar is never highlighted for selection.
			BarButton upperBar = (BarButton)top.getComponent( i );
			Box lowerBox = (Box)bottom.getComponent( i );
			BarButton lower1Bar = (BarButton)lowerBox.getComponent( 0 );

			// If the oligo index is not the specified one, deselect the oligo.
			// Otherwise, select it.
			if( i != index ) {
				if( upperBar.getHeight() > 0.0 ) { upperBar.deselect(); }
				else { lower1Bar.deselect(); }
			} else {
				if( upperBar.getHeight() > 0.0 ) { upperBar.select(); }
				else { lower1Bar.select(); }
			}
		}

		// Reset the name of the results panel.
		String oldIndex =
			name.substring( name.indexOf( ";" ) + 1, name.indexOf( "/" ) ); 
		panelResults.setName( name.replaceFirst( oldIndex, indexString ) );
	}

	/**
	 * An inner class that encapsulates a bar button.
	 *
	 * @author Jessica Reuter
	 */
	class BarButton
		extends JButton {
		private static final long serialVersionUID = 20120802;

		/**
		 * The current color painted on this button.
		 */
		private Color currentColor;

		/**
		 * The default color forthis button.
		 */
		private Color defaultColor;

		/**
		 * The selection color for this button.
		 */
		private Color selectionColor;

		/**
		 * Constructor.
		 */
		public BarButton() {
			currentColor = Color.white;
			defaultColor = Color.white;
			selectionColor = Color.red;
		}

		@Override
		public void paintComponent( Graphics g ) {
			super.paintComponent( g );
			g.setColor( currentColor );
			((Graphics2D)g).fill( getVisibleRect() );
		}

		/**
		 * Deselect this button.
		 */
		public void deselect() {
			currentColor = defaultColor;
			repaint();
		}

		/**
		 * Select this button.
		 */
		public void select() {
			currentColor = selectionColor;
			repaint();
		}

		/**
		 * Set the default color that will be painted on this button.
		 *
		 * @param color   The color for this button.
		 */
		public void setBarColor( Color color ) {
			currentColor = color;
			defaultColor = color;
		}
	}

	/**
	 * An inner class that handles creation of the oligo graph menu for the
	 * OligoWalk results window.
	 *
	 * @author Jessica S. Reuter
	 */
	class OligoWalkResultsMenu
		extends MainMenu {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 */
		public OligoWalkResultsMenu() {
			super( "Oligo Graph" );

			// Build the menu.
			addItem(
				"Overall and Duplex",
				"Show overall and duplex oligo results." );
			addSeparator();
			addItem( "Overall", "Show overall oligo results." );
			addItem( "Duplex", "Show duplex oligo results." );
			addItem(
				"Broken Target Structure",
				"Show broken target structure oligo results." );
			addItem(
				"Oligomer Unimolecular",
				"Show unimolecular oligo results." );
			addItem(
				"Oligomer Bimolecular",
				"Show bimolecular oligo results." );

			// Set the action commands.
			for( JMenuItem mi: getMenuItems() )
				mi.setActionCommand("Menu " + mi.getName());
		}

		@Override
		protected void onMenuAction(String command, final ActionEvent ev) {
			invokeCommand(command, this);
		}
	}

	/**
	 * An inner class which creates a dialog that selects a particular oligo
	 * to move to in the results window.
	 *
	 * @author Jessica S. Reuter
	 */
	public class OligoSelectionDialog
		extends ValueSelectionDialog {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 */
		public OligoSelectionDialog() {
			super( "Select Oligo" );

			// Get the selected oligo and the total number of oligos, and use
			// those values to build the number field.
			String name = panelResults.getName();
			int selectedOligo = Integer.parseInt(
					name.substring( name.indexOf( ";" ) + 1,
					name.indexOf( "/" ) ) );
			int numOligos =
				Integer.parseInt( name.substring( name.indexOf( "/" ) + 1 ) );
			IntegerField field = new IntegerField(
				"Oligo Number", selectedOligo, 1, numOligos );

			// Build the dialog.
			buildDialog( "OK", field );
		}

		@Override
		public JPanel createAdditionalControls() {

			// Create the button that selects the most stable oligo.
			JButton stable = new JButton( "Most Stable" );
			stable.setPreferredSize( new Dimension( 125, 30 ) );
			stable.addActionListener( new ActionListener() {
				public void actionPerformed( ActionEvent e ) {
					int mostStable = backend.getMostStableOligo();
					fields[0].setText( Integer.toString( mostStable ) );
				}
			});

			// Add the button to a panel and return it.
			JPanel panel = new JPanel();
			new BorderBuilder().makeBottomBorder( 5, panel );
			panel.add( stable );
			return panel;
		}

		@Override
		public ActionListener createSelectionAction() {
			return new ActionListener() {
				public void actionPerformed( ActionEvent e ) {
					setOligoData( Integer.parseInt( fields[0].getText() ) );
					dispose();
				}
			};
		}
	}

	/**
	 * An class that creates a panel containing a single spinner with a label.
	 *
	 * @author Jessica S. Reuter
	 */
	public class SpinnerPanel
		extends JPanel {
		private static final long serialVersionUID = 20120802;

		/**
		 * The spinner in this panel.
		 */
		private JSpinner spinner;

		/**
		 * Constructor.
		 *
		 * @param text       The label text.
		 * @param initial    The initial value of the spinner.
		 * @param min        The minimum value of the spinner.
		 * @param max        The maximum value of the spinner.
		 */
		public SpinnerPanel( String text, int initial, int min, int max ) {

			// Set the panel layout and create a border builder for use later.
			setLayout( new BorderLayout() );
			BorderBuilder borderBuilder = new BorderBuilder();

			// Create the label and put it in place.
			JLabel label = new JLabel( text );
			borderBuilder.makeTopBorder( 5, label );
			label.setHorizontalAlignment( JLabel.RIGHT );
			label.setVerticalAlignment( JLabel.TOP );
			add( label, BorderLayout.WEST );

			// Create the spinner and put it in place.
			JPanel mini = new JPanel();
			borderBuilder.makeTopBorder( -5, mini );
			spinner = new JSpinner();
			((JSpinner.DefaultEditor)spinner.getEditor())
				.getTextField().setHorizontalAlignment( JTextField.LEFT );
			resetModel( initial, min, max );
			spinner.setName(text);
			mini.add( spinner );
			add( mini );

			// Add a key listener to the spinner text field so it can only
			// accept numbers (1234567890) or backspacing.
			((DefaultEditor)spinner.getEditor()).getTextField()
				.addKeyListener( new KeyAdapter() {

				@Override
				public void keyPressed( KeyEvent e ) {
					if( !Character.isDigit( e.getKeyChar() ) ) { e.consume(); }
				}

				@Override
				public void keyTyped( KeyEvent e ) {
					if( !Character.isDigit( e.getKeyChar() ) ) { e.consume(); }
				}
			});
		}

		/**
		 * Get the spinner.
		 *
		 * @return   The spinner.
		 */
		public JSpinner getSpinner() { return spinner; }

		/**
		 * Get the spinner value.
		 *
		 * @return   The spinner value.
		 */
		public Integer getValue() {
			return (Integer)spinner.getValue();
		}

		/**
		 * Reset the model for the spinner.
		 *
		 * @param initial   The initial value of the spinner.
		 * @param min       The minimum value of the spinner.
		 * @param max       The maximum value of the spinner.
		 */
		public void resetModel( int initial, int min, int max ) {
			spinner.setModel( new SpinnerNumberModel( initial, min, max, 1 ) );
			((JSpinner.DefaultEditor)spinner.getEditor())
				.getTextField().setColumns( 5 );
		}
	}

	/**
	 * An inner class that creates a combo box holding a string list.
	 * <br><br>
	 * Note that this class is compatible with Java 6, but not fully compatible
	 * with Java 7. Java 6 does not require JComboBox to be parameterized, but
	 * Java 7 does. Raw type warnings and unchecked type warnings are
	 * suppressed in this class to ensure backward compatibility.
	 *
	 * @author Jessica S. Reuter
	 */
	@SuppressWarnings("rawtypes")
	public class StringComboBox extends JComboBox {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 *
		 * @param choices   The choices available in this box.
		 */
		@SuppressWarnings("unchecked")
		public StringComboBox( String... choices ) {
			for( String element: choices ) { addItem( element ); }
		}

		@Override
		public String getSelectedItem() {
			return super.getSelectedItem().toString();
		}
	}
}
