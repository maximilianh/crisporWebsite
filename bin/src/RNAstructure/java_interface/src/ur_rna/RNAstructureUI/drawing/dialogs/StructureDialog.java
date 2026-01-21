/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.drawing.dialogs;

import ur_rna.RNAstructureUI.drawing.menus.AnnotationMenu;
import ur_rna.RNAstructureUI.drawing.menus.StructureMenu;
import ur_rna.RNAstructureUI.drawing.proxy.StructureBackend;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.utilities.FileType;
import ur_rna.Utilities.Strings;

import javax.swing.*;
import java.awt.*;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A dialog that holds structure images.
 *
 * @author Jessica S. Reuter
 */
public class StructureDialog
	extends ImageDialog {
	private static final long serialVersionUID = 20120802;
	private int drawOptions;

	/**
	 * Enum that provides choices for backbone rendering style: Normal (Radial)
	 * flat (linear) or circular.
	 */
	public enum DrawOptions {
		NORMAL, // Normal (Radial) backbone, circled nucleotides, clockwise.
		CIRCULAR, // draw backbone circular
		FLAT, // draw backbone linear
		UNCIRCLED, // Nucleotides circled/uncircled
		FLIPPED // counterclockwise
	}

	/**
	 * Determine whether the specified option is set
	 */
	public boolean getOption(DrawOptions option) {
		return (drawOptions & (1 << option.ordinal())) != 0;
	}

	public boolean toggleOption(DrawOptions option) {
		boolean newVal = !getOption(option);
		setOption(option, newVal);
		return newVal;
	}
	public void setOption(DrawOptions option) {
		setOption(option, true);
	}
	public void setOption(DrawOptions option, boolean value) {
		if (value == getOption(option)) return;
		if (value)
			drawOptions |= (1 << option.ordinal());
		else
			drawOptions &= ~(1 << option.ordinal());

		switch(option) {
			case CIRCULAR:
				if (value) {
					drawOptions &= ~(1<<DrawOptions.FLAT.ordinal()); // mutually exclusive
					backend.setBackboneStyle(2);
				} else
					backend.setBackboneStyle(0);
				break;
			case FLAT:
				if (value) {
					drawOptions &= ~(1<<DrawOptions.CIRCULAR.ordinal()); // mutually exclusive
					backend.setBackboneStyle(1);
				} else
					backend.setBackboneStyle(0);
				break;
			case FLIPPED:
				backend.flip(); //.setFlipped(value);
				break;
			case UNCIRCLED:
				backend.setNucleotidesCircled(!value);
				break;
		}
		redrawStructure();
	}

	/**
	 * The back end object that reads data for structures.
	 */
	private StructureBackend backend;

	/**
	 * The array of data for a particular structure.
	 */
	private String[] data;

	/**
	 * Whether the current structure has pairs.
	 */
	private boolean hasPairs = false;

	/**
	 * The structure number viewed on this panel.
	 */
	private int number;

	/**
	 * The array of regexes used to pull data about a structure.
	 */
	private String[] regex;

	/**
	 * Constructor.
	 * 
	 * @param file   The structure file to draw.
	 */
	public StructureDialog( String file ) {

		// Call the superclass and set the menu bar.
		super( file );
		if (!Strings.isEmpty(error)) return;
		customMenus.setSubItemMergePos(200);
		customMenus.add(new AnnotationMenu( this ), new StructureMenu( this ) );
		buildMenuBar();

		// Create the pattern strings.
		regex = new String[6];
		regex[0] =
			"Placed at\\: " +
			"\\(([0-9]+\\.*[0-9]*)\\,([0-9]+\\.*[0-9]*)\\)";
		regex[1] =
			"Paired to\\: " +
			"([0-9]+)";
		regex[2] =
			"Control point\\: " +
			"\\(([0-9]+\\.*[0-9]*)\\,([0-9]+\\.*[0-9]*)\\)";
		regex[3] =
			"Backbone stretches to: " +
			"\\(([0-9]+\\.*[0-9]*)\\,([0-9]+\\.*[0-9]*)\\)";
		regex[4] =
			"Label placed at: " +
			"\\(([0-9]+\\.*[0-9]*)\\,([0-9]+\\.*[0-9]*)\\)";
		regex[5] =
			"Character\\: ([A-Za-z]).*" +
			"Placed at\\: " +
			"\\(([0-9]+\\.*[0-9]*)\\,([0-9]+\\.*[0-9]*)\\).*" +
			"Color \\(RGB\\)\\: " +
			"([0-9]+)\\,([0-9]+)\\,([0-9]+)";

		// Add a key listener that allows the user to move between structures.
		addKeyListener( new KeyAdapter() {
			@Override
			public void keyReleased( KeyEvent e ) {
				int keyCode = e.getKeyCode();
				int dir = keyCode == KeyEvent.VK_UP ? 1 : keyCode == KeyEvent.VK_DOWN ? -1 : 0;
				int max = getStructureCount();

				boolean isCmdKey = 0 != (e.getModifiers() & Toolkit.getDefaultToolkit().getMenuShortcutKeyMask());

				if( isCmdKey && dir != 0 ) {
					int moveTo = number + dir;
					if( moveTo < 1 )
						moveTo=1;
					else if( moveTo > max )
						moveTo = max;
					setStructureNumber( moveTo );
				}
			}
		});

		// Set the structure number to the first.
		setStructureNumber( 1 );
	}

	/**
	 * Clear annotation from this dialog.
	 */
	public void clearAnnotation() {
		backend.removeAnnotation();
		buildLegend();
		refresh();
	}

	@Override
	protected void createDrawnImage( Graphics2D g2 ) {

		// Set the drawing color to black.
		g2.setColor( Color.black );

		// If the structure doesn't have pairs, paint a message stating this
		// and return.
		if(!hasPairs) {
			g2.setFont( g2.getFont().deriveFont( Font.BOLD, 20 ) );
			g2.drawString( "This structure contains no pairs.", 30, 60 );
			return;
		}

		// Set thick lines for the pairs.
		g2.setStroke( new BasicStroke( (float)3 ) );

		// Draw the pairs, where necessary.
		for( int i = 1; i < data.length; i++ ) {

			// Determine whether the next nucleotide is paired.
			// If it is, get the pair data.
			Matcher m = Pattern.compile( regex[1] ).matcher( data[i] );
			int pair = ( m.find() ) ? Integer.parseInt( m.group( 1 ) ) : 0;
			if( pair != 0 ) {

				// Use matchers to find the coordinates and control point of
				// the pair.
				Matcher m1 = Pattern.compile( regex[0] ).matcher( data[i] );
				m1.find();
				Matcher m2 = Pattern.compile( regex[2] ).matcher( data[i] );
				m2.find();
				Matcher m3 = Pattern.compile( regex[0] ).matcher( data[pair] );
				m3.find();

				// Convert the coordinate values to doubles.
				Double x1 = Double.parseDouble( m1.group( 1 ) );
				Double y1 = Double.parseDouble( m1.group( 2 ) );
				Double centerX = Double.parseDouble( m2.group( 1 ) );
				Double centerY = Double.parseDouble( m2.group( 2 ) );
				Double x2 = Double.parseDouble( m3.group( 1 ) );
				Double y2 = Double.parseDouble( m3.group( 2 ) );

				// Draw the pair as a path.
				Path2D arc = new Path2D.Double();
				arc.moveTo( x1, y1 );
				arc.curveTo( x1, y1, centerX, centerY, x2, y2 );
				g2.draw( arc );
			}
		}

		// Set lines of normal thickness after pairs are done being drawn.
		g2.setStroke( new BasicStroke( (float)1 ) );

		// Draw backbone and nucleotide label lines, where necessary.
		for( int i = 1; i < data.length; i++ ) {

			// Check to make sure coordinates exist for the next data point.
			// If they do exist, draw lines.
			Matcher m1 = Pattern.compile( regex[0] ).matcher( data[i] );
			if( m1.find() ) {

				// Get the next nucleotide's location coordinates, as doubles,
				// to use as an anchor for potential backbone or label lines.
				Double nucX = Double.parseDouble( m1.group( 1 ) );
				Double nucY = Double.parseDouble( m1.group( 2 ) );

				// Use matchers to search for coordinates for backbone or
				// label lines.
				Matcher m2 = Pattern.compile( regex[3] ).matcher( data[i] );
				Matcher m3 = Pattern.compile( regex[4] ).matcher( data[i] );

				// Draw a backbone line, if one should exist.
				if( m2.find() ) {
					Double backX = Double.parseDouble( m2.group( 1 ) );
					Double backY = Double.parseDouble( m2.group( 2 ) );
					g2.draw( new Line2D.Double( nucX, nucY, backX, backY ) );
				}

				// Draw a label line, if one should exist.
				if( m3.find() ) {
					Double labelX = Double.parseDouble( m3.group( 1 ) );
					Double labelY = Double.parseDouble( m3.group( 2 ) );
					g2.draw( new Line2D.Double( nucX, nucY, labelX, labelY ) );
				}
			}
		}

		// Set the font bold and big to draw labels.
		g2.setFont( g2.getFont().deriveFont( Font.BOLD, 28 ) );

		// Draw the nucleotide number labels, where they exist.
		for( int i = 1; i < data.length; i++ ) {
			Matcher m = Pattern.compile( regex[4] ).matcher( data[i] );
			if( m.find() ) {

				// Get the label data.
				Integer x =
					new Double( Double.parseDouble( m.group( 1 ) ) ).intValue();
				Integer y =
					new Double( Double.parseDouble( m.group( 2 ) ) ).intValue();
				String text = Integer.toString(i);

				// Get the bounds of the string to draw, and the visual bounds
				// of the text itself.
				Rectangle stringBounds =
					g2.getFontMetrics().getStringBounds( text, g2 ).getBounds();
				Rectangle visualBounds =
					g2.getFont()
					.createGlyphVector( g2.getFontRenderContext(), text )
					.getVisualBounds().getBounds();

				// Draw the number label.
				g2.setColor( Color.white );
				g2.fillRect(
					x - ( stringBounds.width / 2 ) - 1,
					y - ( visualBounds.height / 2 ) - 1,
					stringBounds.width + 2,
					visualBounds.height + 2 );
				g2.setColor( Color.black );
				g2.drawString( text,
					x - ( stringBounds.width / 2 ),
					y - ( visualBounds.height / 2 ) - visualBounds.y );
			}
		}

		// Draw the nucleotide letter labels.
		int nucCharRadius = 30;
		boolean circled = !getOption(DrawOptions.UNCIRCLED);
		int halfCharRadius = nucCharRadius / 2;
		for( int i = 1; i < data.length; i++ ) {

			// Check to make sure label data exists for the next data point.
			// If it does exist, draw the label.
			Matcher m =
				Pattern.compile( regex[5], Pattern.DOTALL ).matcher( data[i] );
			if( m.find() ) {

				// Get the label data.
				Integer x =
					new Double( Double.parseDouble( m.group( 2 ) ) ).intValue();
				Integer y =
					new Double( Double.parseDouble( m.group( 3 ) ) ).intValue();
				String text = m.group( 1 );

				// Get the bounds of the string to draw, and the visual bounds
				// of the text itself.
				Rectangle stringBounds =
					g2.getFontMetrics().getStringBounds( text, g2 ).getBounds();
				Rectangle visualBounds =
					g2.getFont()
					.createGlyphVector( g2.getFontRenderContext(), text )
					.getVisualBounds().getBounds();

				// Determine the color of the nucleotide.
				Integer red = Integer.parseInt( m.group( 4 ) );
				Integer green = Integer.parseInt( m.group( 5 ) );
				Integer blue = Integer.parseInt( m.group( 6 ) );

				// Draw the letter label.
				g2.setColor( Color.white );
				g2.fillOval(
					x - halfCharRadius, y - halfCharRadius,
					nucCharRadius, nucCharRadius );
				if( circled ) {
					g2.setColor( Color.black );
					g2.drawOval(
						x - halfCharRadius, y - halfCharRadius,
						nucCharRadius, nucCharRadius );
				}
				g2.setColor( new Color( red, green, blue ) );
				g2.drawString( text,
					x - ( stringBounds.width / 2 ),
					y - ( visualBounds.height / 2 ) - visualBounds.y );
			}
		}
	}

	@Override
	protected void readImageData() {
		backend = new StructureBackend();
		error = backend.readStructureData(file);
	}

	/**
	 * Set the annotation for this structure drawing.
	 *
	 * @param file   The data file to get structure information from.
	 */
	public void setAnnotation( String file ) { setAnnotation(file, null); }

	public void setAnnotation( String file, FileType type ) {
		if (type==null) type = FileType.guessTypeFromFilePath(file);
		if (type == null)		{
			Dialogs.showWarning("Could not determine the file type from the file name: " + file);
			return;
		}

		// Read the appropriate annotation.
		boolean correctRead;
		switch (type) {
			case PFS:
				correctRead = backend.addAnnotationProbability( file );
				break;
			case SHAPE:
				correctRead = backend.addAnnotationSHAPE(file);
				break;
			default:
				Dialogs.showError("Cannot add annotation from this type of file: " + type.toString());
				return;
		}

		// If annotation was not read correctly, show an error and return.
		if(!correctRead) {
			String error = "error adding annotation.";
			Dialogs.showError( error );
			return;
		}

		// Build the legend.
		String data = backend.getStructureData( 1 );
		if (Strings.isEmpty(data)) data = "ERROR: Unknown error loading data for structure " + 1;
		if (data.startsWith("ERROR")) {
			Dialogs.showError( data );
			return;
		}
		data = data.substring( data.indexOf( "Legend" ) ).trim();
		data = data.substring( data.indexOf( ":" ) + 2 );
		data = data.replaceAll( "rgb\\(", "" )
			.replaceAll( "\\)", "" )
			.replaceAll( ",", " " )
			.replaceAll( "\"", "" )
			.replaceAll( " -- ", " " );
		String[] legend = data.split( "\\n" );

		// Place the legend and repaint the dialog.
		buildLegend( legend );
		
		
		//DHM ADDED to FIX BUG where annotations don't appear.
		//Note that keys do not appear until windows is moved or structure number is changed...
		
		//Note how current is coming from the structure label!
		JLabel label = (JLabel)getContentPane().getComponent( 0 );
		String[] words = label.getText().split( " " );
		Integer current = Integer.parseInt( words[1] );
		setStructureNumber(current);
		
		//End added by DHM
		
		refresh();
	}

	/**
	 * Set the structure number.
	 *
	 * @param number   The structure number.
	 */
	public void setStructureNumber( int number ) {
		// Set the structure number and get its data.
		// If an error occurred, show a general error and return.
		this.number = number;
		String dataString = backend.getStructureData( number );
		if( Strings.isEmpty(dataString)) dataString = "ERROR: failed to retrieve data for structure " + number + ".";
		if( dataString.startsWith( "ERROR" ) ) {
			Dialogs.showError( dataString );
			return;
		}
		data = dataString.split( "Nucleotide" );

		// Set the bounds.
		String bounds =
			dataString.substring( dataString.indexOf( "Max Bounds" ) );
		int xBound1 = bounds.indexOf( "(" ) + 1;
		int xBound2 = bounds.indexOf( "," );
		int yBound1 = xBound2 + 1;
		int yBound2 = bounds.indexOf( ")" );
		imgSizeX = Double.parseDouble( bounds.substring( xBound1, xBound2 ) );
		imgSizeY = Double.parseDouble( bounds.substring( yBound1, yBound2 ) );

		// Set the label text.
		String structureID =
			dataString.substring( 0, dataString.indexOf( "Nuc" ) ).trim();
		String description =
			dataString.substring( dataString.indexOf( "Description" ) );
		description = description.substring( description.indexOf( ":" ) + 2 );
		if( description.contains( "Legend:" ) ) {
			int legendIndex = description.indexOf( "Legend:" );
			description = description.substring( 0, legendIndex );
		}
		setTopCaption( structureID + " -- " + description.trim() );

		// Determine if the structure has pairs.
		// If not, set max bounds equal to the window size and scale to 50%.
		hasPairs = Pattern.compile( regex[1] ).matcher( dataString ).find();
		if(!hasPairs) {
			setScale( 50 );
			imgSizeY = imgSizeX = preferredDialogWidth;
		}

		// Calculate scale so the whole image fits in the window at the
		// first glance, and zoom.
		double viewWidth = (this.isVisible() ? this.getWidth() :  preferredDialogWidth) - 5;
		double viewHeight = (this.isVisible() ? this.getHeight() :  preferredDialogWidth) - 5;
		double xScale = viewWidth / imgSizeX;
		double yScale = viewHeight / imgSizeY;
		scale = Math.min( xScale, yScale );
		zoomImage();
	}


	void redrawStructure() {
		setStructureNumber(getStructureNumber());
	}

	/**
	 * Write a dot bracket file of all structures on this dialog.
	 *
	 * @param outFile   The file to write.
	 */
	public void writeDotBracketFile( String outFile ) {
		String result = backend.writeDotBracketFile( file, outFile );
		if( result.equals( "" ) ) {
			Dialogs.showMessage( "Dot bracket file written." );
		} else {
			Dialogs.showError( "error writing dot bracket file." );
		}
	}

	/**
	 * Write a helix file from the current structure on this dialog.
	 *
	 * @param outFile   The file to write.
	 */
	public void writeHelixFile( String outFile ) {
		String result = backend.writeHelixFile( file, outFile, number );
		if( result.equals( "" ) ) {
			Dialogs.showMessage( "Helix file written." );
		} else {
			Dialogs.showError( "error writing helix file." );
		}
	}

	/**
	 * Write a Postscript file from the current structure on this dialog.
	 *
	 * @param outFile   The file to write.
	 */
	public void writePostscriptFile( String outFile ) {
		backend.writePostscriptFile( outFile, number );
		Dialogs.showMessage(
			"Postscript file written." );
	}

	/**
	 * Write an SVG file from the current structure on this dialog.
	 *
	 * @param outFile   The file to write.
	 */
	public void writeSVGFile( String outFile ) {
		backend.writeSVGFile( outFile, number );
		Dialogs.showMessage( "SVG file written." );
	}
	public int getStructureNumber() {
		return number;
	}
	public int getStructureCount() {
		return backend.getStructureCount();
	}
}
