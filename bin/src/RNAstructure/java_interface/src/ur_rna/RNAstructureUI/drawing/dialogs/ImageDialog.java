/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.drawing.dialogs;

import ur_rna.RNAstructureUI.menus.MainMenu;
import ur_rna.RNAstructureUI.menus.MenuList;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.ui.NumberField.IntegerField;
import ur_rna.RNAstructureUI.ui.ScrollerPane;
import ur_rna.RNAstructureUI.ui.ValueSelectionDialog;
import ur_rna.Utilities.OSInfo;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;

/**
 * A printable and zoomable dialog that holds drawn images.
 * 
 * @author Jessica S. Reuter
 */
public abstract class ImageDialog
	extends JInternalFrame implements Printable {
	private static final long serialVersionUID = 20120802;
	protected MenuList customMenus = new MenuList();

	/**
	 * The error message generated while building this panel, if any.
	 */
	protected String error;

	/**
	 * The file that is drawn in some way on this panel.
	 */
	protected String file;

	/**
	 * The main label that describes the panel.
	 */
	private JLabel label;

	/**
	 * The panel that contains the legend grid.
	 */
	private JPanel legendGrid;

	/**
	 * The main panel that the drawing resides on.
	 */
	protected JPanel mainPanel;

	/**
	 * The maximum X bound of the image.
	 */
	protected double imgSizeX;

	/**
	 * The maximum Y bound of the image.
	 */
	protected double imgSizeY;

	/**
	 * The job that prints out this panel
	 */
	private PrinterJob job = PrinterJob.getPrinterJob();

	/**
	 * The preferred dialog width.
	 */
	protected final int preferredDialogWidth = 450;

	/**
	 * The scale at which the image is seen.
	 */
	protected double scale;
	/**
	 * Constructor.
	 *
	 * @param file    The file that is drawn in some way.
	 */
	@SuppressWarnings("serial")
	protected ImageDialog( String file ) {
		// Set the default bounds and error state.
		imgSizeY = imgSizeX = -Double.MAX_VALUE;
		error = "";

		// Set the dialog defaults.
		setBackground( Color.WHITE );
		setFocusable( true );
		setLayout( new BorderLayout() );
		setDefaultCloseOperation( JDialog.DISPOSE_ON_CLOSE );
		setTitle( file );

		// Add a key listener to the dialog.
		// This will enable the image on the dialog to zoom.
		addKeyListener( new KeyAdapter() {
			private boolean isControl = false;

			@Override
			public void keyPressed( KeyEvent e ) {
				isControl = e.isControlDown();
				if( OSInfo.isMac() ) {
					if( isControl == false ) { isControl = e.isMetaDown(); }
				}
			}

			@Override
			public void keyReleased( KeyEvent e ) {
				int keyCode = e.getKeyCode();
				boolean enter =
					( keyCode == KeyEvent.VK_RIGHT ) ||
					( keyCode == KeyEvent.VK_LEFT );

				if( enter && isControl ) {
					if( keyCode == KeyEvent.VK_RIGHT ) { scale += 0.05; }
					else { scale -= 0.05; }
					setScale( new Double( scale * 100 ).intValue() );
					zoomImage();
				}

				isControl = false;
			}
		});

		// Create the information label.
		label = new JLabel( " " );
		label.setHorizontalAlignment( JLabel.CENTER );
		label.setVerticalAlignment( JLabel.CENTER );
		label.setPreferredSize( new Dimension( preferredDialogWidth, 20 ) );
		label.setBackground( Color.LIGHT_GRAY );
		label.setOpaque( true );
		add( label, BorderLayout.NORTH );

		// Create the image panel.
		mainPanel = new JPanel() {
			public void paintComponent( Graphics g ) {
				super.paintComponent( g );
		
				// Get the graphics context.
				Graphics2D g2 = (Graphics2D)g;
				g2.setRenderingHint( RenderingHints.KEY_INTERPOLATION,
					RenderingHints.VALUE_INTERPOLATION_BILINEAR );

				// Prepare the canvas.
				g2.setColor( Color.white );
				g2.fillRect( 0, 0, getWidth(), getHeight() );
				g2.setColor( Color.black );
				g2.scale( scale, scale );

				// Draw the image.
				createDrawnImage( g2 );
			}
		};
		mainPanel.setOpaque( true );
		mainPanel.setDoubleBuffered( true );

		// Add the image panel to a scroll pane.
		// If the OS is a type of mac, add its Mac-specific key listener.
		final ScrollerPane scroller = new ScrollerPane(
			mainPanel, preferredDialogWidth, preferredDialogWidth );
		scroller.getViewport().setViewPosition( new Point( 0, 0 ) );
		add( scroller, BorderLayout.CENTER );
		if( OSInfo.isMac() ) {
			addKeyListener( scroller.getKeyListeners()[0] );
		}

		// Create the legend pane.
		legendGrid = new JPanel( new GridLayout( 0, 1 ) );
		legendGrid.setOpaque( true );
		legendGrid.setBackground( Color.white );
		ScrollerPane legendScroll =
			new ScrollerPane( legendGrid, preferredDialogWidth, 105 );
		legendScroll.setBarPolicies( ScrollerPane.NEVER, ScrollerPane.ALWAYS );
		buildLegend();
		add( legendScroll, BorderLayout.SOUTH );

		// Read the data necessary to back this dialog.
		this.file = file;
		readImageData();

		// Calculate scale so the whole image fits in the window at the
		// first glance, and zoom.
		double scaleWidth = preferredDialogWidth - 5;
		double xScale = scaleWidth / imgSizeX;
		double yScale = scaleWidth / imgSizeY;
		scale = Math.min( xScale, yScale );
		zoomImage();

		customMenus.add(new PrintMenu());
	}

	private Integer SafeParseInt(String sval, Integer defaultValue) {
		//try {
			return Integer.parseInt(sval);
		//} catch (Throwable ex) {
			//System.out.println("error formatting legend value as Integer: '" + sval + "': " + ex.getMessage());
		//}
//		try {
//			return (int)Double.parseDouble(sval);
//		} catch (Throwable ex) {
//			System.out.println("error formatting legend value as double: '" + sval + "': " + ex.getMessage());
//		}
//		return defaultValue;
	}
	/**
	 * Build the legend grid for this panel.
	 *
	 * @param data   The array of data for the legend.
	 */
	public void buildLegend( String... data ) {
		//System.out.println("Building Legend");

		// Remove all legend components and validate the empty layout.
		legendGrid.removeAll();
		legendGrid.revalidate();

		// If the data array is zero size, specifying an empty legend, return
		// without building a new legend.
		if( data.length == 0 ) {
			legendGrid.getParent().getParent().setVisible( false );
			return;
		}

		// Place each entry in the legend.
		for( String entry: data ) {

			// Get the data from the string.
			String[] entryData = entry.trim().split( " " );
			int entryDataLength = entryData.length;
			Integer r, g, b;
			r = SafeParseInt(entryData[entryDataLength - 3], 0);
			g = SafeParseInt(entryData[entryDataLength - 2], 0);
			b = SafeParseInt(entryData[entryDataLength - 1], 0);
			String text = "";
			for( int i = 1; i <= entryDataLength - 3; i++ ) {
				text += ( entryData[i-1] + " " );
			}
			text.trim();

			// Create a label of the text with the appropriate color, then add
			// it to the legend grid.
			JLabel label = new JLabel( text );
			label.setOpaque( true );
			label.setBackground( Color.WHITE );
			label.setForeground( new Color( r, g, b ) );
			label.setFont( label.getFont().deriveFont( Font.BOLD ) );
			label.setHorizontalAlignment( JLabel.CENTER );
			label.setPreferredSize( new Dimension( preferredDialogWidth, 20 ) );
			label.setMaximumSize( new Dimension( preferredDialogWidth, 20 ) );
			legendGrid.add( label );
		}

		// Show the new legend.
		legendGrid.getParent().getParent().setVisible(true);
	}

	/**
	 * Create the image dialog's menu bar.
	 */
	protected void buildMenuBar() {
		// Create and set an empty menu bar.
		//JMenuBar bar = new JMenuBar();
		//setJMenuBar(bar);

		// TODO: see if some type of menu building is necessary here.
		//customMenus.buildMenu(bar);
		// Add the unique menus first, then the print menu.
		//for( MergeItem menu: customMenus.merged() ) { bar.add( menu.toJMenu() ); }
	}

	public MenuList getCustomMenus() {
		return customMenus;
	}

	/**
	 * Create the image drawn on this panel.
	 *
	 * @param g2   The graphics context.
	 */
	protected abstract void createDrawnImage( Graphics2D g2 );

	/**
	 * Get the file drawn on this dialog.
	 *
	 * @return   The file name.
	 */
	public String getFile() { return file; }

	/**
	 * Get the scale at which the image is zoomed.
	 *
	 * @return   The scale.
	 */
	public double getScale() { return scale; }

	/**
	 * Get whether this panel ran into an error while being built.
	 *
	 * @return   True if an error occurred, false if not.
	 */
	public boolean isError() {
		if( !error.equals( "" ) ) {
			Dialogs.showError( error );
			return true;
		}
		return false;
	}

	/**
	 * Render the panel for the print job.
	 *
	 * @param graphics   The graphics object that is printed.
	 * @param format     The page format.
	 * @param index      The page index
	 */
	public int print( Graphics graphics, PageFormat format, int index ) {

		// The panel only takes up one page, so if the index says we're on page
		// 2, something went wrong.
		if ( index > 0 ) { return NO_SUCH_PAGE; }

		// Print the panel.
		Graphics2D g2d = (Graphics2D)graphics;
		g2d.translate( format.getImageableX(), format.getImageableY() );
		mainPanel.printAll( graphics );
		return PAGE_EXISTS;
	}

	/**
	 * Print if the user decides to, then reset the printer job properties.
	 */
	public void printPanel() {
		job.setPrintable( this);
		if( job.printDialog() ) {
			try { 
				job.print();
				job = PrinterJob.getPrinterJob();
			} catch( PrinterException e ) {
				e.printStackTrace();
			}
		}
	}

	/**
	 * Read image data from a particular file.
	 */
	protected abstract void readImageData();

	/**
	 * Refresh this dialog.
	 */
	public void refresh() {
		Component[] components = getComponents();
		for( Component component: components ) { component.validate(); }
		repaint();
	}

	/**
	 * Set the image scale.
	 *
	 * @param scale   The percent by which the image is scaled,
	 *                as an integer (ex. 5 = 5%).
	 */
	public void setScale( int scale ) { this.scale = (double)scale / 100.0; }

	/**
	 * Set the text on the info label.
	 *
	 * @param text   The text to set on the label.
	 */
	public void setTopCaption( String text ) {
		label.setText(text);
		repaint();
	}

	/**
	 * View the dialog.
	 */
	public void viewDialog() {
		pack();
		//setLocationRelativeTo(null);
		setVisible(true);
	}

	/**
	 * View a zooming dialog.
	 */
	public void viewZoomDialog() {
		new ZoomDialog();
	}

	/**
	 * Zoom the image on the panel.
	 */
	public void zoomImage() {

		// Make sure the scale size isn't too small or too big.
		if( scale < 0.05 ) { scale = 0.05; }
		else if( scale > 10 ) { scale = 10; }

		// Figure out the what the scaled size of the image is.
		Dimension scaled = new Dimension(
			(int)(imgSizeX * scale),
			(int)(imgSizeY * scale) );

		// Set the scaled panel size and repaint.
		mainPanel.setPreferredSize( scaled );
		mainPanel.setMinimumSize( scaled );
		mainPanel.setMaximumSize( scaled );
		mainPanel.setSize( scaled );
		repaint();
		refresh();
	}

	/**
	 * An inner class which creates a print menu.
	 *
	 * @author Jessica S. Reuter
	 */
	protected class PrintMenu
		extends MainMenu {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 */
		public PrintMenu() {
			super( "Print" , 500);
			addItem("Print Image", "Print the current image." );
		}

		@Override
		protected void onMenuAction(String command, final ActionEvent ev) {
			printPanel();
		}
	}

	/**
	 * An inner class which creates a dialog that zooms a structure.
	 *
	 * @author Jessica S. Reuter
	 */
	private class ZoomDialog
		extends ValueSelectionDialog {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 */
		public ZoomDialog() {

			// Call the superclass.
			super( "Zoom" );

			// Create the bounded zoom field and build the dialog with it.
			Integer intScale = new Double( scale * 100.0 ).intValue();
			IntegerField zoomField =
				new IntegerField( "Percent Magnification", intScale, 5, 1000 );
			buildDialog( "OK", zoomField );
		}

		@Override
		public ActionListener createSelectionAction() {
			return new ActionListener() {
				public void actionPerformed( ActionEvent e ) {
					setScale( Integer.parseInt( fields[0].getText() ) );
					zoomImage();
					dispose();
				}
			};
		}
	}
}
