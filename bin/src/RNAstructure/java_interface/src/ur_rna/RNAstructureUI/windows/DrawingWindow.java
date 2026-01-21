/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.drawing.dialogs.ImageDialog;
import ur_rna.RNAstructureUI.drawing.dialogs.PlotDialog;
import ur_rna.RNAstructureUI.drawing.dialogs.StructureDialog;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.ui.ScrollerPane;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.*;
import ur_rna.Utilities.OSInfo;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;

/**
 * A class responsible for displaying a drawing window.
 * <br><br>
 * This class is a wrapper for the self-contained drawing framework dialogs,
 * which essentially converts them into internal frames but does not take over
 * or duplicate any of their functions.
 * <br><br>
 * Most drawings, menus, or listeners are transplanted from the drawing
 * framework dialogs.
 *
 * @author Jessica S. Reuter
 */
public class DrawingWindow
	extends InternalWindow {
	private static final long serialVersionUID = 20120802;
	/**
	 * The image dialog connected to this window.
	 */
	private ImageDialog dialog;

	/**
	 * True if an error occurred building this window, false if not.
	 */
	private boolean error;

	/**
	 * The strand being drawn.
	 */
	private int strand;

	/**
	 * File Constructor.
	 *
	 * @param file   The file whose data should be drawn on the window.
	 */
	public DrawingWindow( String file ) {
		this( file, 1 );
	}

	/**
	 * Strand Constructor.
	 * 
	 * @param file     The dot plot data file to draw.
	 * @param strand   The strand from the data file to draw.
	 */
	public DrawingWindow( String file, int strand ) { this(file, strand, null); }
	public DrawingWindow( String file, int strand, FileType drawingType ) {
		// Set the title of this window and the main frame to reflect the
		// creation of a new drawing. Also, set the frame to be resizable.

		setCaption( file );
		setResizable( true );


		// Set the strand.
		this.strand = strand;

		if (drawingType == null) drawingType = FileType.guessTypeFromFilePath(file);


		// If the file is a CT file, get information from a structure dialog.
		// Otherwise, get information from a plot dialog.
		switch (drawingType) {
			case CT:
			case DBN:
				dialog = new StructureDialog( file );
				break;
			case PFS:
			case SAV:
				dialog = new PlotDialog( file, strand );
				break;
			default:
				Dialogs.showWarning("Cannot draw file (due to Unknown file type): " + file);
				return;
		}
		if (dialog.isError())
			return;

		getDialogContent();

		// Add a component listener to the legend so every time it changes, the
		// window is revalidated.
		ComponentListener cListener = new ComponentListener() {

			@Override
			public void componentHidden( ComponentEvent e ) { doRevalidate(); }

			@Override
			public void componentMoved( ComponentEvent e ) { doRevalidate(); }

			@Override
			public void componentResized( ComponentEvent e ) { doRevalidate(); }

			@Override
			public void componentShown( ComponentEvent e ) { doRevalidate(); }
		};
		getContentPane().getComponent( 2 ).addComponentListener( cListener );

		// If no error occurred and the dialog is a structure dialog, grab its
		// structure movement key listener.
		if( (!error) && ( drawingType == FileType.CT ) ) {
			int index = dialog.getKeyListeners().length - 1;
			addKeyListener( dialog.getKeyListeners()[index] );
		}

		// If no error occurred and the dialog is a probability plot dialog,
		// replace the listener on the "Write Probable Structures" menu item
		// with this class, so structures can pop up in DrawingWindows instead
		// of StructureDialogs.
		if( (!error) && ( drawingType == FileType.PFS ) ) {
			JMenuItem m = getCustomMenus().findByTreePath("Output Plot->Write Probable Structures File");
			if (m == null)
				Dialogs.showError("Failed to find the menu item \"Output Plot->Write Probable Structures File\".");
			else {
				ActionListener aListener = m.getActionListeners()[0];
				m.removeActionListener(aListener);
				m.addActionListener(this);
			}
		}
	}

	/**
	 * {@inheritDoc}
	 * <br><br>
	 * Note that most of the actions for this window are handled by the menus
	 * in the self-contained drawing framework.
	 */
	@Override
	protected void processCommand(final CommandInfo ci) {
		// If the command is to write a probable structures file, do so.
		if( ci.command.equals( "Write Probable Structures File" ) ) {

			// Choose and write a structure file.
			String file = StandardFileChooser.getSaveName(FileFilters.CT);
			boolean drawOK = false;
			if (file != null) {
				drawOK = ((PlotDialog) dialog).writeStructuresFile(
						file, false);
			}

			// If the user wants to draw structures, draw them.
			if (drawOK && ModuleWindow.promptDraw()) {
				DrawingWindow window = new DrawingWindow(file);
				if (!window.isError()) { window.showWindow(); }
			}
		} else
			super.processCommand(ci);
	}

	/**                                                                                                                                          
	 * Revalidate the drawing window.                                                                                                            
	 */
	private void doRevalidate() {

		// Get the content pane.
		Container pane = getContentPane();

		// Get the main window components.
		Container label = (Container)pane.getComponent( 0 );
		ScrollerPane mainScroll = (ScrollerPane)pane.getComponent( 1 );
		ScrollerPane legendScroll = (ScrollerPane)pane.getComponent( 2 );
		Container legend = (Container)(legendScroll.getViewport().getView());

		// Calculate the container pane size and set it.
		Dimension size = mainScroll.getSize();
		size.height += label.getPreferredSize().height;
		if( legend.getComponents().length != 0 ) {
			size.height += legendScroll.getPreferredSize().height;
		}
		getContentPane().setPreferredSize( size );

		// Repaint the window.                                                                                                                  
		repaint();
	}

	/**
	 * Get the data from a dialog to put in this window.
	 */
	private void getDialogContent() {

		// If the dialog encountered an error while being built,
		// set the error flag and return.
		if( dialog.isError() ) {
			error = true;
			return;
		}

		// Set the dialog's content pane as the frame's content pane.
		setContentPane( dialog.getContentPane() );

		// Convert the dialog's menus into rollover menus and set them in the
		// dynamic menu bar.
		setCustomMenus(dialog.getCustomMenus());
		getCustomMenus().enableMenus();

		// Grab the dialog's zooming key listener.
		addKeyListener( dialog.getKeyListeners()[0] );

		// Add the dialog's scroll bar listener, if the OS is a type of Mac.
		if( OSInfo.isMac() ) {
			addKeyListener( dialog.getKeyListeners()[1] );
		}
	}

	/**
	 * Get whether this window encountered an error.
	 *
	 * @return   True if an error occurred, false if not.
	 */
	public boolean isError() { return error; }

	public PlotDialog getPlot() { return dialog instanceof PlotDialog ? (PlotDialog)dialog : null; }
	public StructureDialog getStructure() { return dialog instanceof StructureDialog ? (StructureDialog)dialog : null; }
}
