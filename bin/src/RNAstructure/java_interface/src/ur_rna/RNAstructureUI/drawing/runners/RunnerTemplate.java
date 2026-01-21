/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.drawing.runners;

import javax.swing.*;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

/**
 * A class that templates the way in which the RNAstructure Java drawing
 * framework is run independently.
 *
 * @author Jessica S. Reuter
 */
public abstract class RunnerTemplate {

	/**
	 * Constructor.
	 *
	 * @param libraryName   The name of the native library connected to this
	 *                      template.
	 */
	public RunnerTemplate( String libraryName ) {
		try {
			UIManager.setLookAndFeel(
				UIManager.getSystemLookAndFeelClassName() );
			System.loadLibrary( libraryName );
		} catch( Exception e ) { e.getStackTrace(); }
	}

	/**
	 * Add a listener to the dialog spawned by this template that exits the
	 * application when the dialog is closed.
	 *
	 * @param dialog   The dialog to close.
	 */
	public void addClosingListener( JDialog dialog ) {
		dialog.addWindowListener( new WindowAdapter() {
			public void windowClosed( WindowEvent e ) { System.exit( 0 ); }
		});
	}

	public void addClosingListener( JInternalFrame dialog ) {
		dialog.addInternalFrameListener( new InternalFrameAdapter() {
			@Override
			public void internalFrameClosed(final InternalFrameEvent e) {
				System.exit( 0 );
			}
		});
	}

	/**
	 * Select an image file.
	 *
	 * @return   The image file.
	 */
	protected abstract String selectFile();
}
