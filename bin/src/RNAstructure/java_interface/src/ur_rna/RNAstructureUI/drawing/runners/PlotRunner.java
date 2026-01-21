/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.drawing.runners;

import ur_rna.RNAstructureUI.drawing.dialogs.PlotDialog;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.io.Serializable;

/**
 * A class that runs RNAstructure Java dot plot drawing independently of the
 * main RNAstructure application.
 *
 * @author Jessica S. Reuter
 */
public class PlotRunner
	extends RunnerTemplate implements Serializable {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 */
	public PlotRunner() { super( "DrawingPlots" ); }

	@Override
	protected String selectFile() {

		// Show a message dialog that tells the user a file must be selected.
		JOptionPane.showMessageDialog( null,
			"Please select a dot plot file to draw using the " +
			"following file chooser." );

		// Create the file chooser and set its defaults.
		JFileChooser chooser = new JFileChooser();
		chooser.setFileSelectionMode( JFileChooser.FILES_ONLY );
		chooser.setAcceptAllFileFilterUsed( false );

		// Add the file filters.
		chooser.addChoosableFileFilter( new FileNameExtensionFilter(
			"Partition Function Save Files (*.pfs)", "pfs" ) );
		chooser.addChoosableFileFilter( new FileNameExtensionFilter(
			"Folding Save Files (*.sav)", "sav" ) );
		chooser.addChoosableFileFilter( new FileNameExtensionFilter(
			"Dynalign Save Files (*.dsv)", "dsv" ) );
		chooser.addChoosableFileFilter( new FileNameExtensionFilter(
			"Dot Plot Files (*.dp)", "dp" ) );

		// Show the file dialog, and then if no file was selected, return the
		// empty string.
		int approvalValue = chooser.showOpenDialog( null );
		if( approvalValue != JFileChooser.APPROVE_OPTION ) { return ""; }

		// Get the file, and make sure it ends with the proper extension.
		String file = chooser.getSelectedFile().getAbsolutePath().trim();
		String description = chooser.getFileFilter().getDescription();
		int start = description.lastIndexOf( "." ) + 1;
		int end = description.lastIndexOf( ")" );
		if( ( start != -1 ) && ( end != -1 ) ) {
			String extension = description.substring( start, end ).trim();
			if( !file.endsWith( extension ) ) {
				file = file.concat( "." + extension );
			}
		}

		// Return the file name.
		return file;
	}

	/**
	 * The main method.
	 *
	 * @param args   The command line arguments (ignored).
	 */
	public static void main( String args[] ) {
		PlotRunner template = new PlotRunner();
		String file = template.selectFile();
		if( !file.equals( "" ) ) {
			PlotDialog dialog = new PlotDialog( file );
			template.addClosingListener( dialog );
			dialog.viewDialog();
		}
	}
}
