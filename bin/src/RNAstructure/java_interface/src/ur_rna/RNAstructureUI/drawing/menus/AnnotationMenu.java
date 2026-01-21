/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.drawing.menus;

import ur_rna.RNAstructureUI.drawing.dialogs.StructureDialog;
import ur_rna.RNAstructureUI.menus.MainMenu;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.FileFilters;
import ur_rna.RNAstructureUI.utilities.FileType;

import java.awt.event.ActionEvent;

/**
 * A class that handles a menu which can annotate a drawn structure.
 *
 * @author Jessica S. Reuter
 */
public class AnnotationMenu
	extends MainMenu {
	private static final long serialVersionUID = 20120802;

	/**
	 * The structure panel connected to this menu.
	 */
	private StructureDialog dialog;

	/**
	 * Constructor.
	 *
	 * @param dialog   The structure dialog connected to this menu.
	 */
	public AnnotationMenu( StructureDialog dialog ) {

		// Create the menu.
		super( "Annotations" );
		addItem(
			"Add Probability Annotation",
			"Add probability color annotation to the structure." );
		addItem(
			"Add SHAPE Annotation",
			"Add SHAPE color annotation to the structure." );
		addItem(
			"Remove Annotation",
			"Remove color annotation from this structure." );

		// Connect the structure panel.
		this.dialog = dialog;
	}

	@Override
	protected void onMenuAction(String command, final ActionEvent ev) {

		// Select a probability data file for annotation, if necessary.
		if( command.contains( "Probability" ) ) {
			String file = StandardFileChooser.getOpenName(FileFilters.PartitionSav);
			if( file != null ) { dialog.setAnnotation( file, FileType.PFS ); }
		}

		// Remove annotation, if necessary.
		else if( command.contains( "Remove" ) ) { dialog.clearAnnotation(); }

		// Select a SHAPE data file for annotation, if necessary.
		else if( command.contains( "SHAPE" ) ) {
			String file = StandardFileChooser.getOpenName(FileFilters.SHAPE);
			if( file != null ) { dialog.setAnnotation( file, FileType.SHAPE ); }
		}
	}
}
