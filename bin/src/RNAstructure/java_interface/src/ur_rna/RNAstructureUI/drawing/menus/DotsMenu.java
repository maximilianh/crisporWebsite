/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.drawing.menus;

import ur_rna.RNAstructureUI.drawing.dialogs.PlotDialog;
import ur_rna.RNAstructureUI.menus.MainMenu;
import ur_rna.RNAstructureUI.ui.BorderBuilder;
import ur_rna.RNAstructureUI.ui.FieldPanel;
import ur_rna.RNAstructureUI.ui.NumberField.DoubleField;
import ur_rna.RNAstructureUI.ui.NumberField.IntegerField;
import ur_rna.RNAstructureUI.ui.ValueSelectionDialog;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * A class that handles a menu which can manipulate a drawn dot plot.
 *
 * @author Jessica S. Reuter
 */
public class DotsMenu
	extends MainMenu {
	private static final long serialVersionUID = 20120615;

	/**
	 * The structure dialog connected to this menu.
	 */
	private PlotDialog dialog;

	/**
	 * Constructor.
	 *
	 * @param dialog   The dot plot dialog connected to this menu.
	 */
	public DotsMenu( PlotDialog dialog ) {

		// Create the menu.
		super( "Draw" );
		addItem(
			"Zoom",
			"Zoom this structure." );
		addItem(
			"Choose Colors",
			"Select the number of colors in a dot plot." );
		addItem(
			"Plot Range",
			"Select the range of values displayed in the plot." );

		// Connect the plot dialog.
		this.dialog = dialog;
	}

	@Override
	protected void onMenuAction(String command, final ActionEvent ev) {

		// If the command comes from the color item, select the number of
		// colors that should be in the plot.
		if( command.equals( "Choose Colors" ) ) { new ColorDialog(); }

		// If the command is to change the plot range, do so.
		else if( command.equals( "Plot Range" ) ) { new RangeDialog(); }

		// If the command is to zoom the image, show a zooming dialog.
		else if( command.equals( "Zoom" ) ) { dialog.viewZoomDialog(); }
	}

	/**
	 * An inner class which creates a dialog that sets the number of colors in
	 * a dot plot.
	 *
	 * @author Jessica S. Reuter
	 */
	private class ColorDialog
		extends ValueSelectionDialog {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 */
		public ColorDialog() {

			// Call the superclass.
			super( "Choose Colors" );

			// Get the colors from the back end.
			Integer[] colors = dialog.getColors();

			// Create the colors field and build the dialog with it.
			IntegerField colorsField = new IntegerField(
				"Colors", colors[0], colors[1], colors[2] );
			buildDialog( "OK", colorsField );
		}

		@Override
		public ActionListener createSelectionAction() {
			return new ActionListener() {
				public void actionPerformed( ActionEvent e ) {
					dialog.setColors(
						Integer.parseInt( fields[0].getText() ) );
					dialog.zoomImage();
					dialog.repaint();
					dispose();
				}
			};
		}
	}

	/**
	 * An inner class which creates a dialog that sets the range of values in
	 * a dot plot.
	 *
	 * @author Jessica S. Reuter
	 */
	private class RangeDialog
		extends ValueSelectionDialog {
		private static final long serialVersionUID = 20120802;

		/**
		 * The number panel that holds the range data.
		 */
		private FieldPanel range;

		/**
		 * Constructor.
		 */
		public RangeDialog() {

			// Call the superclass and build the dialog.
			// Note that this dialog doesn't have the usual fields; all input
			// controls are handled outside of this constructor.
			super( "Set Plot Range" );
			buildDialog( "OK" );
		}

		@Override
		public JPanel createAdditionalControls() {

			// Create the range panel.
			Double[] currentBounds = dialog.getCurrentBounds();
			Double[] defaultBounds = dialog.getDefaultBounds();
			DoubleField minField = new DoubleField(
				"Minimum", currentBounds[0],
				defaultBounds[0], defaultBounds[1] );
			DoubleField maxField = new DoubleField(
				"Maximum", currentBounds[1],
				defaultBounds[0], defaultBounds[1] );
			range = new FieldPanel( minField, maxField );
			range.setPanelWidth( 300 );
			range.makePanel();

			// Create the button that resets the plot range.
			JButton reset = new JButton( "Reset" );
			reset.setPreferredSize( new Dimension( 125, 30 ) );
			reset.addActionListener( new ActionListener() {
				public void actionPerformed( ActionEvent e ) {
					Double[] bounds = dialog.getDefaultBounds();
					((DoubleField)range.getField( 1 )).resetField( bounds[0] );
					((DoubleField)range.getField( 2 )).resetField( bounds[1] );
				}
			});
			JPanel resetPanel = new JPanel();
			new BorderBuilder().makeBottomBorder( 5, resetPanel );
			resetPanel.add( reset );

			// Add the components to a panel and return it.
			JPanel panel = new JPanel( new BorderLayout() );
			panel.add( range );
			panel.add( resetPanel, BorderLayout.SOUTH );
			return panel;
		}

		@Override
		public ActionListener createSelectionAction() {
			return new ActionListener() {
				public void actionPerformed( ActionEvent e ) {
					dialog.setRange(
						((DoubleField)range.getField( 1 )).getValue(),
						((DoubleField)range.getField( 2 )).getValue() );
					dialog.zoomImage();
					dialog.repaint();
					dispose();
				}
			};
		}
	}
}
