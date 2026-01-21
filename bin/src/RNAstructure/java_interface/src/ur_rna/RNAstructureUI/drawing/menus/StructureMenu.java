/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.drawing.menus;

import ur_rna.RNAstructureUI.drawing.dialogs.StructureDialog;
import ur_rna.RNAstructureUI.menus.MainMenu;
import ur_rna.RNAstructureUI.ui.NumberField.IntegerField;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.ui.ValueSelectionDialog;
import ur_rna.RNAstructureUI.utilities.FileFilters;
import ur_rna.Utilities.swing.CheckMergeItem;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * A class that handles a menu which can manipulate a drawn structure.
 *
 * @author Jessica S. Reuter
 */
public class StructureMenu
	extends MainMenu {
	private static final long serialVersionUID = 20120802;

	/**
	 * The structure dialog connected to this menu.
	 */
	private StructureDialog dialog;
	private CheckMergeItem circular, flat, flipped, circled;

	/**
	 * Constructor.
	 *
	 * @param dialog   The structure dialog connected to this menu.
	 */
	public StructureMenu( StructureDialog dialog ) {

		// Create the menu.
		super( "Draw" );
		this.dialog = dialog;

		addItem(
			"Go to Structure...",
			"Switch to a selected structure." );
		addItem(
			"Zoom",
			"Zoom this structure." );
		addSeparator();
		flipped = addCheckItem(
			"Render Counterclockwise/Flipped",
			"Checked if a structure is rendered Counterclockwise (Flipped Horizontally)." );

		circled = addCheckItem(
			"Render Circles around Nucleotides",
			"Checked if nucleotides are surrounded by circles." );

		circular = addCheckItem(
				"Render Circular Structure",
				"Draw the structure with its backbone stretched around a circle." );
		flat = addCheckItem(
				"Render Linear (Flat) Structure",
				"Draw the structure with a linear backbone." );
		addSeparator();
		addItem(
			"Write Dot Bracket File",
			"Write a dot bracket file of all structures in this window." );
		addItem(
			"Write Helix (Text) File",
			"Write a text file of the helices in this structure." );
		addItem(
			"Write Postscript File",
			"Write the current structure to a Postscript image file." );
		addItem(
			"Write SVG File",
			"Write the current structure to an SVG image file." );
		updateMenuUI();
	}

	@Override
	protected void onMenuAction(String command, final ActionEvent ev) {

		// If the command comes from the circled nucleotides item, switch
		// circling of nucleotides on or off.
		if( command.contains( "Circles" ) ) { dialog.toggleOption(StructureDialog.DrawOptions.UNCIRCLED); updateMenuUI(); }
		else if( command.contains( "Circular" ) ) { dialog.toggleOption(StructureDialog.DrawOptions.CIRCULAR); updateMenuUI();}
		else if( command.contains( "Flat" ) ) { dialog.toggleOption(StructureDialog.DrawOptions.FLAT); updateMenuUI();}
		else if( command.contains( "Flipped" ) ) { dialog.toggleOption(StructureDialog.DrawOptions.FLIPPED); updateMenuUI();}

		// If the command is to write a particular type of output file, do so.
		else if( command.contains( "File" ) ) {

			// Write a dot bracket file, if necessary.
			if( command.contains( "Dot Bracket" ) ) {
				String file = StandardFileChooser.getSaveName(FileFilters.Bracket);
				if( file != null ) { dialog.writeDotBracketFile( file ); }
			}

			// Write a helix file, if necessary.
			else if( command.contains( "Helix" ) ) {
				String file = StandardFileChooser.getSaveName(FileFilters.Helix);
				if( file != null ) { dialog.writeHelixFile( file ); }
			}

			// Write a Postscript file, if necessary.
			else if( command.contains( "Postscript" ) ) {
				String file = StandardFileChooser.getSaveName(FileFilters.Postscript);
				if( file != null ) { dialog.writePostscriptFile( file ); }
			}

			// Write an SVG file, if necessary.
			else if( command.contains( "SVG" ) ) {
				String file = StandardFileChooser.getSaveName(FileFilters.SVG);
				if( file != null ) { dialog.writeSVGFile( file ); }
			}
		}

		// If the command comes from the "Go To Structure..." item, switch
		// structures in the structure window.
		else if( command.equals( "Go to Structure..." ) ) { new ChooseStructureDialog(); }

		// If the command is to zoom the image, show a zooming dialog.
		else if( command.equals( "Zoom" ) ) { dialog.viewZoomDialog(); }
	}
	private void updateMenuUI() {
		circular.setSelected(dialog.getOption(StructureDialog.DrawOptions.CIRCULAR));
		flat.setSelected(dialog.getOption(StructureDialog.DrawOptions.FLAT));
		flipped.setSelected(dialog.getOption(StructureDialog.DrawOptions.FLIPPED));
		circled.setSelected(!dialog.getOption(StructureDialog.DrawOptions.UNCIRCLED));
	}

	/**
	 * An inner class which creates a dialog that moves to a structure.
	 *
	 * @author Richard Watson
	 */
	private class ChooseStructureDialog
		extends ValueSelectionDialog {
		private static final long serialVersionUID = 20120802;
		private final int max;
		private int prevValue = -1;
		private JButton first = new JButton("|<<"), last = new JButton(">>|"), next=new JButton(">"), prev = new JButton("<");
		private JLabel display = new JLabel("Structure 1000 of 1000"); // Structure X of Y
		private IntegerField valueField;

		/**
		 * Constructor.
		 */
		public ChooseStructureDialog() {
			// Call the superclass.
			super(  "Go to Structure" );

			max = dialog.getStructureCount();
			// Create the structure field and build the dialog with it.
			valueField = new IntegerField( "Structure Number", dialog.getStructureNumber(), 1, max );
			valueField.getDocument().addDocumentListener(new DocumentListener() {
				public void insertUpdate(final DocumentEvent e) {					valueChanged();				}
				public void removeUpdate(final DocumentEvent e) {					valueChanged();				}
				public void changedUpdate(final DocumentEvent e) {					valueChanged();				}
			});
			buildDialog( "OK", valueField);
			updateNavUI();
			updateDisplay();
		}
		private void valueChanged() {
			updateNavUI();
			int cur = getCurrent(-1);
			if (cur != -1)
				updateStructure();
		}
		private void updateNavUI() {
			Integer current = valueField.getValue();
			if (current == null) {
				first.setEnabled(true);
				last.setEnabled(true);
				next.setEnabled(false);
				prev.setEnabled(false);
			} else {
				int i = current;
				first.setEnabled(i > 1);
				last.setEnabled(i < max);
				next.setEnabled(i < max);
				prev.setEnabled(i > 1);
			}
		}

		private void setCurrent(int value) {
			if (value < 1) value = 1;
			if (value > max) value = max;
			valueField.setValue(value);
		}
		private int getCurrent(int defaultValue) {
			Integer c = valueField.getValue();
			if (c == null) return defaultValue;
			if (c < 1) return 1;
			if (c > max) return max;
			return c;
		}
		private int getCurrent() { return getCurrent(1); }
		/**
		 * Create any additional controls that need to show up on this dialog.
		 * <br><br>
		 * This method does not need to be overridden; it's only to provide an
		 * ability to add more to this dialog when the situation calls for it.
		 */
		@Override
		public JPanel createAdditionalControls() {
			JPanel pnl = new JPanel();
			pnl.add(first , BorderLayout.WEST);
			pnl.add(prev, BorderLayout.WEST);
			pnl.add(display, BorderLayout.CENTER);
			pnl.add(next, BorderLayout.EAST);
			pnl.add(last, BorderLayout.EAST);
			first.addActionListener(new ActionListener() { public void actionPerformed(final ActionEvent e) { setCurrent(1); }});
			last.addActionListener(new ActionListener() { public void actionPerformed(final ActionEvent e) { setCurrent(max); }	});
			next.addActionListener(new ActionListener() { public void actionPerformed(final ActionEvent e) { setCurrent(getCurrent(0)+1); }});
			prev.addActionListener(new ActionListener() { public void actionPerformed(final ActionEvent e) { setCurrent(getCurrent()-1); }});
			return pnl;
		}
		@Override
		public ActionListener createSelectionAction() {
			return new ActionListener() {
				public void actionPerformed( ActionEvent e ) {
					updateStructure();
					dispose();
				}
			};
		}
		private void updateStructure() {
			if (dialog.getStructureNumber() != getCurrent())
				dialog.setStructureNumber(getCurrent());
			dialog.zoomImage();
			dialog.repaint();
			updateDisplay();
		}
		private void updateDisplay() {
			display.setText(String.format("Structure %s of %s",  dialog.getStructureNumber(), max));
		}
	}
}
