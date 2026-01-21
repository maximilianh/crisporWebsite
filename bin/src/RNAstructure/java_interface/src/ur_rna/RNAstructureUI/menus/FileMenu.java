/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.menus;

import ur_rna.RNAstructureUI.RNAstructure;
import ur_rna.RNAstructureUI.RNAstructureBackendCalculator;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.FileFilters;
import ur_rna.RNAstructureUI.utilities.UserOptions;
import ur_rna.RNAstructureUI.windows.*;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.swing.MergeMenu;

import javax.swing.*;
import java.awt.event.ActionEvent;

/**
 * A class that creates a "File" menu.
 *
 * Most File menus are identical, but items can be added to or subtracted from
 * them in specific contexts.
 *
 * @author Jessica S. Reuter
 * @author Richard M. Watson -- significant revisions 2015-2017
 */
public class FileMenu extends MainMenu {
	private static final long serialVersionUID = 20120802;
	private final MergeMenu drawOptionMenu;

	public static final int SequenceSection = -100;
	public static final int DrawSection = -80;
	public static final int RefoldSection = -60;

	/**
	 * Constructor.
	 */
	public FileMenu() {
		super( "File" );
		this.setSubItemMergePos(SequenceSection);
		addItem( "&New Sequence", "Create a new sequence.", 'N' );
		addItem( "&Open Sequence", "Open an existing sequence.", 'O' );
		addSeparator();
		this.setSubItemMergePos(DrawSection);
		addItem( "O&ligoScreen", "Calculate thermodynamic parameters for a set of oligonucleotides." );
		addSeparator();
		addItem( "&Draw", "Draw a secondary structure." );
		addItem( "Do&t Plot", "Display the energy dot plot for a sequence that was previously folded." );
		addItem( "Dot Plot Pa&rtition Function",
			"Display base pairing probabilities for a previously calculated sequence." );
		addItem( "Dot Plot Dyn&align", "Generate a Dynalign dot plot for two sequences." );
		addItem( "Dot Plot From T&ext File", "Draw a dot plot from a text file." );
		addSeparator();
		add(drawOptionMenu = createDrawStructuresMenu());
		addSeparator();
		this.setSubItemMergePos(RefoldSection);
		addItem( "&Refold From Save File", "Refold a sequence from its save file." );
		addItem( "Refold From Dyn&align Save File", "Refold from a Dynalign calculation." );
		this.setSubItemMergePos(0);
		//addSeparator();
		this.setSubItemMergePos(Integer.MAX_VALUE);
		addSeparator();
        if (AppLog.getDefault().isDebugEnabled()) {
            addItem("Test", "Create a new sequence.", 'T');
            addSeparator();
        }
		addItem( "E&xit", "Exit the RNAstructure application." );
	}
	private MainMenu createDrawStructuresMenu() {
		MainMenu drawMenu = new MainMenu("Draw structures after calculations...");
		drawMenu.addCheckItem("Always Draw.", "After each calculation, structures will automatically be drawn without prompting.");
		drawMenu.addCheckItem("Never Draw.", "Structures will not be drawn after calculations (but output files will still be saved.");
		drawMenu.addCheckItem("Prompt to Draw after each Calculation.", "After each calculation, a prompt will be shown asking whether or not to draw structures.");
		writeDrawPref(drawMenu, null); // updates the check state to match current preference
		return drawMenu;
	}

	@Override
	protected void onMenuAction(String command, final ActionEvent ev) {
        if (command.equals("Test")) {
            Dialogs.showMessage("DATAPATH=" + RNAstructureBackendCalculator.getEnvVar("DATAPATH"));
			//Dialogs.showMessage("drawOptionMenu: " + drawOptionMenu.getPopupMenu().getComponents().length);
            return;
        }

		// If the command is a drawing command, prepare a drawing window.
		if( command.startsWith( "D" ) ) {
			// Show an energy dot plot window, if possible.
			if( command.equals( "Dot Plot" ) ) {
				String file = StandardFileChooser.getOpenName(FileFilters.FoldingSav);
				if( file != null ) {
					DrawingWindow window = new DrawingWindow( file );
					if(!window.isError()) { window.showWindow(); }
				}
			}

			// Show a Dynalign dot plot window, if possible.
			else if( command.equals( "Dot Plot Dynalign" ) ) {
				String file = StandardFileChooser.getOpenName(FileFilters.DynalignSav);
				if( file != null ) {
					DrawingWindow window1 = new DrawingWindow( file, 1 );
					DrawingWindow window2 = new DrawingWindow( file, 2 );
					boolean validWindows =
						(!window1.isError()) &&
						(!window2.isError());
					if( validWindows ) {
						window1.showWindow();
						window2.showWindow();
					}
				}
			}

			// Show a text dot plot window, if possible.
			else if( command.equals( "Dot Plot From Text File" ) ) {
				String file = StandardFileChooser.getOpenName(FileFilters.DotPlot);
				if( file != null ) {
					DrawingWindow window = new DrawingWindow( file );
					if(!window.isError()) { window.showWindow(); }
				}
			}

			// Show a probability dot plot window, if possible.
			else if( command.equals( "Dot Plot Partition Function" ) ) {
				String file = StandardFileChooser.getOpenName(FileFilters.PartitionSav);
				if( file != null ) {
					DrawingWindow window = new DrawingWindow( file );
					if(!window.isError()) { window.showWindow(); }
				}
			}

			// Initialize a structure window, if possible.
			else if( command.equals( "Draw" ) ) {
				String file = StandardFileChooser.getOpenName(FileFilters.CT);
				if( file != null ) {
					DrawingWindow window = new DrawingWindow( file );
					if(!window.isError()) { window.showWindow(); }
				}
			}
		}

		else if( command.equals( "Always Draw." ) ) {
			writeDrawPref(null, UserOptions.DRAW_ALWAYS);
		} else if( command.equals( "Never Draw." ) ) {
			writeDrawPref(null, UserOptions.DRAW_NEVER);
		}else if( command.startsWith( "Prompt to Draw" ) ) {
			writeDrawPref(null, UserOptions.DRAW_PROMPT);
		}

		// If the command is to exit the application, do so.
		else if( command.equals( "Exit" ) ) { System.exit( 0 ); }

		// If the command is to build a new sequence, show a blank sequence
		// display window.
		else if( command.equals( "New Sequence" ) ) {
			new SequenceDisplayWindow().showWindow();
		}

		// If the command is to run OligoScreen, show an OligoScreen window.
		else if( command.equals( "OligoScreen" ) ) {
			new OligoScreenWindow().showWindow();
		}

		// If the command is to open an existing sequence, show a filled
		// sequence display window.
		else if( command.equals( "Open Sequence" ) ) {
			String file = StandardFileChooser.getOpenName(FileFilters.SequenceExtended);
			if( file != null ) {
				new SequenceDisplayWindow( file ).showWindow();
			}
		}

		// If the command is to refold from a Dynalign save file, open a
		// Dynalign refolding window.
		else if( command.equals( "Refold From Dynalign Save File" ) ) {
			new DynalignRefoldWindow().showWindow();
		}

		// If the command is to refold from a folding save file, open a
		// refolding window.
		else if( command.equals( "Refold From Save File" ) ) {
			new RefoldWindow().showWindow();
		}
	}

	private void writeDrawPref(MergeMenu drawMenu, String value) {
		if (value == null)
			value = RNAstructure.options.drawStructures;
		else {
			RNAstructure.options.drawStructures = value;
			RNAstructure.options.save();
		}
		if (drawMenu == null) drawMenu = drawOptionMenu;
		java.util.List<? extends JMenuItem> items = drawMenu.getMenuItems();
		items.get(0).setSelected(UserOptions.DRAW_ALWAYS.equals(value));
		items.get(1).setSelected(UserOptions.DRAW_NEVER.equals(value));
		items.get(2).setSelected(UserOptions.DRAW_PROMPT.equals(value));
	}
//	private void select(JComponent parent, int index, boolean value) {
//		JMenuItem mi = item(parent, index);
//		if (mi == null ) throw new NullPointerException("Missing Menu Item " + index);
//		mi.setSelected(value);
//	}
//	private JMenuItem item(JComponent parent, int index) {
//		int count = 0;
//		for(Component c : parent.getComponents())
//			if (c instanceof JMenuItem && (count++ == index))
//				return (JMenuItem)c;
//		return null;
//	}
}
