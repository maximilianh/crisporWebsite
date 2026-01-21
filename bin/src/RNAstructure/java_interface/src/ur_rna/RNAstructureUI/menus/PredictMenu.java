/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.menus;

import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.windows.*;

import java.awt.event.ActionEvent;

/**
 * A class that creates a "File" menu.
 * <br><br>
 * Most File menus are identical, but items can be added to or subtracted from
 * them in specific contexts.
 *
 * @author Richard M. Watson
 */
public class PredictMenu
	extends MainMenu {
	private static final long serialVersionUID = 20160121;

	/**
	 * Constructor.
	 *
	 */
	public PredictMenu() {
		super( "Structure &Prediction" );
		addItem( "&Single Sequence", "Predict a Secondary Structure for a single RNA or DNA sequence.", 'P', "fold");
		addItem( "Structure &Common to Two Sequences", "Align two sequences and predict a secondary structure common to both.", 'A', "align");
		addItem( "Structure Common to &Multiple Sequences", "Align three or more sequences and predict a secondary structure common to all.", 'M', "multi-align");
		addItem( "&Bimolecular Structure", "Folds two sequences, either RNA or DNA, into their lowest hybrid free energy conformation.", 'B', "bifold");
	}

	@Override
	protected void onMenuAction(String command, final ActionEvent ev) {
		switch(command) {
			case "fold": new PredictSingleStructureWindow().showWindow(); break;
			case "align": new PredictTwoStructureWindow().showWindow(); break;
			case "multi-align": new PredictMultiStructureWindow().showWindow(); break;
			case "bifold": new PredictBimolStructureWindow().showWindow(); break;
			default:
				Dialogs.showMessage("Unknown menu action: " + command);
				break;
		}
	}
}
