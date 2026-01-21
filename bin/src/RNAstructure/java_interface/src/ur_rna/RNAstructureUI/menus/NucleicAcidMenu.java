/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.menus;

import ur_rna.RNAstructureUI.windows.*;
import ur_rna.Utilities.AppLog;

import javax.swing.*;
import java.awt.event.ActionEvent;

/**
 * A class that creates a menu which holds commands specific for a nucleic
 * acid type.
 * <br><br>
 * Most File menus are identical, but items can be added to or subtracted from
 * them in specific contexts.
 *
 * @author Jessica S. Reuter
 */
public class NucleicAcidMenu
	extends MainMenu {
	private static final long serialVersionUID = 20120802;
	private final String acidName;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	private NucleicAcidMenu( String acid ) {
		super( acid );
		acidName = acid;
		addAcidMenuItem("FoldSingleStrand", "Fold {acid} &Single Strand", "Predict a secondary structure for a single strand using {acid} thermodynamics.");
		addAcidMenuItem("FoldBimolecular", "Fold {acid} &Bimolecular", "Predict a secondary structure for two strands using {acid} thermodynamics.");
		addAcidMenuItem("AccessFold", "{acid} &AccessFold", "Predict a secondary structure for two strands using the AccessFold algorithm.");
		addAcidMenuItem("DuplexFold", "{acid} &DuplexFold", "Predict a secondary structure for two strands using the DuplexFold algorithm.");
		addSeparator();
		addAcidMenuItem("Partition", "&Partition Function {acid}", "Predict base pairing probabilities for all {acid} pairs.");
		addAcidMenuItem("PartitionBimolecular", "Partition Function {acid} Bi&molecular", "Predict base pairing probabilities for all {acid} pairs for two strands.");
		addSeparator();
		addAcidMenuItem("AllSub", "Generate &All Suboptimal {acid} Structures", "Generate all {acid} structures within a given increment of the lowest free energy structure.");
		addAcidMenuItem("Stochastic", "Stoc&hastic {acid} Sampling", "Predict {acid} structures using stochastic sampling.");
		addAcidMenuItem("MaxExpect", "Max&Expect: Predict {acid} MEA Structure", "Predict the maximum expected accuracy {acid} structure.");
		addAcidMenuItem("ProbKnot", "Prob&Knot: Predict {acid} Structures Including Pseudoknots", "Predict pseudoknots in a {acid} structure.");
		addAcidMenuItem("Efn2", "E&fn2 {acid}", "Calculate the free energy of a {acid} structure.");
		addSeparator();
		addAcidMenuItem("Dynalign", "{acid} &Dynalign", "Find a common secondary structure for two {acid} sequences.");
		addAcidMenuItem("Multilign", "{acid} M&ultilign", "Predict a secondary structure common to multiple {acid} sequences.");
		addSeparator();
		addAcidMenuItem("BreakPseudoknots", "&Break {acid} Pseudoknots", "Break pseudoknots in a structure, leaving the lowest free energy pseudoknot-free structure.");

		if (acid.equals("RNA")) {
			// RNA-only items
			addAcidMenuItem("Oligo&Walk", "{acid} OligoWalk", "RNA OligoWalk calculation.");
			addAcidMenuItem("&TurboFold", "{acid} TurboFold", "Find common structures for multiple sequences using base pair probabilities.");
		}
	}

	private void addAcidMenuItem(String action, String text, String tooltip) {
		super.addItem(text.replace("{acid}", acidName),
				tooltip.replace("{acid}", acidName),
				(KeyStroke)null, action.replace("{acid}", acidName));
	}

	/**
	 * Create a DNA menu.
	 *
	 * @return   The DNA menu.
	 */
	public static NucleicAcidMenu createDNA() {
		return new NucleicAcidMenu( "DNA" );
	}

	/**
	 * Create a RNA menu.
	 *
	 * @return   The RNA menu.
	 */
	public static NucleicAcidMenu createRNA() { return new NucleicAcidMenu( "RNA" ); }

	@Override
	protected void onMenuAction(String command, final ActionEvent ev) {
		String type = this.acidName;
		switch (command) {
			case "FoldSingleStrand":
				new FoldSingleWindow(type).showWindow();
				break;
			case "FoldBimolecular":
				new FoldDoubleWindow(type).showWindow();
				break;
			case "AccessFold":
				new AccessFoldWindow(type).showWindow();
				break;
			case "DuplexFold":
				new DuplexFoldWindow(type).showWindow();
				break;
			case "Partition":
				new PartitionSingleWindow(type).showWindow();
				break;
			case "PartitionBimolecular":
				new PartitionDoubleWindow(type).showWindow();
				break;
			case "AllSub":
				new FoldSuboptimalWindow(type).showWindow();
				break;
			case "Stochastic":
				new StochasticWindow(type).showWindow();
				break;
			case "MaxExpect":
				new MaxExpectWindow(type).showWindow();
				break;
			case "ProbKnot":
				new ProbKnotWindow(type).showWindow();
				break;
			case "Efn2":
				new Efn2Window(type).showWindow();
				break;
			case "Dynalign":
				new DynalignFoldWindow(type).showWindow();
				break;
			case "Multilign":
				new MultilignWindow(type).showWindow();
				break;
			case "BreakPseudoknots":
				new PseudoknotWindow(type).showWindow();
				break;
			case "OligoWalk":
				new OligoWalkWindow().showWindow();
				break;
			case "TurboFold":
				new TurboFoldWindow().showWindow();
				break;
			default:
				AppLog.getDefault().warn("Unknown menu command in " + this.getClass().getName() + "#onMenuAction");
				break;
		}
	}
}
