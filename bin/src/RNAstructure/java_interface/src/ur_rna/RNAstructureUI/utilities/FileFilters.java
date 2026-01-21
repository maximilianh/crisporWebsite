/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.utilities;

/**
 * A list of constants for common file extensions.
 *
 * Note that most file-specs end with a bar character. This is meant to allow easy concatenation of filters.
 *
 * @author Richard M. Watson
 */
public class FileFilters {
	private FileFilters() {} //class is not creatable

	/*	 * Select an alignment file. */
	public static String Alignment="Alignment|ali|";

	/*	 * Select a dot bracket file. */
	public static String Bracket="Bracket|bracket|";

	/*	 * Select a constraints file. */
	public static String Constraints="Constraint|con|";

	/*	 * Select a CT file. */
	public static String CT="CT|ct|";

	/*	 * Select a dot plot file. */
	public static String DotPlot="Dot Plot|dp|";

	/*	 * Select a Dynalign save file. */
	public static String DynalignSav ="Dynalign Save|dsv|";

	/*	 * Select a folding save file. */
	public static String FoldingSav ="Folding Save|sav;save|";

	/*	 * Select a structure helix file. */
	public static String Helix="Helix|txt|";

	/*	 * Select an oligo list file. */
	public static String List="List|lis|";

	/*	 * Select a general OUT file. */
	public static String OUT="OUT|out|";

	/*	 * Select a partition function save file. */
	public static String PartitionSav ="Partition Function Save|pfs|";

	/*	 * Select a Postscript image file. */
	public static String Postscript="Postscript|ps|";

	/*	 * Select an oligo report file. */
	public static String Report="Report|rep|";

	/*	 * Select a sequence file. */
	public static String Sequence="FASTA|fasta|FA|fa|Sequence|seq|txt";

	/*	 * Select a multi-sequence FASTA file. */
	public static String MultiSequence="FASTA|fasta|FA|fa";

	/**
	 * Select a sequence file in its extended context, with more than its
	 * usual possible file types.
	 */
	public static String SequenceExtended =Sequence+"Genbank|gen|Text|txt|";

	/**
	 * Select a SHAPE file.
	 */
	public static String SHAPE = "SHAPE|shape|";

	/**
	 * Select an SVG file.
	 */
	public static String SVG  = "SVG - Scalable Vector Graphics|svg|";
}
