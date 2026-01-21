/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.ui;

import javax.swing.*;
import java.awt.*;

/**
 * A class that creates a text field which holds a file name.
 *
 * @author Jessica S. Reuter
 */
public class FileField
	extends JTextField {
	private static final long serialVersionUID = 20120802;
	private boolean isOutputFile,isInputFile;
	private String filters;


	/**
	 * Constructor.
	 * 
	 * @param text      The name of the field.
	 * @param enabled   Whether the field is enabled.
	 */
	private FileField( String text, boolean enabled ) {
		setName( text );
		//System.out.println("FileField Role: " + this.getAccessibleContext().getAccessibleRole().toDisplayString());
		//System.out.println("FileField RoleClass: " + this.getAccessibleContext().getAccessibleRole().getClass().getTypeName());
		//System.out.println("FileField Name: " + this.getAccessibleContext().getAccessibleName());
		//System.out.println("FileField Desc: " + this.getAccessibleContext().getAccessibleDescription());
		this.getAccessibleContext().setAccessibleName("Input " + text);
		this.getAccessibleContext().setAccessibleDescription("text");
		//System.out.println("FileField Name: " + this.getAccessibleContext().getAccessibleName());
		//System.out.println("FileField Desc: " + this.getAccessibleContext().getAccessibleDescription());
		//System.out.println("FileField role-after: " + this.getAccessibleContext().getAccessibleRole().toDisplayString());

		if( !enabled ) {
			setEditable( false );
			setBackground( Color.WHITE );
		}
	}

	/**
	 * Create a disabled file input field.
	 *
	 * @return   The field.
	 */
	public static FileField createDisabled( String text ) {
		return new FileField( text, false );
	}


	public boolean isInputFile() {
		return isInputFile;
	}
	public boolean isOutputFile() {
		return isOutputFile;
	}

	/**
	 * Indicate whether or not this field is used to accept an input file (as opposed to an output file)
	 * @param value Whether or not this field is an input file.
	 * @return This FileField (for chaining)
     */
	public FileField inputFile(boolean value) { 		isInputFile = value;		return this;	}
	public FileField inputFile() { return inputFile(true); }
	public FileField inputFile(String filters) { return inputFile(true).setFilters(filters); }

	public FileField outputFile(boolean value) {		isOutputFile = value;		return this; }
	public FileField outputFile() { return outputFile(true); }
	public FileField outputFile(String filters) { return outputFile(true).setFilters(filters); }

	/**
	 * Create an enabled file input field.
	 *
	 * @return   The field.
	 */
	public static FileField createEnabled( String text ) {
		return new FileField( text, true );
	}

	/**
	 * Get the file name in the field.
	 *
	 * @return   The file name.
	 */
	public String getFile() { return getText().trim(); }

	/**
	 * Gets the file filters associated with this field.
     */
	public String getFilters() {
		return filters;
	}
	/**
	 * Sets the file filters associated with this field.
	 */
	public FileField setFilters(final String filters) {
		this.filters = filters;
		return this;
	}
}
