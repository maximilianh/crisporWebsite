/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.ui;

import javax.swing.*;
import javax.swing.border.Border;
import java.io.Serializable;

/**
 * A class that places a border type on a specified component type.
 *
 * @author Jessica S. Reuter
 */
public class BorderBuilder
	implements Serializable {
	private static final long serialVersionUID = 20120802;

	/**
	 * Create an empty border on a component.
	 *
	 * @param type   The border type.
	 * @param gap    The gap in pixels.
	 * @param comp   The component on which to create the border.
	 */
	private void makeEmptyBorder( String type, int gap, JComponent comp ) {

		// Check to see which type of border is being made.
		boolean isEqual = ( type.equals( "Equal" ) );
		boolean isBottom = ( type.equals( "Bottom" ) ) || isEqual;
		boolean isLeft = ( type.equals( "Left" ) ) || isEqual;
		boolean isRight = ( type.equals( "Right" ) ) || isEqual;
		boolean isTop = ( type.equals( "Top" ) ) || isEqual;

		// Add padding to the proper sides of the border.
		int bottom = ( isBottom ) ? gap : 0;
		int left = ( isLeft ) ? gap : 0;
		int top = ( isTop ) ? gap : 0;
		int right = ( isRight ) ? gap : 0;

		// Create the border on the component.
		Border border =
			BorderFactory.createEmptyBorder( top, left, bottom, right );
		comp.setBorder( border );
	}

	/**
	 * Create an empty bottom border on a component.
	 *
	 * @param gap    The gap in pixels.
	 * @param comp   The component on which to create the border.
	 */
	public void makeBottomBorder( int gap, JComponent comp ) {
		makeEmptyBorder( "Bottom", gap, comp );
	}

	/**
	 * Create an empty equidistant border on a component.
	 *
	 * @param gap    The gap in pixels.
	 * @param comp   The component on which to create the border.
	 */
	public void makeEqualBorder( int gap, JComponent comp ) {
		makeEmptyBorder( "Equal", gap, comp );
	}

	/**
	 * Create an empty left border on a component.
	 *
	 * @param gap    The gap in pixels.
	 * @param comp   The component on which to create the border.
	 */
	public void makeLeftBorder( int gap, JComponent comp ) {
		makeEmptyBorder( "Left", gap, comp );
	}

	/**
	 * Create an empty right border on a component.
	 *
	 * @param gap    The gap in pixels.
	 * @param comp   The component on which to create the border.
	 */
	public void makeRightBorder( int gap, JComponent comp ) {
		makeEmptyBorder( "Right", gap, comp );
	}

	/**
	 * Create an empty top border on a component.
	 *
	 * @param gap    The gap in pixels.
	 * @param comp   The component on which to create the border.
	 */
	public void makeTopBorder( int gap, JComponent comp ) {
		makeEmptyBorder( "Top", gap, comp );
	}

	/**
	 * Create a padded titled border on a component.
	 * @param title   The title of the border.
	 * @param comp    The component on which to create the border.
	 */
	public void makeTitledBorder( String title, JComponent comp ) {
		Border border = BorderFactory.createCompoundBorder(
			BorderFactory.createEmptyBorder( 2, 2, 2, 2 ),
			BorderFactory.createTitledBorder( title ) );
		comp.setBorder( border );
	}
}
