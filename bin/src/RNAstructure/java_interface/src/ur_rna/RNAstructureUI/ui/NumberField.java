/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.ui;

import javax.swing.*;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

/**
 * A class that creates a text field whose contents may be bounded by a
 * minimum or maximum value.
 *
 * @author Jessica S. Reuter
 */
public abstract class NumberField
	extends JTextField {
	private static final long serialVersionUID = 20120802;

	/**
	 * The initial default value of the field.
	 */
	protected Number initial;

	/**
	 * The maximum value of the field.
	 */
	protected Number maximum;

	/**
	 * The minimum value of the field.
	 */
	protected Number minimum;

	/**
	 * Constructor.
	 * 
	 * @param text    The name of the field.
	 * @param value   The initial numeric value of the field.
	 * @param min     The minimum value of the field.
	 * @param max     The maximum value of the field.
	 */
	protected NumberField(
		String text, Number value, Number min, Number max ) {

		// Set the field name and values.
		setName( text );
		resetField( value, min, max );

		// Add a key listener so the field can only accept valid characters.
		addKeyListener( new KeyAdapter() {

			/**
			 * Check to see if the next character is valid, and if not,
			 * consume it.
			 *
			 * @param e   The key event.
			 */
			private void checkKey( KeyEvent e ) {

				// Get the next character.
				char key = e.getKeyChar();

				// If the next character is a digit, allow it.
				if( Character.isDigit( key ) ) { return; }

				// If the next character is a backspace character, allow it.
				if( key == KeyEvent.VK_BACK_SPACE ) { return; }

				// If the character is a negative sign, and a negative sign is
				// acceptable right now, allow it.
				boolean negativeSign = ( key == '-' );
				boolean isEmptyText = getText().trim().equals( "" );
				boolean negativeBound = ( minimum.doubleValue() < 0.0 );
				if( negativeSign && isEmptyText && negativeBound ) { return; }

				// Do any extra checking specific to a particular field type,
				// and if the character was valid by that check, allow it.
				boolean valid = doAdditionalCheck( key );
				if( valid ) { return; }

				// If none of the previous validity checks worked, consume the
				// character so it isn't outputted.
				e.consume();
			}

			@Override
			public void keyPressed( KeyEvent e ) { checkKey( e ); }

			@Override
			public void keyTyped( KeyEvent e ) { checkKey( e ); }
		});

		// Add a focus listener so the field will check itself as soon as the
		// focus is lost.
		addFocusListener( new FocusAdapter() {

			@Override
			public void focusLost( FocusEvent e ) {
				boolean ok = checkValue();
				if( ok == false ) { setText( initial.toString() ); }
			}
		});
	}

	/**
	 * Check the value of a field to make sure it's valid.
	 *
	 * @return   True if the field value is valid, false if not.
	 */
	protected abstract boolean checkValue();

	/**
	 * Do any additional necessary checks on a key character to make sure it's
	 * valid to be output.
	 *
	 * @param key   The key character.
	 */
	protected boolean doAdditionalCheck( char key ) { return false; }

	/**
	 * Get the text field's numeric value.
	 *
	 * @return   The numeric value.
	 */
	public abstract Number getValue();

	/**
	 * Reset a field value, if it's within bounds.
	 *
	 * @param value   The initial numeric value of the field.
	 */
	public void resetField( Number value ) {

		// If the range of numbers given isn't valid, return.
		boolean goodMin = ( minimum.doubleValue() <= value.doubleValue() );
		boolean goodMax = ( maximum.doubleValue() >= value.doubleValue() );
		if( !( goodMin && goodMax ) ) { return; }

		// Set the field value.
		setText( value.toString() );
		initial = value;
	}

	/**
	 * Reset a field value and bounds.
	 *
	 * @param value   The initial numeric value of the field.
	 * @param min     The minimum value of the field.
	 * @param max     The maximum value of the field.
	 */
	public void resetField( Number value, Number min, Number max ) {

		// If the range of numbers given isn't valid, return.
		boolean goodMin = ( min.doubleValue() <= value.doubleValue() );
		boolean goodMax = ( max.doubleValue() >= value.doubleValue() );
		if( !( goodMin && goodMax ) ) { return; }

		// Set the field value and bounds.
		setText( value.toString() );
		initial = value;
		minimum = min;
		maximum = max;
	}

	/**
	 * An inner class that handles a field which contains a double value.
	 *
	 * @author Jessica S. Reuter
	 */
	public static class DoubleField
		extends NumberField {
		private static final long serialVersionUID = 20120802;

		/**
		 * Unbounded Constructor.
		 * 
		 * @param text    The name of the field.
		 * @param value   The initial numeric value of the field.
		 */
		public DoubleField( String text, Number value ) {
			super( text, value, Double.MAX_VALUE * -1, Double.MAX_VALUE );
		}

		/**
		 * Minimum Bounds Constructor.
		 * 
		 * @param text    The name of the field.
		 * @param value   The initial numeric value of the field.
		 * @param min     The minimum value of the field.
		 */
		public DoubleField( String text, Number value, Number min ) {
			super( text, value, min, Double.MAX_VALUE );
		}

		/**
		 * Minimum and Maximum Bounds Constructor.
		 * 
		 * @param text    The name of the field.
		 * @param value   The initial numeric value of the field.
		 * @param min     The minimum value of the field.
		 * @param max     The maximum value of the field.
		 */
		public DoubleField(
			String text, Number value, Number min, Number max ) {
			super( text, value, min, max );
		}

		@Override
		protected boolean checkValue() {
			try {
				Double value = Double.parseDouble( getText() );
				if( value < minimum.doubleValue() ) { return false; }
				if( value > maximum.doubleValue() ) { return false; }
			} catch( Exception e ) { return false; }
			return true;
		}

		@Override
		protected boolean doAdditionalCheck( char key ) {
			boolean decimal = ( key == '.' );
			boolean noPreviousDecimal = ( getText().contains( "." ) == false );
			return ( decimal && noPreviousDecimal );
		}

		@Override
		public Double getValue() {
			Double value = null;
			try { value = Double.parseDouble( getText() ); }
			catch( NumberFormatException e ) { /* Will never be thrown */ }
			return value;
		}
	}

	/**
	 * An inner class that handles a field which contains a float value.
	 *
	 * @author Jessica S. Reuter
	 */
	public static class FloatField
		extends NumberField {
		private static final long serialVersionUID = 20120802;

		/**
		 * Unbounded Constructor.
		 * 
		 * @param text    The name of the field.
		 * @param value   The initial numeric value of the field.
		 */
		public FloatField( String text, Number value ) {
			super( text, value, Float.MAX_VALUE * -1, Float.MAX_VALUE );
		}

		/**
		 * Minimum Bounds Constructor.
		 * 
		 * @param text    The name of the field.
		 * @param value   The initial numeric value of the field.
		 * @param min     The minimum value of the field.
		 */
		public FloatField( String text, Number value, Number min ) {
			super( text, value, min, Float.MAX_VALUE );
		}

		/**
		 * Minimum and Maximum Bounds Constructor.
		 * 
		 * @param text    The name of the field.
		 * @param value   The initial numeric value of the field.
		 * @param min     The minimum value of the field.
		 * @param max     The maximum value of the field.
		 */
		public FloatField(
			String text, Number value, Number min, Number max ) {
			super( text, value, min, max );
		}

		@Override
		protected boolean checkValue() {
			try {
				Float value = Float.parseFloat( getText() );
				if( value < minimum.floatValue() ) { return false; }
				if( value > maximum.floatValue() ) { return false; }
			} catch( Exception e ) { return false; }
			return true;
		}

		@Override
		protected boolean doAdditionalCheck( char key ) {
			boolean decimal = ( key == '.' );
			boolean noPreviousDecimal = ( getText().contains( "." ) == false );
			return ( decimal && noPreviousDecimal );
		}

		@Override
		public Float getValue() {
			Float value = null;
			try { value = Float.parseFloat( getText() ); }
			catch( NumberFormatException e ) { /* Will never be thrown */ }
			return value;
		}
	}

	/**
	 * An inner class that handles a field which contains an integer value.
	 *
	 * @author Jessica S. Reuter
	 */
	public static class IntegerField
		extends NumberField {
		private static final long serialVersionUID = 20120802;

		/**
		 * Unbounded Constructor.
		 * 
		 * @param text    The name of the field.
		 * @param value   The initial numeric value of the field.
		 */
		public IntegerField( String text, Integer value ) {
			super( text, value, Integer.MAX_VALUE * -1, Integer.MAX_VALUE );
		}

		/**
		 * Minimum Bounds Constructor.
		 * 
		 * @param text    The name of the field.
		 * @param value   The initial numeric value of the field.
		 * @param min     The minimum value of the field.
		 */
		public IntegerField( String text, Integer value, Integer min ) {
			super( text, value, min, Integer.MAX_VALUE );
		}

		/**
		 * Minimum and Maximum Bounds Constructor.
		 * 
		 * @param text    The name of the field.
		 * @param value   The initial numeric value of the field.
		 * @param min     The minimum value of the field.
		 * @param max     The maximum value of the field.
		 */
		public IntegerField(
			String text, Integer value, Integer min, Integer max ) {
			super( text, value, min, max );
		}

		@Override
		protected boolean checkValue() {
			try {
				Integer value = Integer.parseInt( getText() );
				if( value < minimum.intValue() ) { return false; }
				if( value > maximum.intValue() ) { return false; }
			} catch( Exception e ) { return false; }
			return true;
		}

		@Override
		public Integer getValue() {
			Integer value = null;
			try { value = Integer.parseInt( getText() ); }
			catch( NumberFormatException e ) { /* Will never be thrown */ }
			return value;
		}
		public void setValue(final int value) {
			setText(""+value);
		}
	}
}
