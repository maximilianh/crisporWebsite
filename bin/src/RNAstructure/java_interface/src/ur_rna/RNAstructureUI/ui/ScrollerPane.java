/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */
package ur_rna.RNAstructureUI.ui;

import ur_rna.Utilities.OSInfo;

import javax.swing.*;
import java.awt.*;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
/**
 * A class that creates a scroll pane that can be controlled via the keyboard.
 *
 * @author Jessica S. Reuter
 */
public class ScrollerPane
	extends JScrollPane {
	private static final long serialVersionUID = 20120802;
	/**
	 * A generic constant symbolizing that a scroll bar must always be used.
	 */
	public static final int ALWAYS = 1;
	/**
	 * A generic constant symbolizing that a scroll bar must never be used.
	 */
	public static final int NEVER = 2;
	/**
	 * Constructor.
	 *
	 * @param comp     The component to put in the scroll pane.
	 * @param width    The width of the scroll pane.
	 * @param height   The height of the scroll pane.
	 */
	public ScrollerPane( Component comp, int width, int height ) {
		super( comp );
		// Set the scroll pane size and set it focusable.
		setSizes( width, height );
		setFocusable( true );
		// Bind actions to the scroll pane so the user can scroll with arrow
		// keys. If the OS is a type of Mac, add a key listener as well as the
		// key bindings.
		JScrollBar horizontal = getHorizontalScrollBar();
		InputMap horizontalMap =
			horizontal.getInputMap( JComponent.WHEN_IN_FOCUSED_WINDOW );
		horizontalMap.put(
			KeyStroke.getKeyStroke( "RIGHT" ), "positiveUnitIncrement" );
		horizontalMap.put(
			KeyStroke.getKeyStroke( "LEFT" ), "negativeUnitIncrement" );
		JScrollBar vertical = getVerticalScrollBar();
		InputMap verticalMap =
			vertical.getInputMap( JComponent.WHEN_IN_FOCUSED_WINDOW );
		verticalMap.put(
			KeyStroke.getKeyStroke( "DOWN" ), "positiveUnitIncrement" );
		verticalMap.put(
			KeyStroke.getKeyStroke( "UP" ), "negativeUnitIncrement" );
		if( OSInfo.isMac() ) {
			final JScrollPane source = this;
			addKeyListener( new KeyAdapter() {
				@Override
				public void keyPressed( KeyEvent e ) {
					int code = e.getKeyCode();
					JScrollBar horizontalBar = source.getHorizontalScrollBar();
					JScrollBar verticalBar = source.getVerticalScrollBar();
					int horizontal = horizontalBar.getValue();
					int vertical = verticalBar.getValue();
					if( code == KeyEvent.VK_UP ) {
						verticalBar.setValue( vertical - 1 );
					} else if( code == KeyEvent.VK_DOWN ) {
						verticalBar.setValue( vertical + 1 );
					} else if( code == KeyEvent.VK_LEFT ) {
						horizontalBar.setValue( horizontal - 1 );
					} else if( code == KeyEvent.VK_RIGHT ) {
						horizontalBar.setValue( horizontal + 1 );
					}
				}
			});
		}
	}
	/**
	 * Set the scroll bar policies for this scroll pane.
	 *
	 * @param horizontal   The horizontal policy.
	 * @param vertical     The vertical policy.
	 */
	public void setBarPolicies( int horizontal, int vertical ) {
		// Set the horizontal scroll bar policy.
		int horizontalPolicy =
			( horizontal == ALWAYS ) ?
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS :
			( horizontal == NEVER ) ?
				JScrollPane.HORIZONTAL_SCROLLBAR_NEVER :
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED;
		setHorizontalScrollBarPolicy( horizontalPolicy );
		// Set the vertical scroll bar policy.
		int verticalPolicy =
			( vertical == ALWAYS ) ?
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS :
			( vertical == NEVER ) ?
				JScrollPane.VERTICAL_SCROLLBAR_NEVER :
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED;
		setVerticalScrollBarPolicy( verticalPolicy );
	}
	/**
	 * Set the size of the scroller pane by setting preferred, minimum, and
	 * maximum sizes all at once.
	 *
	 * @param width    The width of the scroll pane.
	 * @param height   The height of the scroll pane.
	 */
	public void setSizes( int width, int height ) {
		Dimension size = new Dimension( width, height );
		setPreferredSize( size );
		setMinimumSize( size );
		setMaximumSize( size );
	}
}
