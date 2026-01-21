package ur_rna.Utilities.swing;

import java.awt.*;
import java.awt.event.*;

/**
 * Funnels 10 input events (from MouseListener, MouseMotionListener, and KeyListener) into
 * two events -- one for key input and one for mouse input.
 */
public abstract class InputAdapter {
    public enum InputType {
        KeyDown, KeyUp, KeyTyped,
        MouseDown, MouseUp, Drag, Move, Wheel, Enter, Exit, Click
    }

    //public abstract void listen(Component c);
    //protected abstract void onInput(InputType type, InputEvent e);

    public static abstract class Keyboard implements KeyListener { //extends InputAdapter
        /**
         * Invoked when a key has been typed.
         * See the class description for {@link KeyEvent} for a definition of
         * a key typed event.
         *
         * @param e
         */
        @Override
        public void keyTyped(final KeyEvent e) {
            onKeyInput(InputType.KeyTyped, e);
        }
        /**
         * Invoked when a key has been pressed.
         * See the class description for {@link KeyEvent} for a definition of
         * a key pressed event.
         *
         * @param e
         */
        @Override
        public void keyPressed(final KeyEvent e) {
            onKeyInput(InputType.KeyDown, e);
        }
        /**
         * Invoked when a key has been released.
         * See the class description for {@link KeyEvent} for a definition of
         * a key released event.
         *
         * @param e
         */
        @Override
        public void keyReleased(final KeyEvent e) {
            onKeyInput(InputType.KeyUp, e);
        }

        //@Override
        public void listen(Component c) {
            c.addKeyListener(this);
        }
        protected abstract void onKeyInput(InputType type, KeyEvent e);
    }

    public static abstract class Mouse implements MouseListener, MouseMotionListener, MouseWheelListener {
        /**
         * Invoked when the mouse button has been clicked (pressed
         * and released) on a component.
         *
         * @param e
         */
        @Override
        public void mouseClicked(final MouseEvent e) {
            onMouseInput(InputType.Click, e);
        }
        /**
         * Invoked when a mouse button has been pressed on a component.
         *
         * @param e
         */
        @Override
        public void mousePressed(final MouseEvent e) {
            onMouseInput(InputType.MouseDown, e);
        }
        /**
         * Invoked when a mouse button has been released on a component.
         *
         * @param e
         */
        @Override
        public void mouseReleased(final MouseEvent e) {
            onMouseInput(InputType.MouseUp, e);
        }
        /**
         * Invoked when the mouse enters a component.
         *
         * @param e
         */
        @Override
        public void mouseEntered(final MouseEvent e) { onMouseInput(InputType.Enter, e);}
        /**
         * Invoked when the mouse exits a component.
         *
         * @param e
         */
        @Override
        public void mouseExited(final MouseEvent e) {onMouseInput(InputType.Exit, e);}
        /**
         * Invoked when a mouse button is pressed on a component and then
         * dragged.  <code>MOUSE_DRAGGED</code> events will continue to be
         * delivered to the component where the drag originated until the
         * mouse button is released (regardless of whether the mouse position
         * is within the bounds of the component).
         * <p>
         * Due to platform-dependent Drag&amp;Drop implementations,
         * <code>MOUSE_DRAGGED</code> events may not be delivered during a native
         * Drag&amp;Drop operation.
         *
         * @param e
         */
        @Override
        public void mouseDragged(final MouseEvent e) {
            onMouseInput(InputType.Drag, e);
        }
        /**
         * Invoked when the mouse cursor has been moved onto a component
         * but no buttons have been pushed.
         *
         * @param e
         */
        @Override
        public void mouseMoved(final MouseEvent e) {
            onMouseInput(InputType.Move, e);
        }

        /**
         * Invoked when the mouse wheel is rotated.
         *
         * @param e
         * @see MouseWheelEvent
         */
        @Override
        public void mouseWheelMoved(final MouseWheelEvent e) {
            onMouseInput(InputType.Wheel, e);
        }

        //@Override
        public void listen(Component c) { listen(c, true, true, true);  }
        public void remove(Component c) {
            c.removeMouseListener(this);
            c.removeMouseMotionListener(this);
            c.removeMouseWheelListener(this);
        }
        public void listen(Component c, boolean includeMouse, boolean includeMotion, boolean includeWheel) {
            if (includeMouse) c.addMouseListener(this);
            if (includeMotion) c.addMouseMotionListener(this);
            if (includeWheel) c.addMouseWheelListener(this);
        }

        /**
         * Override this to respond to mouse input.
         * @param type The type of mouse input (e.g. click, move, etc)
         * @param e A MouseEvent with additional information about the event.
         */
        protected abstract void onMouseInput(InputType type, MouseEvent e);
    }
}
