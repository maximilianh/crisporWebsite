import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.Vector;

public class ScrollDemo2 extends JPanel
        implements MouseListener {
    private Dimension area; //indicates area taken up by graphics
    private Vector<Rectangle> circles; //coordinates used to draw graphics
    private JPanel drawingPane;

    private final Color colors[] = {
            Color.red, Color.blue, Color.green, Color.orange,
            Color.cyan, Color.magenta, Color.darkGray, Color.yellow};
    private final int color_n = colors.length;

    public ScrollDemo2() {
        super(new BorderLayout());

        area = new Dimension(0,0);
        circles = new Vector<Rectangle>();

        //Set up the instructions.
        JLabel instructionsLeft = new JLabel(
                "Click left mouse button to place a circle.");
        JLabel instructionsRight = new JLabel(
                "Click right mouse button to clear drawing area.");
        JPanel instructionPanel = new JPanel(new GridLayout(0,1));
        instructionPanel.setFocusable(true);
        instructionPanel.add(instructionsLeft);
        instructionPanel.add(instructionsRight);

        //Set up the drawing area.
        drawingPane = new DrawingPane();
        drawingPane.setBackground(Color.white);
        drawingPane.addMouseListener(this);

        //Put the drawing area in a scroll pane.
        JScrollPane scroller = new JScrollPane(drawingPane);
        scroller.setPreferredSize(new Dimension(200,200));

        //Lay out this demo.
        add(instructionPanel, BorderLayout.PAGE_START);
        add(scroller, BorderLayout.CENTER);
    }

    /** The component inside the scroll pane. */
    public class DrawingPane extends JPanel {
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);

            Rectangle rect;
            for (int i = 0; i < circles.size(); i++) {
                rect = circles.elementAt(i);
                g.setColor(colors[(i % color_n)]);
                g.fillOval(rect.x, rect.y, rect.width, rect.height);
            }
        }
    }

    //Handle mouse events.
    public void mouseReleased(MouseEvent e) {
        final int W = 100;
        final int H = 100;
        boolean changed = false;
        if (SwingUtilities.isRightMouseButton(e)) {
            //This will clear the graphic objects.
            circles.removeAllElements();
            area.width=0;
            area.height=0;
            changed = true;
        } else {
            int x = e.getX() - W/2;
            int y = e.getY() - H/2;
            if (x < 0) x = 0;
            if (y < 0) y = 0;
            Rectangle rect = new Rectangle(x, y, W, H);
            circles.addElement(rect);
            drawingPane.scrollRectToVisible(rect);

            int this_width = (x + W + 2);
            if (this_width > area.width) {
                area.width = this_width; changed=true;
            }

            int this_height = (y + H + 2);
            if (this_height > area.height) {
                area.height = this_height; changed=true;
            }
        }
        if (changed) {
            //Update client's preferred size because
            //the area taken up by the graphics has
            //gotten larger or smaller (if cleared).
            drawingPane.setPreferredSize(area);

            //Let the scroll pane know to update itself
            //and its scrollbars.
            drawingPane.revalidate();
        }
        drawingPane.repaint();
    }
    public void mouseClicked(MouseEvent e){}
    public void mouseEntered(MouseEvent e){}
    public void mouseExited(MouseEvent e){}
    public void mousePressed(MouseEvent e){}

    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event-dispatching thread.
     */
    private static void createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame("ScrollDemo2");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Create and set up the content pane.
        JComponent newContentPane = new ScrollDemo2();
        newContentPane.setOpaque(true); //content panes must be opaque
        frame.setContentPane(newContentPane);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }

    public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                createAndShowGUI();
            }
        });
    }
}