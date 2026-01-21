package ur_rna.Utilities.swing;

import ur_rna.Utilities.annotation.Nullable;

import javax.swing.*;
import java.awt.*;

/**
 * Common "MessageBox" dialogs.
 */
public class Dialogs {
    public static String programTitle;
    public static Component owningComponent;

    private static Font defaultFont;
    private static Font monospaceFont;


//    private static int ERROR=JRootPane.ERROR_DIALOG,
//            QUESTION=JRootPane.QUESTION_DIALOG,
//            WARNING=JRootPane.WARNING_DIALOG,
//            INFO = JRootPane.INFORMATION_DIALOG,
//            PLAIN = JRootPane.PLAIN_DIALOG;

    public static int ERROR=JOptionPane.ERROR_MESSAGE,
            QUESTION=JOptionPane.QUESTION_MESSAGE,
            WARNING=JOptionPane.WARNING_MESSAGE,
            INFO = JOptionPane.INFORMATION_MESSAGE,
            PLAIN = JOptionPane.PLAIN_MESSAGE;

    public static int PROMPT_YES_NO = JOptionPane.YES_NO_OPTION,
            PROMPT_YES_NO_CANCEL = JOptionPane.YES_NO_CANCEL_OPTION,
            PROMPT_OK_CANCEL = JOptionPane.OK_CANCEL_OPTION,
            PROMPT_NONE = JOptionPane.DEFAULT_OPTION;;

    public static void showWarning(String message) { showMessage(message, programTitle, WARNING); }
    public static void showWarning(String message, String title) { showMessage(message, title, WARNING); }
    public static void showInfo(String message) { showMessage(message, programTitle, INFO); }
    public static void showInfo(String message, String title) { showMessage(message, title, INFO); }
    public static void showMessage(String message, String title, int dialogType) { showMessage(message, title, dialogType, false); }
    public static void showMessage(String message, String title, int dialogType, boolean monospace) { showMessage(message, title, dialogType, monospace ? getMonospaceFont() : null); }
    public static void showMessage(String message, String title, int dialogType, Font textFont) {
        JOptionPane.showMessageDialog(owningComponent, createTextPanel(message, textFont), title, dialogType);
    }
    public static int prompt(String message) { return prompt(message, programTitle, PROMPT_OK_CANCEL); }
    public static int prompt(String message, String title, int promptType) {
        return JOptionPane.showConfirmDialog(owningComponent, createTextPanel(message, null), title, promptType, QUESTION);
    }
    public static int prompt(String message, String title, Object[] options, Object initialValue) {
        return JOptionPane.showOptionDialog(owningComponent, createTextPanel(message, null), title, PROMPT_NONE, QUESTION, null, options, initialValue);
    }

    public static String input(String message, String title, String initialValue) {
        return (String)JOptionPane.showInputDialog(owningComponent, createTextPanel(message, null), title, QUESTION, null, null, initialValue);
    }

    public static JPanel createTextPanel(String content, @Nullable Font useFont) {
        JPanel p = new JPanel(new GridBagLayout());
        //JTextArea text = new JTextArea();
        JTextPane text = new JTextPane();
        text.setFont(useFont == null ? p.getFont() : useFont);
        //text.setLineWrap(true);
        //text.setWrapStyleWord(true);
        text.setEditable(false);
        JScrollPane scroll = new JScrollPane(text);
        scroll.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.weighty = gbc.weightx = 1;
        gbc.fill = GridBagConstraints.VERTICAL;
        p.add(scroll, gbc);
        //p.setPreferredSize(new Dimension(560, 300));

        final int TARGET_WIDTH = 600; // width of ScrollPane. TextArea may be smaller, if scrollbars are need. Dialog will be larger.
        Rectangle screen = getScreenWorkingArea(null);

        // Get the maximum width and height for our text
        final int wmax = Math.min((int)(screen.width * 0.6f), TARGET_WIDTH);
        final int hmax = (int)(screen.height * 0.7f);

        // Set the text width and query its preferred height.  (for now, assume no vertical scroll)
        text.setSize(new Dimension(wmax, 10));
        if (content.startsWith("<html>"))
            text.setContentType("text/html");
        text.setBackground(p.getBackground());
        text.setForeground(p.getForeground());
        text.setText(content);
        scroll.setPreferredSize(new Dimension(wmax, Math.min(text.getPreferredSize().height, hmax)));

        javax.swing.SwingUtilities.invokeLater(() -> {
            scroll.getVerticalScrollBar().setValue(0);
            scroll.getHorizontalScrollBar().setValue(0);
        });
        return p;
    }
    /**
     * getScreenWorkingArea, This returns the working area of the screen. (The working area excludes
     * any task bars.) This function accounts for multi-monitor setups. If a window is supplied,
     * then the the monitor that contains the window will be used. If a window is not supplied, then
     * the primary monitor will be used.
     */
    static public Rectangle getScreenWorkingArea(Window windowOrNull) {
        Insets insets;
        Rectangle bounds;
        if (windowOrNull == null) {
            GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
            insets = Toolkit.getDefaultToolkit().getScreenInsets(ge.getDefaultScreenDevice()
                    .getDefaultConfiguration());
            bounds = ge.getDefaultScreenDevice().getDefaultConfiguration().getBounds();
        } else {
            GraphicsConfiguration gc = windowOrNull.getGraphicsConfiguration();
            insets = windowOrNull.getToolkit().getScreenInsets(gc);
            bounds = gc.getBounds();
        }
        bounds.x += insets.left;
        bounds.y += insets.top;
        bounds.width -= (insets.left + insets.right);
        bounds.height -= (insets.top + insets.bottom);
        return bounds;
    }
    public static Font getDefaultFont() {
        if (defaultFont == null) defaultFont = Font.decode(null);
        return defaultFont;
    }
    public static Font getMonospaceFont() {
        if (monospaceFont == null) {
            Font f = getDefaultFont(); // temporary, to get default size etc.
            monospaceFont = new Font(Font.MONOSPACED, f.getStyle(), f.getSize());
        }
        return monospaceFont;
    }
}
