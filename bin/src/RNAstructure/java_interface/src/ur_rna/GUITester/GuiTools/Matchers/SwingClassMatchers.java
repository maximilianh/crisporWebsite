package ur_rna.GUITester.GuiTools.Matchers;

import javax.swing.*;
import javax.swing.text.JTextComponent;
import java.awt.*;

public abstract class SwingClassMatchers {
    public static final ClassMatcher menuBar = ClassMatcher.of(JMenuBar.class);
    public static final ClassMatcher menu = ClassMatcher.of(JMenu.class);
    //public static final ClassMatcher menu = ClassMatcher.of(MenuElement.class);
    public static final ClassMatcher menuItem = ClassMatcher.of(JMenuItem.class);
    public static final ClassMatcher button = ClassMatcher.of(JButton.class);
    public static final ClassMatcher text = ClassMatcher.of(JTextComponent.class);
    public static final ClassMatcher field = ClassMatcher.of(JTextField.class);
    public static final ClassMatcher textArea = ClassMatcher.of(JTextArea.class);
    public static final ClassMatcher check = ClassMatcher.of(JCheckBox.class);
    public static final ClassMatcher toggle = ClassMatcher.of(JToggleButton.class);
    public static final ClassMatcher radio = ClassMatcher.of(JRadioButton.class);
    public static final ClassMatcher list = ClassMatcher.of(JList.class);
    public static final ClassMatcher spinner = ClassMatcher.of(JSpinner.class);
    public static final ClassMatcher combo = ClassMatcher.of(JComboBox.class);
    public static final ClassMatcher label = ClassMatcher.of(JLabel.class);
    public static final ClassMatcher panel = ClassMatcher.of(JPanel.class);
    public static final ClassMatcher dialog = ClassMatcher.of(Dialog.class);
    public static final ClassMatcher window = ClassMatcher.of(Window.class);
}
