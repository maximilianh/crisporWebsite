package ur_rna.StructureEditor.windows;

import ur_rna.RNAstructure.backend.RNA;
import ur_rna.RNAstructure.backend.RNABackend;
import ur_rna.StructureEditor.FileType;
import ur_rna.StructureEditor.Program;
import ur_rna.StructureEditor.models.RnaScene;
import ur_rna.StructureEditor.models.RnaSceneGroup;
import ur_rna.StructureEditor.services.SceneColorizer;
import ur_rna.StructureEditor.services.SceneColorizer.ColorMode;
import ur_rna.StructureEditor.services.fileIO.RnaFileIO;
import ur_rna.Utilities.Colors;
import ur_rna.Utilities.PathTools;
import ur_rna.Utilities.Strings;
import ur_rna.Utilities.annotation.Nullable;
import ur_rna.Utilities.swing.Dialogs;
import ur_rna.Utilities.swing.FileDialog;
import ur_rna.Utilities.swing.TextDocListener;

import javax.swing.*;
import javax.swing.Timer;
import javax.swing.text.BadLocationException;
import javax.swing.text.JTextComponent;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.function.Consumer;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static ur_rna.Utilities.ObjTools.indexOf;

/**
 * Allows the user to choose how to color nucleobases.
 */
public class ColorizeDialog extends JDialog implements ActionListener {
    private JPanel pnlMain;
    private JButton btnClose;
    private JButton btnApplyText;
    private JRadioButton optStructuresAll;
    private JTabbedPane tabs;
    private JTextArea txtBase;
    private JCheckBox chkUse;
    private JRadioButton optBasesSelected;
    private JButton btnChooseColor;
    private JTextArea txtColor;
    private JButton btnAddSymbolRule;
    private JPanel pnlSymbolColors;
    private JLabel lblColor;
    private JRadioButton optStructuresCurrent;
    private JRadioButton optBasesAll;
    private JButton btnResetSymbolRules;
    private JTextField txtSingleColor;
    private JLabel lblSingleColor;
    private JButton btnSingleColor;
    private JCheckBox chkIgnoreCase;
    private JTextArea txtModData;
    private JButton btnLoadModFile;
    private JButton btnHelpMod;
    private JTextField txtValueLimit;
    private JTextField txtValueColor;
    private JLabel lblValueColor;
    private JButton btnResetValueRules;
    private JButton btnValueChooseColor;
    private JPanel pnlValueColors;
    private JButton btnAddValueRule;
    private JButton btnRemValueRule;
    private JButton btnClearValueText;
    private JComboBox cmbValueComparison;
    private JLabel lblValueColorHeader;
    private JTextArea txtProbData;
    private JLabel lblProbColorHeader;
    private JPanel pnlProbColors;
    private JButton btnAddProbRule;
    private JButton btnRemProbRule;
    private JButton btnResetProbRules;
    private JButton btnLoadProbFile;
    private JButton btnClearProbText;
    private JButton btnHelpProb;
    private JButton btnRemSymbolRule;
    private JComboBox cmbProbComparison;
    private JButton btnApplyFill;
    private JButton btnApplyOutline;
    private JButton btnApplyBond;
    private JButton btnResetText;
    private JButton btnResetFill;
    private JButton btnResetOutline;
    private JButton btnResetBond;
    private JButton btnResetAll;
    private Timer updateTimer = new Timer(100, this::updateTimerElapsed);
    private TextDocListener docListener = new TextDocListener();
    private Preferences prefs;
    private SimpleColorControl singleColor;
    private String lastValueDataPath, lastProbDataPath;
    private Font buttonFont, buttonFontBold;

    public ColorizeDialog(DrawWindow targetWindow, @Nullable Preferences userPreferencesNode, boolean modal) {
        super(SwingUtilities.getWindowAncestor(targetWindow), "Color Annotations", modal ? ModalityType.APPLICATION_MODAL : ModalityType.MODELESS);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        prefs = userPreferencesNode;

        addActionListeners(btnApplyText, btnApplyFill, btnApplyOutline, btnApplyBond,
                btnResetText, btnResetFill, btnResetOutline, btnResetBond, btnResetAll,
                btnClose, btnAddSymbolRule, btnAddValueRule, btnAddProbRule);

        cmbValueComparison.addItemListener(this::valueComparisonChanged);
        cmbProbComparison.addItemListener(this::probComparisonChanged);

        btnClearValueText.addActionListener(e -> txtModData.setText(""));
        btnClearProbText.addActionListener(e -> txtProbData.setText(""));
        btnResetSymbolRules.addActionListener(e -> resetSymbolColorRules());
        btnResetProbRules.addActionListener(e -> resetProbColorRules());
        btnResetValueRules.addActionListener(e -> resetValueColorRules());

        btnRemSymbolRule.addActionListener(e -> removeLastRule(symbolColorRules));
        btnRemValueRule.addActionListener(e -> removeLastRule(valueColorRules));
        btnRemProbRule.addActionListener(e -> removeLastRule(probColorRules));

        btnSingleColor.addActionListener(e -> showColorChooser(singleColor));
        btnLoadModFile.addActionListener(e -> browseLoadValueData());
        btnLoadProbFile.addActionListener(e -> browseLoadProbData());
        btnHelpMod.addActionListener(e -> showHelp("mod"));

        singleColor = new SimpleColorControl(txtSingleColor, lblSingleColor);
        docListener.listen(txtSingleColor);
        singleColor.addBoxClickHandler(this);

        resetProbColorRules();
        resetValueColorRules();
        resetSymbolColorRules();

//        if (targetWindow.getScenes().size() == 1) {
//            optStructuresCurrent.setSelected(true);
//            optStructuresCurrent.setEnabled(false);
//            optStructuresAll.setEnabled(false);
//        } else {
//            optStructuresAll.setSelected(true);
//        }

        optStructuresAll.setSelected(true);
        optBasesAll.setSelected(true);
        //optStructuresAll.addItemListener(e -> updateUI());

        setLayout(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets.set(5, 5, 5, 5);
        gbc.fill = GridBagConstraints.BOTH;
        gbc.weightx = gbc.weighty = 1;
        add(pnlMain, gbc);
        pack();
        setSize(new Dimension(700, getHeight()));

        buttonFont = btnApplyText.getFont();
        buttonFontBold = buttonFont.deriveFont(Font.BOLD);

        loadSettings();
        updateUI();
        updateTimer.start();
        SwingUtilities.invokeLater(() -> {
            txtProbData.select(0, 0);
            txtModData.select(0, 0);
        });
    }

    private void addActionListeners(AbstractButton... buttons) {
        for (AbstractButton b : buttons) b.addActionListener(this);
    }
    private void showHelp(final String section) {
        Dialogs.showMessage(Program.getResourceText("colorize-help", section), "Colorize Drawing Help", Dialogs.INFO);
    }

    private void removeNonHeaders(JComponent panel) {
        for (int i = panel.getComponentCount() - 1; i > -1; i--) {
            String name = panel.getComponent(i).getName();
            if (name == null || !name.startsWith("header"))
                panel.remove(i);
        }
    }
    private void resetSymbolColorRules() {
        removeNonHeaders(pnlSymbolColors);
        symbolColorRules.clear();
//        for (Component c : new Component[]{lblColor, txtBase, txtColor, chkUse, btnChooseColor})
//            pnlSymbolColors.remove(c);
        addSymbolColorRule("A a", Colors.Red);
        addSymbolColorRule("U u T t", Colors.Orange);
        addSymbolColorRule("G g", Colors.Blue);
        addSymbolColorRule("C c", Colors.Green);
        addSymbolColorRule("X x N n", Colors.Gray);
        addSymbolColorRule("PpRrZz", Colors.Purple);
        addSymbolColorRule("XxYy", Colors.Teal);
        pnlSymbolColors.revalidate();
    }

    private void resetValueColorRules() {
        correctPreferredHeight(txtValueColor, txtValueLimit);
        lblValueColor.setPreferredSize(new Dimension(lblValueColor.getPreferredSize().width, txtValueColor.getPreferredSize().height));
        removeNonHeaders(pnlValueColors);
        valueColorRules.clear();
//        for (Component c : new Component[]{lblColor, txtBase, txtColor, chkUse, btnChooseColor})
//            pnlSymbolColors.remove(c);
        addValueColorRule("0.85", Colors.Red);
        addValueColorRule("0.40", Colors.Orange);
        addValueColorRule("-INF", Colors.Black);
        setValueComparison(">=");
        pnlValueColors.revalidate();
    }

    private void resetProbColorRules() {
        correctPreferredHeight(txtValueColor, txtValueLimit);
        lblValueColor.setPreferredSize(new Dimension(lblValueColor.getPreferredSize().width, txtValueColor.getPreferredSize().height));
        removeNonHeaders(pnlProbColors);
        probColorRules.clear();
//        for (Component c : new Component[]{lblColor, txtBase, txtColor, chkUse, btnChooseColor})
//            pnlSymbolColors.remove(c);
        addProbColorRule(".99", Colors.Red);
        addProbColorRule(".95", Colors.Orange);
        addProbColorRule(".90", Colors.Yellow);
        addProbColorRule(".80", Colors.DarkGreen);
        addProbColorRule(".70", Colors.LimeGreen);
        addProbColorRule(".60", Colors.DodgerBlue);
        addProbColorRule(".50", Colors.DarkBlue);
        addProbColorRule("0", Colors.HotPink);
        setProbComparison(">=");
        pnlValueColors.revalidate();
    }

    private void updateUI() {
//        if (optStructuresAll.isSelected()) {
//            optBasesAll.setSelected(true);
//            optBasesAll.setEnabled(false);
//            optBasesSelected.setEnabled(false);
//        } else { //if (target.controller.getSelected().length != 0) {
//            optBasesAll.setEnabled(true);
//            optBasesSelected.setEnabled(true);
//        }

        if (_lastColorEffect == null) _lastColorEffect = "<none>";
        JRootPane rootPane = SwingUtilities.getRootPane(this);
        for (JButton b : new JButton[]{btnApplyText, btnApplyFill, btnApplyOutline, btnApplyBond}) {
            if (_lastColorEffect.equals(b.getActionCommand())) {
                b.setFont(buttonFontBold);
                if (rootPane != null) rootPane.setDefaultButton(b);
            } else
                b.setFont(buttonFont);
        }
    }

    private void updateTimerElapsed(final ActionEvent event) {
        if (!docListener.resetIfChanged()) return;
        docListener.suspendUpdates();
        try {
            for (ColorEditControl r : probColorRules)
                r.setColorFromText();
            for (ColorEditControl r : valueColorRules)
                r.setColorFromText();
            for (ColorEditControl r : symbolColorRules)
                r.setColorFromText();
            singleColor.setColorFromText();
        } finally {
            docListener.resumeUpdates();
        }
    }

    private static String[] valueComparisonSymbols = "< <= > >=".split(" ");
    private void setValueComparison(String symbol) { cmbValueComparison.setSelectedIndex(indexOf(valueComparisonSymbols, symbol)); }
    private String getValueComparisonSymbol() {
        int pos = cmbValueComparison.getSelectedIndex();
        return valueComparisonSymbols[pos == -1 ? 0 : pos];
    }
    private void valueComparisonChanged(final ItemEvent event) {
        if (event.getStateChange() == ItemEvent.SELECTED) {
            lblValueColorHeader.setText(String.format("Color (if value %s Limit)", getValueComparisonSymbol()));
        }
    }
    private void setProbComparison(String symbol) {
        cmbProbComparison.setSelectedIndex(indexOf(valueComparisonSymbols, symbol));
    }
    private String getProbComparisonSymbol() {
        int pos = cmbProbComparison.getSelectedIndex();
        return valueComparisonSymbols[pos == -1 ? 0 : pos];
    }
    private void probComparisonChanged(final ItemEvent event) {
        if (event.getStateChange() == ItemEvent.SELECTED) {
            lblProbColorHeader.setText(String.format("Color (if prob. %s Limit)", getProbComparisonSymbol()));
        }
    }

    private void browseLoadValueData() {
        String path = FileDialog.getOpenName(FileType.SHAPE.getFilterString() + ";txt;chem;mod|*", lastValueDataPath, "Open Data File", "chem-data", this);
        if (path != null) {
            lastValueDataPath = path;
            try {
                txtModData.setText(Strings.readAll(new FileInputStream(path)));
                SwingUtilities.invokeLater(() -> txtModData.select(0, 0));
            } catch (IOException ex) {
                Dialogs.showWarning("The data file could not be read: " + ex.getMessage(), "Error Reading File");
            }
        }
    }
    private void browseLoadProbData() {
        String path = FileDialog.getOpenName(FileType.PartitionSav.getFilterString() + "|" + "Probability Text File|txt;prob;part|*", lastProbDataPath, "Open Probability Data File", "prob-data", this);
        if (path != null) {
            lastProbDataPath = path;
            if (PathTools.getExt(path, false, false).toUpperCase().equals("PFS")) {
                RNA rna = new RNA(path, RNABackend.FILE_PFS, RnaFileIO.ALPHABET_RNA, true, true);
                if (rna.GetErrorCode() == 0) {
                    int size = rna.GetSequenceLength();
                    double[] probs = new double[(size - 1) * size / 2];
                    int result = rna.GetPairProbabilities(probs, probs.length);
                    if (result < 0) {
                        Dialogs.showWarning("Error reading probabilities: " + RNA.GetErrorMessage(-result) + "\nFile: " + path, "Data File Error");
                        return;
                    }
                    StringBuilder sb = new StringBuilder();
                    File pfs = new File(path);
                    sb.append("# PFS File: ").append(pfs.getName()).append('\n')
                            .append("# Modified: ").append(new Date(pfs.lastModified()).toString()).append('\n')
                            .append("# Sequence-Length: ").append(size).append('\n');
                    int pos = 0;
                    String format = "%" + (Integer.toString(size).length() + 1) + "d"; // e.g. "%4d" if size < 999
                    for (int i = 0; i < size; i++)
                        for (int j = i + 1; j < size; j++)
                            sb.append(String.format(format, i + 1)).append(String.format(format, j + 1)).append(' ').append(probs[pos++]).append('\n');
                    txtProbData.setText(sb.toString());
                    SwingUtilities.invokeLater(() -> txtProbData.select(0, 0));
                } else {
                    Dialogs.showWarning("The PFS file could not be read: " + rna.GetFullErrorMessage(), "Error Reading PFS File");
                }
            } else {
                try {
                    txtProbData.setText(Strings.readAll(new FileInputStream(path)));
                    SwingUtilities.invokeLater(() -> txtProbData.select(0, 0));
                } catch (IOException ex) {
                    Dialogs.showWarning("The data file could not be read: " + ex.getMessage(), "Error Reading File");
                }
            }
        }
    }

    public void showDialog() {
        setVisible(true);
    }
    {
// GUI initializer generated by IntelliJ IDEA GUI Designer
// DO NOT EDIT OR ADD ANY CODE HERE!
        $$$setupUI$$$();
    }

    /**
     * Method generated by IntelliJ IDEA GUI Designer
     * DO NOT edit this method OR call it in your code!
     *
     * @noinspection ALL
     */
    private void $$$setupUI$$$() {
        pnlMain = new JPanel();
        pnlMain.setLayout(new GridBagLayout());
        tabs = new JTabbedPane();
        GridBagConstraints gbc;
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        gbc.insets = new Insets(0, 0, 5, 0);
        pnlMain.add(tabs, gbc);
        tabs.setBorder(BorderFactory.createTitledBorder(BorderFactory.createRaisedBevelBorder(), null));
        final JPanel panel1 = new JPanel();
        panel1.setLayout(new GridBagLayout());
        panel1.setName("prob");
        tabs.addTab("Basepair Probability", panel1);
        final JLabel label1 = new JLabel();
        label1.setAutoscrolls(false);
        label1.setText("<html>This tool sets the color of each base pair according the probability of it being basepaired (or the probability that unpaired nucleotides are unpaired).  \nData can be loaded from a <b>Partition Save File (<i>pfs</i>)</b> obtained from running <b><tt>partition</tt></b> or from a text file that lists pairwise basepairing probabilities in three columns: <tt>BASE1  BASE2  PROB(BASE1:BASE2) </tt><br>\nClick the <b>Help</b> button for more information.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(3, 3, 8, 3);
        panel1.add(label1, gbc);
        final JSplitPane splitPane1 = new JSplitPane();
        splitPane1.setContinuousLayout(true);
        splitPane1.setDividerLocation(230);
        splitPane1.setOneTouchExpandable(false);
        splitPane1.setResizeWeight(0.5);
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        panel1.add(splitPane1, gbc);
        final JPanel panel2 = new JPanel();
        panel2.setLayout(new GridBagLayout());
        splitPane1.setLeftComponent(panel2);
        final JLabel label2 = new JLabel();
        label2.setFont(new Font(label2.getFont().getName(), Font.BOLD, label2.getFont().getSize()));
        label2.setName("headerData");
        label2.setText("Probability Data");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 4, 0, 0);
        panel2.add(label2, gbc);
        final JScrollPane scrollPane1 = new JScrollPane();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        panel2.add(scrollPane1, gbc);
        txtProbData = new JTextArea();
        txtProbData.setText("1  1\t0.949369084\n1  2\t0.008434194\n1  3\t0.115617068\n1  4\t0.182823643\n2  2\t0.656706731\n2  3\t0.761605327\n2  4\t0.115617068\n...");
        scrollPane1.setViewportView(txtProbData);
        final JPanel panel3 = new JPanel();
        panel3.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.VERTICAL;
        panel2.add(panel3, gbc);
        btnLoadProbFile = new JButton();
        btnLoadProbFile.setText("Load File...");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel3.add(btnLoadProbFile, gbc);
        btnHelpProb = new JButton();
        btnHelpProb.setText("Help");
        btnHelpProb.setMnemonic('H');
        btnHelpProb.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel3.add(btnHelpProb, gbc);
        btnClearProbText = new JButton();
        btnClearProbText.setActionCommand("clearModText");
        btnClearProbText.setText("Clear");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel3.add(btnClearProbText, gbc);
        final JPanel panel4 = new JPanel();
        panel4.setLayout(new GridBagLayout());
        splitPane1.setRightComponent(panel4);
        final JScrollPane scrollPane2 = new JScrollPane();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        panel4.add(scrollPane2, gbc);
        final JPanel panel5 = new JPanel();
        panel5.setLayout(new GridBagLayout());
        scrollPane2.setViewportView(panel5);
        pnlProbColors = new JPanel();
        pnlProbColors.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 3;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        panel5.add(pnlProbColors, gbc);
        final JLabel label3 = new JLabel();
        label3.setFont(new Font(label3.getFont().getName(), Font.BOLD, label3.getFont().getSize()));
        label3.setName("headerLimit");
        label3.setText("Probability Limit");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(5, 5, 5, 5);
        pnlProbColors.add(label3, gbc);
        lblProbColorHeader = new JLabel();
        lblProbColorHeader.setFont(new Font(lblProbColorHeader.getFont().getName(), Font.BOLD, lblProbColorHeader.getFont().getSize()));
        lblProbColorHeader.setName("headerColor");
        lblProbColorHeader.setRequestFocusEnabled(true);
        lblProbColorHeader.setText("Color (if prob. < Limit)");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.gridwidth = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(5, 5, 5, 5);
        pnlProbColors.add(lblProbColorHeader, gbc);
        final JButton button1 = new JButton();
        button1.setActionCommand("color");
        button1.setText("...");
        button1.setToolTipText("Choose Color...");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlProbColors.add(button1, gbc);
        final JLabel label4 = new JLabel();
        label4.setAutoscrolls(false);
        label4.setBackground(new Color(-16711936));
        label4.setOpaque(true);
        label4.setPreferredSize(new Dimension(40, 24));
        label4.setText("");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlProbColors.add(label4, gbc);
        final JTextField textField1 = new JTextField();
        textField1.setPreferredSize(new Dimension(100, 24));
        textField1.setText("INF");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(2, 2, 2, 2);
        pnlProbColors.add(textField1, gbc);
        final JTextField textField2 = new JTextField();
        textField2.setPreferredSize(new Dimension(100, 24));
        textField2.setText("Lime");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlProbColors.add(textField2, gbc);
        final JLabel label5 = new JLabel();
        label5.setMaximumSize(new Dimension(200, 300));
        label5.setText("INF = infinity;  -INF = negative infinity");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.gridwidth = 3;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(4, 2, 4, 2);
        panel5.add(label5, gbc);
        final JLabel label6 = new JLabel();
        label6.setText("Probability must be ");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        panel5.add(label6, gbc);
        cmbProbComparison = new JComboBox();
        cmbProbComparison.setActionCommand("valueComparisonChanged");
        final DefaultComboBoxModel defaultComboBoxModel1 = new DefaultComboBoxModel();
        defaultComboBoxModel1.addElement("<  less than");
        defaultComboBoxModel1.addElement("<= less than or equal to");
        defaultComboBoxModel1.addElement(">  greater than");
        defaultComboBoxModel1.addElement(">= greater than or equal to");
        cmbProbComparison.setModel(defaultComboBoxModel1);
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        panel5.add(cmbProbComparison, gbc);
        final JLabel label7 = new JLabel();
        label7.setText(" each limit.");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        panel5.add(label7, gbc);
        final JLabel label8 = new JLabel();
        label8.setFont(new Font(label8.getFont().getName(), Font.BOLD, label8.getFont().getSize()));
        label8.setName("headerValues");
        label8.setText("Probability Ranges and Colors:");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 4, 0, 0);
        panel4.add(label8, gbc);
        final JPanel panel6 = new JPanel();
        panel6.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.VERTICAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel4.add(panel6, gbc);
        btnAddProbRule = new JButton();
        btnAddProbRule.setActionCommand("addProbRule");
        btnAddProbRule.setText("Add Rule");
        btnAddProbRule.setMnemonic('A');
        btnAddProbRule.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel6.add(btnAddProbRule, gbc);
        btnResetProbRules = new JButton();
        btnResetProbRules.setText("Reset List");
        btnResetProbRules.setMnemonic('R');
        btnResetProbRules.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel6.add(btnResetProbRules, gbc);
        btnRemProbRule = new JButton();
        btnRemProbRule.setActionCommand("remProbRule");
        btnRemProbRule.setText("Remove Last Rule");
        btnRemProbRule.setMnemonic('A');
        btnRemProbRule.setDisplayedMnemonicIndex(8);
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel6.add(btnRemProbRule, gbc);
        final JPanel panel7 = new JPanel();
        panel7.setLayout(new GridBagLayout());
        panel7.setName("mod");
        tabs.addTab("Chemical Modification (e.g. SHAPE)", panel7);
        final JLabel label9 = new JLabel();
        label9.setAutoscrolls(false);
        label9.setText("<html>This tool sets the color of each base according to a data table arranged in two columns.<br>\nThe data can represent chemical modification reactvities (e.g. SHAPE) or any arbitrary criteria.<br>\nIt is <b>not</b> necessary to include a line for each base (unlisted bases will not be colored).<br>\nClick the <b>Help</b> button for more information.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(3, 3, 8, 3);
        panel7.add(label9, gbc);
        final JSplitPane splitPane2 = new JSplitPane();
        splitPane2.setContinuousLayout(true);
        splitPane2.setOneTouchExpandable(false);
        splitPane2.setResizeWeight(0.5);
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        panel7.add(splitPane2, gbc);
        final JPanel panel8 = new JPanel();
        panel8.setLayout(new GridBagLayout());
        splitPane2.setLeftComponent(panel8);
        final JLabel label10 = new JLabel();
        label10.setFont(new Font(label10.getFont().getName(), Font.BOLD, label10.getFont().getSize()));
        label10.setName("headerData");
        label10.setText("Numeric Data");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 4, 0, 0);
        panel8.add(label10, gbc);
        final JScrollPane scrollPane3 = new JScrollPane();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        panel8.add(scrollPane3, gbc);
        txtModData = new JTextArea();
        txtModData.setText("1\t0.949369084\n2\t0.008434194\n3\t0.115617068\n4\t0.182823643\n5\t0.656706731\n6\t0.761605327");
        scrollPane3.setViewportView(txtModData);
        final JPanel panel9 = new JPanel();
        panel9.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.VERTICAL;
        panel8.add(panel9, gbc);
        btnLoadModFile = new JButton();
        btnLoadModFile.setText("Load File...");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel9.add(btnLoadModFile, gbc);
        btnHelpMod = new JButton();
        btnHelpMod.setText("Help");
        btnHelpMod.setMnemonic('H');
        btnHelpMod.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel9.add(btnHelpMod, gbc);
        btnClearValueText = new JButton();
        btnClearValueText.setActionCommand("clearModText");
        btnClearValueText.setText("Clear");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel9.add(btnClearValueText, gbc);
        final JPanel panel10 = new JPanel();
        panel10.setLayout(new GridBagLayout());
        splitPane2.setRightComponent(panel10);
        final JScrollPane scrollPane4 = new JScrollPane();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        panel10.add(scrollPane4, gbc);
        final JPanel panel11 = new JPanel();
        panel11.setLayout(new GridBagLayout());
        scrollPane4.setViewportView(panel11);
        pnlValueColors = new JPanel();
        pnlValueColors.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 3;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        panel11.add(pnlValueColors, gbc);
        final JLabel label11 = new JLabel();
        label11.setFont(new Font(label11.getFont().getName(), Font.BOLD, label11.getFont().getSize()));
        label11.setName("headerLimit");
        label11.setText("Numeric Limit");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(5, 5, 5, 5);
        pnlValueColors.add(label11, gbc);
        lblValueColorHeader = new JLabel();
        lblValueColorHeader.setFont(new Font(lblValueColorHeader.getFont().getName(), Font.BOLD, lblValueColorHeader.getFont().getSize()));
        lblValueColorHeader.setName("headerColor");
        lblValueColorHeader.setRequestFocusEnabled(true);
        lblValueColorHeader.setText("Color (if value < Limit)");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.gridwidth = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(5, 5, 5, 5);
        pnlValueColors.add(lblValueColorHeader, gbc);
        btnValueChooseColor = new JButton();
        btnValueChooseColor.setActionCommand("color");
        btnValueChooseColor.setText("...");
        btnValueChooseColor.setToolTipText("Choose Color...");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlValueColors.add(btnValueChooseColor, gbc);
        lblValueColor = new JLabel();
        lblValueColor.setAutoscrolls(false);
        lblValueColor.setBackground(new Color(-16711936));
        lblValueColor.setOpaque(true);
        lblValueColor.setPreferredSize(new Dimension(40, 24));
        lblValueColor.setText("");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlValueColors.add(lblValueColor, gbc);
        txtValueLimit = new JTextField();
        txtValueLimit.setPreferredSize(new Dimension(100, 24));
        txtValueLimit.setText("INF");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(2, 2, 2, 2);
        pnlValueColors.add(txtValueLimit, gbc);
        txtValueColor = new JTextField();
        txtValueColor.setPreferredSize(new Dimension(100, 24));
        txtValueColor.setText("Lime");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlValueColors.add(txtValueColor, gbc);
        final JLabel label12 = new JLabel();
        label12.setMaximumSize(new Dimension(200, 300));
        label12.setText("INF = infinity;  -INF = negative infinity");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.gridwidth = 3;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(4, 2, 4, 2);
        panel11.add(label12, gbc);
        final JLabel label13 = new JLabel();
        label13.setText("Value must be ");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        panel11.add(label13, gbc);
        cmbValueComparison = new JComboBox();
        cmbValueComparison.setActionCommand("valueComparisonChanged");
        final DefaultComboBoxModel defaultComboBoxModel2 = new DefaultComboBoxModel();
        defaultComboBoxModel2.addElement("<  less than");
        defaultComboBoxModel2.addElement("<= less than or equal to");
        defaultComboBoxModel2.addElement(">  greater than");
        defaultComboBoxModel2.addElement(">= greater than or equal to");
        cmbValueComparison.setModel(defaultComboBoxModel2);
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        panel11.add(cmbValueComparison, gbc);
        final JLabel label14 = new JLabel();
        label14.setText(" each limit.");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        panel11.add(label14, gbc);
        final JLabel label15 = new JLabel();
        label15.setFont(new Font(label15.getFont().getName(), Font.BOLD, label15.getFont().getSize()));
        label15.setName("headerValues");
        label15.setText("Value Ranges and Colors:");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 4, 0, 0);
        panel10.add(label15, gbc);
        final JPanel panel12 = new JPanel();
        panel12.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.VERTICAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel10.add(panel12, gbc);
        btnAddValueRule = new JButton();
        btnAddValueRule.setActionCommand("addValueRule");
        btnAddValueRule.setText("Add Rule");
        btnAddValueRule.setMnemonic('A');
        btnAddValueRule.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel12.add(btnAddValueRule, gbc);
        btnResetValueRules = new JButton();
        btnResetValueRules.setText("Reset List");
        btnResetValueRules.setMnemonic('R');
        btnResetValueRules.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel12.add(btnResetValueRules, gbc);
        btnRemValueRule = new JButton();
        btnRemValueRule.setActionCommand("remValueRule");
        btnRemValueRule.setText("Remove Last Rule");
        btnRemValueRule.setMnemonic('A');
        btnRemValueRule.setDisplayedMnemonicIndex(8);
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel12.add(btnRemValueRule, gbc);
        final JPanel panel13 = new JPanel();
        panel13.setLayout(new GridBagLayout());
        panel13.setName("base");
        tabs.addTab("Nucleobase Identity", panel13);
        final JLabel label16 = new JLabel();
        label16.setAutoscrolls(false);
        label16.setText("<html>This tool sets the color of each base according to its identity (i.e. the symbol/character used to represent it).<br>\nYou can enter as many or as few symbols as desired. Nucleobases that are not listed will not be affected.<br>\nYou can disable a rule below by unchecking the box next to it or by deleting the symbol from the text box.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(3, 3, 3, 3);
        panel13.add(label16, gbc);
        final JScrollPane scrollPane5 = new JScrollPane();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        panel13.add(scrollPane5, gbc);
        final JPanel panel14 = new JPanel();
        panel14.setLayout(new GridBagLayout());
        panel14.setAlignmentX(0.0f);
        panel14.setAlignmentY(0.0f);
        scrollPane5.setViewportView(panel14);
        pnlSymbolColors = new JPanel();
        pnlSymbolColors.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 4;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        panel14.add(pnlSymbolColors, gbc);
        final JLabel label17 = new JLabel();
        label17.setFont(new Font(label17.getFont().getName(), Font.BOLD, 14));
        label17.setName("headerSymbol");
        label17.setText("Nucleobase(s)");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(5, 5, 5, 5);
        pnlSymbolColors.add(label17, gbc);
        txtBase = new JTextArea();
        txtBase.setText("A");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlSymbolColors.add(txtBase, gbc);
        chkUse = new JCheckBox();
        chkUse.setText("");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        pnlSymbolColors.add(chkUse, gbc);
        txtColor = new JTextArea();
        txtColor.setPreferredSize(new Dimension(150, 17));
        txtColor.setText("Red");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlSymbolColors.add(txtColor, gbc);
        final JLabel label18 = new JLabel();
        label18.setFont(new Font(label18.getFont().getName(), Font.BOLD, 14));
        label18.setName("headerColor");
        label18.setText("Color (name or RGB code)");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.gridwidth = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(5, 5, 5, 5);
        pnlSymbolColors.add(label18, gbc);
        btnChooseColor = new JButton();
        btnChooseColor.setActionCommand("color");
        btnChooseColor.setText("Choose...");
        gbc = new GridBagConstraints();
        gbc.gridx = 4;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlSymbolColors.add(btnChooseColor, gbc);
        lblColor = new JLabel();
        lblColor.setAutoscrolls(false);
        lblColor.setBackground(new Color(-16711936));
        lblColor.setOpaque(true);
        lblColor.setText("");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.BOTH;
        gbc.insets = new Insets(4, 4, 4, 4);
        pnlSymbolColors.add(lblColor, gbc);
        btnAddSymbolRule = new JButton();
        btnAddSymbolRule.setActionCommand("addSymbolRule");
        btnAddSymbolRule.setText("Add Rule");
        btnAddSymbolRule.setMnemonic('A');
        btnAddSymbolRule.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weighty = 1.0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(0, 20, 0, 0);
        panel14.add(btnAddSymbolRule, gbc);
        btnResetSymbolRules = new JButton();
        btnResetSymbolRules.setActionCommand("resetSymbolRules");
        btnResetSymbolRules.setText("Reset List");
        btnResetSymbolRules.setMnemonic('R');
        btnResetSymbolRules.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.weighty = 1.0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(0, 4, 0, 0);
        panel14.add(btnResetSymbolRules, gbc);
        chkIgnoreCase = new JCheckBox();
        chkIgnoreCase.setActionCommand("");
        chkIgnoreCase.setSelected(true);
        chkIgnoreCase.setText("Ignore Case (G = g etc.)");
        chkIgnoreCase.setToolTipText("If selected the character case of bases will be ignored. So if \"G\" is red then \"g\" will be also.  Otherwise \"G\" and \"g\" are treated as different entities and they must each be listed separately in order to color both.");
        chkIgnoreCase.setVisible(false);
        chkIgnoreCase.putClientProperty("html.disable", Boolean.FALSE);
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(4, 0, 0, 0);
        panel14.add(chkIgnoreCase, gbc);
        btnRemSymbolRule = new JButton();
        btnRemSymbolRule.setActionCommand("remSymbolRule");
        btnRemSymbolRule.setText("Remove Last Rule");
        btnRemSymbolRule.setMnemonic('A');
        btnRemSymbolRule.setDisplayedMnemonicIndex(8);
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(0, 1, 0, 1);
        panel14.add(btnRemSymbolRule, gbc);
        final JPanel panel15 = new JPanel();
        panel15.setLayout(new GridBagLayout());
        panel15.setName("single");
        tabs.addTab("Single Color", panel15);
        final JLabel label19 = new JLabel();
        label19.setText("<html>This tool can be used to color all bases (or all <i>selected</i> bases) the same color.<br>\nFor example, you can select a domain and apply a single color to all bases in that domain.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(4, 4, 4, 4);
        panel15.add(label19, gbc);
        final JPanel panel16 = new JPanel();
        panel16.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.insets = new Insets(4, 4, 0, 0);
        panel15.add(panel16, gbc);
        btnSingleColor = new JButton();
        btnSingleColor.setActionCommand("chooseSingleColor");
        btnSingleColor.setText("Choose...");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        panel16.add(btnSingleColor, gbc);
        lblSingleColor = new JLabel();
        lblSingleColor.setAutoscrolls(false);
        lblSingleColor.setBackground(new Color(-16711936));
        lblSingleColor.setOpaque(true);
        lblSingleColor.setPreferredSize(new Dimension(100, 20));
        lblSingleColor.setText("");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(4, 4, 4, 4);
        panel16.add(lblSingleColor, gbc);
        txtSingleColor = new JTextField();
        txtSingleColor.setPreferredSize(new Dimension(120, 24));
        txtSingleColor.setText("Red");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        panel16.add(txtSingleColor, gbc);
        final JLabel label20 = new JLabel();
        label20.setText("Color (name or RGB code): ");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        panel16.add(label20, gbc);
        final JPanel panel17 = new JPanel();
        panel17.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlMain.add(panel17, gbc);
        final JPanel panel18 = new JPanel();
        panel18.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.ipady = 3;
        gbc.insets = new Insets(0, 0, 0, 10);
        panel17.add(panel18, gbc);
        panel18.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Scope of Operation"));
        final JLabel label21 = new JLabel();
        label21.setText("Structures Affected:");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        panel18.add(label21, gbc);
        optBasesAll = new JRadioButton();
        optBasesAll.setText("All Bases");
        optBasesAll.setMnemonic('B');
        optBasesAll.setDisplayedMnemonicIndex(4);
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 4, 0, 0);
        panel18.add(optBasesAll, gbc);
        optStructuresAll = new JRadioButton();
        optStructuresAll.setText("All Structures");
        optStructuresAll.setMnemonic('S');
        optStructuresAll.setDisplayedMnemonicIndex(4);
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 4, 0, 0);
        panel18.add(optStructuresAll, gbc);
        optStructuresCurrent = new JRadioButton();
        optStructuresCurrent.setText("Current Structure");
        optStructuresCurrent.setMnemonic('C');
        optStructuresCurrent.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 4, 0, 0);
        panel18.add(optStructuresCurrent, gbc);
        optBasesSelected = new JRadioButton();
        optBasesSelected.setText("Selected Bases Only");
        optBasesSelected.setMnemonic('S');
        optBasesSelected.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 4, 0, 0);
        panel18.add(optBasesSelected, gbc);
        final JLabel label22 = new JLabel();
        label22.setText("Nucleobases Affected:");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        panel18.add(label22, gbc);
        final JPanel panel19 = new JPanel();
        panel19.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.gridwidth = 2;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        panel17.add(panel19, gbc);
        btnClose = new JButton();
        btnClose.setActionCommand("close");
        btnClose.setHorizontalTextPosition(0);
        btnClose.setText("Close");
        btnClose.setMnemonic('C');
        btnClose.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 5;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(2, 2, 2, 2);
        panel19.add(btnClose, gbc);
        final JLabel label23 = new JLabel();
        label23.setText("Apply To:");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        panel19.add(label23, gbc);
        btnApplyText = new JButton();
        btnApplyText.setActionCommand("apply-text");
        btnApplyText.setHorizontalAlignment(0);
        btnApplyText.setPreferredSize(new Dimension(180, 32));
        btnApplyText.setText("Text (symbols)");
        btnApplyText.setMnemonic('T');
        btnApplyText.setDisplayedMnemonicIndex(0);
        btnApplyText.setToolTipText("Apply the color operation to the Text (symbols) of bases.");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(2, 2, 2, 2);
        panel19.add(btnApplyText, gbc);
        btnApplyFill = new JButton();
        btnApplyFill.setActionCommand("apply-fill");
        btnApplyFill.setHorizontalAlignment(0);
        btnApplyFill.setPreferredSize(new Dimension(180, 32));
        btnApplyFill.setText("Fill (background)");
        btnApplyFill.setMnemonic('F');
        btnApplyFill.setDisplayedMnemonicIndex(0);
        btnApplyFill.setToolTipText("Apply the color operation to fill in the background of bases.");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(2, 2, 2, 2);
        panel19.add(btnApplyFill, gbc);
        btnApplyOutline = new JButton();
        btnApplyOutline.setActionCommand("apply-line");
        btnApplyOutline.setHorizontalAlignment(0);
        btnApplyOutline.setPreferredSize(new Dimension(180, 32));
        btnApplyOutline.setText("Outlines");
        btnApplyOutline.setMnemonic('O');
        btnApplyOutline.setDisplayedMnemonicIndex(0);
        btnApplyOutline.setToolTipText("Apply the color operation to the circle outline around bases.");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(2, 2, 2, 2);
        panel19.add(btnApplyOutline, gbc);
        btnApplyBond = new JButton();
        btnApplyBond.setActionCommand("apply-bond");
        btnApplyBond.setHorizontalAlignment(0);
        btnApplyBond.setPreferredSize(new Dimension(180, 32));
        btnApplyBond.setText("Base-Pair Bonds");
        btnApplyBond.setMnemonic('B');
        btnApplyBond.setDisplayedMnemonicIndex(10);
        btnApplyBond.setToolTipText("Apply the color operation to the bonds that connect base-pairs. (This has no effect on unpaired bases.)");
        gbc = new GridBagConstraints();
        gbc.gridx = 4;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(2, 2, 2, 2);
        panel19.add(btnApplyBond, gbc);
        final JPanel panel20 = new JPanel();
        panel20.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        gbc.ipady = 3;
        panel17.add(panel20, gbc);
        panel20.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Reset to Defaults"));
        btnResetText = new JButton();
        btnResetText.setActionCommand("reset-text");
        btnResetText.setHorizontalAlignment(0);
        btnResetText.setPreferredSize(new Dimension(180, 32));
        btnResetText.setText("Reset Text");
        btnResetText.setToolTipText("Reset the Text color of bases to the default color.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(2, 2, 2, 2);
        panel20.add(btnResetText, gbc);
        btnResetFill = new JButton();
        btnResetFill.setActionCommand("reset-fill");
        btnResetFill.setHorizontalAlignment(0);
        btnResetFill.setPreferredSize(new Dimension(180, 32));
        btnResetFill.setText("Reset Fill");
        btnResetFill.setToolTipText("Reset the Fill color of bases to the default color.");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(2, 2, 2, 2);
        panel20.add(btnResetFill, gbc);
        btnResetAll = new JButton();
        btnResetAll.setActionCommand("reset-all");
        btnResetAll.setFont(new Font(btnResetAll.getFont().getName(), btnResetAll.getFont().getStyle(), btnResetAll.getFont().getSize()));
        btnResetAll.setHorizontalAlignment(0);
        btnResetAll.setPreferredSize(new Dimension(180, 32));
        btnResetAll.setText("Reset All");
        btnResetAll.setToolTipText("Reset/Clear all custom Color annotations.");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.gridheight = 2;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(2, 2, 2, 2);
        panel20.add(btnResetAll, gbc);
        btnResetOutline = new JButton();
        btnResetOutline.setActionCommand("reset-line");
        btnResetOutline.setHorizontalAlignment(0);
        btnResetOutline.setPreferredSize(new Dimension(180, 32));
        btnResetOutline.setText("Reset Outlines");
        btnResetOutline.setToolTipText("Reset the Outline color of bases to the default color.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(2, 2, 2, 2);
        panel20.add(btnResetOutline, gbc);
        btnResetBond = new JButton();
        btnResetBond.setActionCommand("reset-bond");
        btnResetBond.setHorizontalAlignment(0);
        btnResetBond.setPreferredSize(new Dimension(180, 32));
        btnResetBond.setText("Reset Bonds");
        btnResetBond.setToolTipText("Reset the Bond color of bases to the default color.");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(2, 2, 2, 2);
        panel20.add(btnResetBond, gbc);
        ButtonGroup buttonGroup;
        buttonGroup = new ButtonGroup();
        buttonGroup.add(optStructuresAll);
        buttonGroup.add(optStructuresCurrent);
        buttonGroup = new ButtonGroup();
        buttonGroup.add(optBasesAll);
        buttonGroup.add(optBasesSelected);
    }
    /** @noinspection ALL */
    public JComponent $$$getRootComponent$$$() { return pnlMain; }

    private interface ColorEditControl {
        String getColorText();
        void setColorText(String text);
        void setColor(Color c);
        Color getColor();
        default void setColorFromText() {
            setColor(Colors.getColor(getColorText(), true));
        }
    }

    public void showColorChooser(ColorEditControl c) {
        Color color = JColorChooser.showDialog(this, "Choose a Color", c.getColor());
        if (color != null)
            c.setColor(color);
    }

    private static class SimpleColorControl implements ColorEditControl {
        protected final JTextComponent txt;
        protected final JComponent box;
        public SimpleColorControl(JTextComponent txt, JComponent box) {
            this.txt = txt;
            this.box = box;
        }
        public void setColor(Color c) {
            // avoid setting text if possible. Allow user to edit without interruption.
            if (c != null && !c.equals(Colors.getColor(getColorText(), true)))
                txt.setText(Colors.getName(c));
            box.setBackground(c == null ? Color.GRAY : c);
        }
        public Color getColor() {
            return box.getBackground();
        }
        public String getColorText() { return txt.getText(); }
        public void setColorText(String text) { txt.setText(text); }
        public void addBoxClickHandler(final ColorizeDialog dialog) {
            box.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            addClickHandler(box, e -> dialog.showColorChooser(SimpleColorControl.this));
        }
        public JComponent getPanel() {
            return (JComponent) txt.getParent();
        }
        public void remove() { }
    }

    private static void addClickHandler(JComponent c, final Consumer<MouseEvent> handler) {
        c.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(final MouseEvent e) {
                handler.accept(e);
            }
        });
    }

    private void correctPreferredHeight(JComponent... components) {
        for (JComponent c : components) {
            Dimension d = c.getPreferredSize();
            c.setPreferredSize(null);
            c.setPreferredSize(new Dimension(d.width, c.getPreferredSize().height));
        }
    }
    private static class ValueColorRule extends SimpleColorControl {
        private JTextField txtValue, txtColor;
        private JButton btnChoose;
        private JLabel lblColor;
        public void setValue(String text) { txtValue.setText(text); }
        public String getValue() { return txtValue.getText().trim(); }

        public ValueColorRule(ColorizeDialog dialog) {
            super(new JTextField(), new JLabel(""));
            txtColor = (JTextField) super.txt;
            lblColor = (JLabel) super.box;
            txtValue = new JTextField();
            btnChoose = new JButton(dialog.btnValueChooseColor.getText());
            btnChoose.setToolTipText(dialog.btnValueChooseColor.getToolTipText());

            txtValue.setPreferredSize(dialog.txtValueLimit.getPreferredSize());
            lblColor.setPreferredSize(dialog.lblValueColor.getPreferredSize());
            txtColor.setPreferredSize(dialog.txtValueColor.getPreferredSize());
            lblColor.setOpaque(true);

            addBoxClickHandler(dialog);
            btnChoose.addActionListener(e -> dialog.showColorChooser(this));
            dialog.docListener.listen(txtColor);
        }

        public void addToPanel(JComponent panel) { addToPanel(panel, GridBagConstraints.RELATIVE); }
        public void addToPanel(JComponent panel, int lineIndex) {
            GridBagConstraints gbc = new GridBagConstraints();
            gbc.insets = new Insets(0, 2, 0, 2);
            gbc.gridy = lineIndex;

            gbc.fill = GridBagConstraints.NONE;
            gbc.anchor = GridBagConstraints.WEST;

            gbc.gridx = 0;
            panel.add(txtValue, gbc);

            gbc.gridx = 1;
            panel.add(txtColor, gbc);

            gbc.gridx = 2;
            panel.add(lblColor, gbc);

            gbc.gridx = 3;
            gbc.fill = GridBagConstraints.HORIZONTAL;
            gbc.anchor = GridBagConstraints.CENTER;
            panel.add(btnChoose, gbc);
        }
        @Override
        public void remove() {
            JComponent panel = (JComponent) btnChoose.getParent();
            if (panel == null) return;
            panel.remove(lblColor);
            panel.remove(btnChoose);
            panel.remove(txtColor);
            panel.remove(txtValue);
        }
    }

    private static class SymbolColorRule extends SimpleColorControl {
        private JCheckBox chkEnable;
        private JTextField txtSymbol, txtColor;
        private JButton btnChoose;
        private JLabel lblColor;

        public boolean isEnabled() { return chkEnable.isSelected(); }
        public void setEnabled(boolean enabled) { chkEnable.setSelected(enabled); }
        public void setSymbol(String text) { txtSymbol.setText(text); }
        public String getSymbol() { return txtSymbol.getText().trim(); }

        public SymbolColorRule(ColorizeDialog dialog) {
            super(new JTextField(), new JLabel(""));
            txtColor = (JTextField) super.txt;
            lblColor = (JLabel) super.box;
            chkEnable = new JCheckBox("");
            txtSymbol = new JTextField();
            btnChoose = new JButton("Choose...");

            chkEnable.setSelected(true);
            lblColor.setOpaque(true);
            chkEnable.addItemListener(this::updateEnabledUI);

            txtColor.setPreferredSize(new Dimension(120, txtColor.getPreferredSize().height));

            addBoxClickHandler(dialog);
            //AcceleratorKey.resetTabTraversalKeys(txtColor, txtSymbol);

            btnChoose.addActionListener(e -> dialog.showColorChooser(this));
            dialog.docListener.listen(txtColor);
        }
        private void updateEnabledUI(final ItemEvent event) {
            boolean enabled = isEnabled();
            txtSymbol.setEnabled(enabled);
            txtColor.setEnabled(enabled);
            btnChoose.setEnabled(enabled);
        }
        public void addToPanel(JComponent panel) { addToPanel(panel, GridBagConstraints.RELATIVE); }
        public void addToPanel(JComponent panel, int lineIndex) {
            GridBagConstraints gbc = new GridBagConstraints();
            gbc.insets = new Insets(0, 2, 0, 2);
            gbc.gridy = lineIndex;

            gbc.gridx = 0;
            gbc.fill = GridBagConstraints.NONE;
            gbc.anchor = GridBagConstraints.WEST;
            panel.add(chkEnable, gbc);

            gbc.fill = GridBagConstraints.HORIZONTAL;
            gbc.anchor = GridBagConstraints.CENTER;

            gbc.gridx = 1;
            panel.add(txtSymbol, gbc);

            gbc.gridx = 2;
            panel.add(txtColor, gbc);

            gbc.gridx = 4;
            panel.add(btnChoose, gbc);

            gbc.gridx = 3;
            gbc.fill = GridBagConstraints.BOTH;
            gbc.insets = new Insets(2, 2, 2, 2);
            panel.add(lblColor, gbc);
        }
        @Override
        public void remove() {
            JComponent panel = (JComponent) btnChoose.getParent();
            if (panel == null) return;
            panel.remove(lblColor);
            panel.remove(btnChoose);
            panel.remove(txtColor);
            panel.remove(txtSymbol);
            panel.remove(chkEnable);
        }
    }

    private ArrayList<SymbolColorRule> symbolColorRules = new ArrayList<>();
    private void addSymbolColorRule(final String text, final Color c) {
        SymbolColorRule r = new SymbolColorRule(this);
        r.setColor(c);
        r.setSymbol(text);
        symbolColorRules.add(r);
        r.addToPanel(pnlSymbolColors, symbolColorRules.size());
    }
    private ArrayList<ValueColorRule> valueColorRules = new ArrayList<>();
    private void addValueColorRule(final String value, final Color c) {
        ValueColorRule r = new ValueColorRule(this);
        r.setValue(value);
        r.setColor(c);
        valueColorRules.add(r);
        r.addToPanel(pnlValueColors, valueColorRules.size());
    }
    private ArrayList<ValueColorRule> probColorRules = new ArrayList<>();
    private void addProbColorRule(final String value, final Color c) {
        ValueColorRule r = new ValueColorRule(this);
        r.setValue(value);
        r.setColor(c);
        probColorRules.add(r);
        r.addToPanel(pnlProbColors, probColorRules.size());
    }
    private void invokeScrollToComponent(final Component c) {
        SwingUtilities.invokeLater(() -> ((JComponent) c.getParent()).scrollRectToVisible(c.getBounds()));
    }

    private boolean removeLastRule(List<? extends SimpleColorControl> list) {
        if (list.size() <= 1) return false;
        SimpleColorControl c = list.get(list.size() - 1);
        JComponent panel = c.getPanel();
        c.remove();
        list.remove(list.size() - 1);
        panel.revalidate();
        return true;
    }

    @Override
    public void actionPerformed(final ActionEvent e) {
        switch (e.getActionCommand()) {
            case "addSymbolRule":
                addSymbolColorRule("", Color.GRAY);
                pnlSymbolColors.revalidate();
                invokeScrollToComponent(btnAddSymbolRule);
                break;
            case "addValueRule":
                addValueColorRule("", Color.GRAY);
                pnlValueColors.revalidate();
                invokeScrollToComponent(pnlValueColors.getComponent(pnlValueColors.getComponentCount() - 1));
                break;
            case "addProbRule":
                addProbColorRule("", Color.GRAY);
                pnlProbColors.revalidate();
                invokeScrollToComponent(pnlProbColors.getComponent(pnlProbColors.getComponentCount() - 1));
                break;
            case "apply-text":
                applyColors(ColorMode.Text);
                break;
            case "apply-fill":
                applyColors(ColorMode.Fill);
                break;
            case "apply-line":
                applyColors(ColorMode.Outline);
                break;
            case "apply-bond":
                applyColors(ColorMode.Bond);
                break;
            case "reset-text":
                resetExistingBaseColors(ColorMode.Text);
                break;
            case "reset-fill":
                resetExistingBaseColors(ColorMode.Fill);
                break;
            case "reset-line":
                resetExistingBaseColors(ColorMode.Outline);
                break;
            case "reset-bond":
                resetExistingBaseColors(ColorMode.Bond);
                break;
            case "reset-all":
                resetExistingBaseColors(ColorMode.Text | ColorMode.Fill | ColorMode.Outline | ColorMode.Bond);
                break;
            case "close":
                close();
                break;
        }
    }
    private void close() {
        this.dispose();
    }
    private boolean applyColors(int colorMode) {
        saveSettings();

        DrawWindow target = getActiveDrawing(true);
        if (target == null) return false;

        SceneColorizer sc = null;
        //int colorMode = getColorMode();

        switch (getTab()) {
            case "prob": {
                SceneColorizer.NumberRangeList colorRules = readNumericRules(probColorRules, "probability", getProbComparisonSymbol());
                if (colorRules == null) return false;
                SceneColorizer.PairProbabilityData data = parseProbData(txtProbData.getText(), calcMaxNucCount(target.getScenes()));
                if (data == null) return false;
                sc = SceneColorizer.ColorByProbability(data, colorRules, colorMode);
            }
            break;
            case "mod": {
                SceneColorizer.NumberRangeList colorRules = readNumericRules(valueColorRules, "numeric", getValueComparisonSymbol());
                if (colorRules == null) return false;
                double[] data = parseModData(txtModData.getText(), calcMaxNucCount(target.getScenes()));
                if (data == null) return false;
                sc = SceneColorizer.ColorByDataValue(data, colorRules, colorMode);
            }
            break;
            case "base":
                // Create a map that ignores case if chkIgnoreCase is selected. (HashMap is case sensitive)
                Map<String, Color> colorMap = new HashMap<>(); //chkIgnoreCase.isSelected() ? new TreeMap<>(String.CASE_INSENSITIVE_ORDER) : new HashMap<>();
                for (SymbolColorRule r : symbolColorRules) {
                    if (r.isEnabled() && !r.getSymbol().trim().isEmpty() || !r.getColorText().trim().isEmpty()) {
                        for (char c : r.getSymbol().toCharArray())
                            if (!Character.isWhitespace(c))
                                colorMap.put(Character.toString(c), r.getColor());
                    }
                }
                sc = SceneColorizer.ColorBySymbol(colorMap, colorMode);
                break;
            case "single":
                sc = SceneColorizer.SingleColor(singleColor.getColor(), colorMode);
                break;
            default:
                Dialogs.showWarning("Unknown tab: " + getTab());
                break;
        }
        if (sc == null)
            return false;

        target.colorizeScenes(sc, optStructuresAll.isSelected(), optBasesSelected.isSelected());
        return true;
    }

    private int calcMaxNucCount(RnaSceneGroup scenes) {
        // Note that an RnaSceneGroup is *supposed* to only contain structures that have the SAME sequence (and so the SAME number of nucs.)
        // but just in case, let's get the maximum count.
        int nucs = 0;
        for (RnaScene s : scenes)
            nucs = Math.max(nucs, s.getNucCount());
        return nucs;
    }
    private SceneColorizer.NumberRangeList readNumericRules(List<ValueColorRule> rules, String ruleType, String comparison) {
        double[] limits = new double[rules.size()];
        Color[] colors = new Color[rules.size()];

        int count = 0;
        for (ValueColorRule r : rules) {
            if (r.getValue().trim().isEmpty() || r.getColorText().trim().isEmpty())
                continue;
            limits[count] = parseDouble(r.getValue());
            colors[count] = r.getColor();
            if (Double.isNaN(limits[count])) {
                Dialogs.showWarning("The " + ruleType + " limit in rule %d is not a valid number. Please correct this to continue.");
                r.txtValue.grabFocus();
                return null;
            }
            count++;
        }
        // resize arrays if some rules were not used.
        if (count < limits.length) {
            limits = Arrays.copyOf(limits, count);
            colors = Arrays.copyOf(colors, count);
        }

        SceneColorizer.NumberRangeList ranges = new SceneColorizer.NumberRangeList(limits, colors, comparison, null);

        // We should sort the rules in case the user hasn't listed them correctly.
        // For < and <= comparisons, we should sort the limits in ascending order (from least to greatest) so that
        //   when a data value is analyzed, it will match the smallest limit that it is less than (or equal to).
        // Conversely, for > and >= comparisons, we should sort the values in descending order so that
        //   a data value will only match the largest value that is is greater than (or equal to).
        ranges.sortForComparison();

        return ranges;
    }

    // TODO: try to combine with parseModData
    private SceneColorizer.PairProbabilityData parseProbData(String text, int nucCount) {
        // parse each line of the data and store the values in the array.
        // Each line must be in the form <BASE1><SPACES><BASE2><SPACES><PROB>
        // Where <SPACES> can be any number of spaces or tabs.
        // In addition, spaces or tabs before or after the content is ignored and any line starting with # or ; is ignored as well.
        String[] lines = text.split("\r?\n");
        Pattern values = Pattern.compile("(\\d+)\\s+(\\d+)\\s+([-\\d.,+eEdD]+)");  // <DIGITS><SPACES><DIGITS><SPACES><FLOAT_CHARS>
        String comments = ";#!/";
        String error = null;

        SceneColorizer.PairProbabilityData data = new SceneColorizer.PairProbabilityData(nucCount);

        // Data is listed in the double[] in upper-triangular order. For N=5 (number of nucleobases) the array looks like this:
        // P(1:2) P(1:3) P(1:4) P(1:5)
        // P(2:3) P(2:4) P(2:5)
        // P(3:4) P(3:5)
        // P(4:5)
        // (i.e. 5*4/2 = 10 elements)
        // So to calculate the position of P(I, J)  (where I < J -- if not, swap I, J)
        // Position of first element in row I (0-based) is:
        // ROW(I) = TRIANGLE(N) - TRIANGLE(N-I)   (where TRIANGLE(X) is the number of elements in triangular array -- i.e. X*(X-1)/2
        //        = (N-1)N/2    - (N-I-1)(N-I)/2
        //        = (N-1)N/2    - (N-I-1)(N-I)/2
        //        = (  N^2 - N  -  (N^2-2NI+I^2 +I -N)  ) / 2
        //        = (  2NI - I^2 - I )/2
        //        = (2N-I-1)I/2
        // The first element in row I is P(I, I+1) so the position of P(I, J) is ROW(I) + J - (I+1)
        // so P(I, J) (where I < J) is at (2*N-I-1)*I/2+J-I-1    (for 0-based I and J)
        //                                (2*N-I)*(I-1)/2+J-I-1  (for 1-based I and J -- as is the case with user input)

        int lineNum;
        for (lineNum = 0; lineNum < lines.length; lineNum++) {
            String line = lines[lineNum].trim();
            if (line.isEmpty() || comments.indexOf(line.charAt(0)) != -1)
                continue;
            Matcher m = values.matcher(line);
            if (m.matches()) {
                int i1 = parseInt(m.group(1), Integer.MIN_VALUE); // returns MIN_VALUE if invalid.
                int i2 = parseInt(m.group(2), Integer.MIN_VALUE); // returns MIN_VALUE if invalid.
                double value = parseDouble(m.group(3)); // returns NaN if invalid.
                if (i1 == Integer.MIN_VALUE || i1 > 50000) {
                    // 50000 is an arbitrary limit to prevent memory over-run due to malformed data.
                    error = "invalid nucleobase index: " + m.group(1) + " (outside of valid range)";
                } else if (i2 == Integer.MIN_VALUE || i2 > 50000) {
                    error = "invalid nucleobase index: " + m.group(2) + " (outside of valid range)";
                } else if (i1 < 1 || i2 < 1) {
                    error = "invalid nucleobase index: " + Math.min(i1, i2) + " (less than 1)";
                } else if (Double.isNaN(value)) {
                    error = "invalid numeric value: \"" + m.group(3) + "\"";
                } else if (i1 == i2 && Math.abs(value) > 1E-7) {
                    error = "invalid pair: " + i1 + " cannot pair with itself";
                } else if (i1 > nucCount || i2 > nucCount) {
                    continue; // ignore data outside the range that will be used by the current drawing window.
                } else {
                    data.setPair(i1 - 1, i2 - 1, value);
                    continue;
                }
            } else  // pattern failed.
                error = "invalid text (malformed number or wrong data format, etc.)";
            break; // if we didn't continue already, there was an error.
        }
        if (error != null) {
            Dialogs.showWarning("The data table contains an error: " + error + " on line " + (lineNum + 1) + ".\nPlease correct the problem before continuing.\n\n\nYou can click the Help button to see details about the expected data format.");
            data = null;
            try {
                int lineStart = txtProbData.getLineStartOffset(lineNum);
                int lineEnd = txtProbData.getLineEndOffset(lineNum);
                txtProbData.select(lineStart, lineEnd - 1);
                txtProbData.grabFocus();
            } catch (BadLocationException ex) {
                // nevermind
            }
        }
        return data;
    }

    private double[] parseModData(String text, int nucCount) {
        // parse each line of the data and store the values in the array.
        // Each line must be in the form <NUC_INDEX><SPACE><VALUE>
        // Where <SPACE> can be any number of spaces or tabs.
        // In addition, spaces or tabs before or after the content is ignored and any line starting with # or ; is ignored as well.
        String[] lines = text.split("\r?\n");
        Pattern values = Pattern.compile("(\\d+)\\s+([-\\d.,+eEdD]+)"); // <DIGITS><SPACES><FLOAT_CHARS>
        String comments = ";#!/";
        String error = null;

        // Start out with the biggest best-guess for the size. i.e. max(lineCount, nucCount, 256);
        double[] data = new double[nucCount];
        Arrays.fill(data, Double.NaN); // Double.NaN; indicates a nucleobase isn't listed.

        int lineNum;
        for (lineNum = 0; lineNum < lines.length; lineNum++) {
            String line = lines[lineNum].trim();
            if (line.isEmpty() || comments.indexOf(line.charAt(0)) != -1)
                continue;
            Matcher m = values.matcher(line);
            if (m.matches()) {
                int index = parseInt(m.group(1), Integer.MIN_VALUE); // returns MIN_VALUE if invalid.
                double value = parseDouble(m.group(2)); // returns NaN if invalid.
                if (index == Integer.MIN_VALUE || index > 50000) {
                    // 50000 is an arbitrary limit to prevent memory over-run due to malformed data.
                    error = "invalid nucleobase index: " + m.group(1) + " (outside of valid range)";
                } else if (index < 1) {
                    error = "invalid nucleobase index: " + index + " (less than 1)";
                } else if (Double.isNaN(value)) {
                    error = "invalid numeric value: \"" + m.group(2) + "\"";
                } else {
//                    if (index >= data.length) {
//                        // increase the size of the array to hold the largest index.
//                        int len = data.length;
//                        data = Arrays.copyOf(data, Math.max(2 * data.length, index + 256));
//                        // fill remainder with NaN
//                        Arrays.fill(data, len, data.length, Double.NaN);
//                    }
                    if (index < data.length) // if it is outside the range, just ignore it. It will not be used in any of the drawing's structures
                        data[index - 1] = value;
                    continue;
                }
            } else  // pattern failed.
                error = "invalid text (malformed number or wrong data format, etc.)";
            break; // if we didn't continue already, there was an error.
        }
        if (error != null) {
            Dialogs.showWarning("The data table contains an error: " + error + " on line " + (lineNum + 1) + ".\nPlease correct the problem before continuing.\n\n\nYou can click the Help button to see details about the expected data format.");
            data = null;
            try {
                int lineStart = txtModData.getLineStartOffset(lineNum);
                int lineEnd = txtModData.getLineEndOffset(lineNum);
                txtModData.select(lineStart, lineEnd - 1);
                txtModData.grabFocus();
            } catch (BadLocationException ex) {
                // nevermind
            }
        }
        return data;
    }
    private int parseInt(String value, int valueOnError) {
        try {
            return Integer.parseInt(value);
        } catch (NumberFormatException ex) {
            return valueOnError;
        }
    }
    private double parseDouble(String value) {
        value = value.trim().toUpperCase();
        if (value.equals("INF") || value.equals("+INF"))
            return Double.POSITIVE_INFINITY;
        if (value.equals("-INF") || value.equals("- INF") || value.equals("NEGINF"))
            return Double.NEGATIVE_INFINITY;
        try {
            return Double.parseDouble(value);
        } catch (NumberFormatException ex) {
            return Double.NaN;
        }
    }

    private DrawWindow getActiveDrawing(boolean showWarning) {
        DrawWindow target = Program.getInstance().getMainFrame().getActiveDrawing();
        if (showWarning && target == null)
            Dialogs.showWarning("No drawing window is currently active.\nPlease make sure a drawing is open and active, and that no other program windows are covering it.", "No Active Drawing");
        return target;
    }

    private void resetExistingBaseColors(int colorMode) {
        DrawWindow target = getActiveDrawing(true);
        if (target != null)
            target.colorizeScenes(SceneColorizer.ColorRemover(colorMode), optStructuresAll.isSelected(), optBasesSelected.isSelected());
    }

// return one of the SceneColorizer.ColorMode constants
//    private int getColorMode() {
//        int mode = optBaseFill.isSelected() ? SceneColorizer.ColorMode.Fill :
//                optBaseText.isSelected() ? SceneColorizer.ColorMode.Text :
//                        SceneColorizer.ColorMode.Outline;
//        if (chkColorBond.isSelected())
//            mode |= SceneColorizer.ColorMode.Bond;
//        return mode;
//    }

    private String getTab() { return tabs.getSelectedComponent().getName(); }
    private String _lastColorEffect = "apply-text";
    private void saveSettings() {
        try {
            if (prefs == null) return;
            prefs.put("tab", getTab());
            prefs.put("bases", optBasesAll.isSelected() ? "all" : "sel");
            prefs.put("structures", optStructuresAll.isSelected() ? "all" : "current");
            //prefs.put("effect", optBaseFill.isSelected() ? "fill" : optBaseText.isSelected() ? "text" : "line");
            prefs.put("effect", _lastColorEffect);
            switch (getTab()) {
                case "prob": {
                    Preferences colors = prefs.node("probColors");
                    if (colors.childrenNames().length > probColorRules.size() + 5)
                        colors.clear();
                    int i = 0;
                    for (ValueColorRule r : probColorRules) {
                        if (r.getValue().trim().isEmpty() && r.getColorText().trim().isEmpty()) continue; //skip empty
                        colors.put("V" + i, r.getValue());
                        colors.put("C" + i, r.getColorText().trim());
                        i++;
                    }
                    colors.putInt("count", i);
                    writeTempFile(txtProbData.getText(), ".rnaProbData");
                    prefs.put("comparison", getProbComparisonSymbol());
                }
                case "mod": {
                    Preferences colors = prefs.node("valueColors");
                    if (colors.childrenNames().length > valueColorRules.size() + 5)
                        colors.clear();
                    int i = 0;
                    for (ValueColorRule r : valueColorRules) {
                        if (r.getValue().trim().isEmpty() && r.getColorText().trim().isEmpty()) continue; //skip empty
                        colors.put("V" + i, r.getValue());
                        colors.put("C" + i, r.getColorText().trim());
                        i++;
                    }
                    colors.putInt("count", i);
                    writeTempFile(txtModData.getText(), ".rnaModData");
                    prefs.put("comparison", getValueComparisonSymbol());
                }
                break;
                case "base": {
                    Preferences colors = prefs.node("symbolColors");
                    if (colors.childrenNames().length > symbolColorRules.size() + 5)
                        colors.clear();
                    int i = 0;
                    for (SymbolColorRule r : symbolColorRules) {
                        if (r.getSymbol().trim().isEmpty() && r.getColorText().trim().isEmpty()) continue; //skip empty
                        colors.put("S" + i, r.getSymbol());
                        colors.put("C" + i, r.getColorText().trim());
                        colors.putBoolean("U" + i, r.isEnabled());
                        i++;
                    }
                    colors.putInt("count", i);
                    colors.putBoolean("ignoreCase", chkIgnoreCase.isSelected());
                }
                break;
                case "single":
                    prefs.put("singleColor", singleColor.getColorText());
                    break;
            }
        } catch (BackingStoreException ex) {
            ex.printStackTrace();
        }
    }

    private void loadSettings() {
//        try {
        if (prefs == null) return;

        setTab(prefs.get("tab", null));

        //-------------- Effect Options ----------------------
        if (optBasesAll.isEnabled())
            ("all".equals(prefs.get("bases", "all")) ? optBasesAll : optBasesSelected).setSelected(true);
        if (optStructuresAll.isEnabled())
            ("all".equals(prefs.get("structures", "all")) ? optStructuresAll : optStructuresCurrent).setSelected(true);

        _lastColorEffect = prefs.get("effect", _lastColorEffect);
        //("text".equals(effect) ? optBaseText : "fill".equals(effect) ? optBaseFill : optBaseOutline).setSelected(true);

        //-------------- Symbol Colors ----------------------
        Preferences colors = prefs.node("symbolColors");
        int colorCount = colors.getInt("count", symbolColorRules.size());
        for (int i = 0; i < colorCount; i++) {
            if (symbolColorRules.size() <= i)
                addSymbolColorRule("", Color.BLACK);
            SymbolColorRule r = symbolColorRules.get(i);
            r.setSymbol(colors.get("S" + i, r.getSymbol()));
            r.setColorText(colors.get("C" + i, r.getColorText()));
            if (!colors.getBoolean("U" + i, true))
                r.setEnabled(false);
        }
        for (int i = symbolColorRules.size() - 1; i >= colorCount; i--) {
            symbolColorRules.get(i).remove();
            symbolColorRules.remove(i);
        }
        pnlSymbolColors.revalidate();
        chkIgnoreCase.setSelected(colors.getBoolean("ignoreCase", chkIgnoreCase.isSelected()));

        //-------------- Value (Chemical Modification) ----------------------
        colors = prefs.node("valueColors");
        colorCount = colors.getInt("count", valueColorRules.size());
        for (int i = 0; i < colorCount; i++) {
            if (valueColorRules.size() <= i)
                addValueColorRule("", Color.BLACK);
            ValueColorRule r = valueColorRules.get(i);
            r.setValue(colors.get("V" + i, r.getValue()));
            r.setColorText(colors.get("C" + i, r.getColorText()));
        }
        for (int i = valueColorRules.size() - 1; i >= colorCount; i--) {
            valueColorRules.get(i).remove();
            valueColorRules.remove(i);
        }
        pnlValueColors.revalidate();
        String tmpData = readTempFile(".rnaModData");
        if (tmpData != null)
            txtModData.setText(tmpData);
        setValueComparison(colors.get("comparison", getValueComparisonSymbol()));

        //-------------- Probability ----------------------
        colors = prefs.node("probColors");
        colorCount = colors.getInt("count", probColorRules.size());
        for (int i = 0; i < colorCount; i++) {
            if (probColorRules.size() <= i)
                addProbColorRule("", Color.BLACK);
            ValueColorRule r = probColorRules.get(i);
            r.setValue(colors.get("V" + i, r.getValue()));
            r.setColorText(colors.get("C" + i, r.getColorText()));
        }
        for (int i = probColorRules.size() - 1; i >= colorCount; i--) {
            probColorRules.get(i).remove();
            probColorRules.remove(i);
        }
        pnlProbColors.revalidate();
        tmpData = readTempFile(".rnaProbData");
        if (tmpData != null)
            txtProbData.setText(tmpData);
        setProbComparison(colors.get("comparison", getProbComparisonSymbol()));

        //-------------- Single Color ----------------------
        singleColor.setColorText(prefs.get("singleColor", singleColor.getColorText()));

//        } catch (BackingStoreException ex) {
//            ex.printStackTrace();
//        }
    }

    private boolean writeTempFile(String data, String fileName) {
        try {
            Path tmpDir = Paths.get(System.getProperty("java.io.tmpdir"));
            Strings.writeAll(data, new FileOutputStream(tmpDir.resolve(fileName).toFile()));
            return true;
        } catch (IOException ex) {
            System.err.println("Failed to write temporary colorize data file.");
            ex.printStackTrace();
            return false;
        }
    }
    private String readTempFile(String fileName) {
        try {
            Path tmpDir = Paths.get(System.getProperty("java.io.tmpdir"));
            File file = tmpDir.resolve(fileName).toFile();
            if (file.exists())
                return Strings.readAll(new FileInputStream(file));
        } catch (IOException ex) {
            System.err.println("Failed to read temporary colorize data file.");
            ex.printStackTrace();
        }
        return null;
    }

    private void setTab(final String tab) {
        if (tab == null)
            return;
        for (int i = 0; i < tabs.getTabCount(); i++) {
            if (tab.equals(tabs.getComponentAt(i).getName())) {
                tabs.setSelectedIndex(i);
                break;
            }
        }
    }

    private void createUIComponents() {
        // TODO: place custom component creation code here
    }
}
