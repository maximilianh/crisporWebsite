/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructure.backend.RNA;
import ur_rna.RNAstructure.backend.RNABackend;
import ur_rna.RNAstructureUI.RNAstructureBackendCalculator;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.utilities.FileFilters;
import ur_rna.RNAstructureUI.utilities.FileType;
import ur_rna.Utilities.PathTools;
import ur_rna.Utilities.swing.FormValidationException;
import ur_rna.Utilities.swing.SimpleDocumentListener;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.text.JTextComponent;
import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Objects;

import static ur_rna.Utilities.PathTools.isFile;
import static ur_rna.Utilities.Strings.isEmpty;

/**
 * @author Richard M. Watson
 */
public class PredictMultiStructureWindow extends PredictWindowBase {
    private static final long serialVersionUID = 20160122;

    private JPanel contentPanel;
    private JTabbedPane tabControls;
    private JSpinner spnTemperature;
    private JButton btnStartCalc;
    private JTextField txtOutputFile;
    private JButton btnResetOptions;
    private JLabel lblOptionWarning;
    private JButton btnCancelTask;
    private JProgressBar prgTaskProgress;
    private JLabel lblTaskStatus;
    private JPanel pnlTaskStatus;
    private JCheckBox chkShowDescription;
    private JLabel lblToolDescription;
    private JComboBox cmbTurboMode;
    private JPanel pnlOutput;
    private JPanel pnlSequences;
    private JPanel pnlOptions;
    private JPanel pnlTurboOptions;
    private JButton btnResetAll;
    private JButton btnResetTurboOptions;
    private JList<SequenceItem> lstSeqs;
    private JSpinner spnTurboGamma;
    private JSpinner spnTurboIterations;
    private JSpinner spnTurboMaxEnergyDiff;
    private JSpinner spnTurboWindowSize;
    private JSpinner spnTurboMEAGamma;
    private JSpinner spnTurboThreshold;
    private JSpinner spnTurboMinHelix;
    private JSpinner spnTurboPKIterations;
    private JSpinner spnMultiIterations;
    private JSpinner spnMultiDsvChange;
    private JSpinner spnMultiGapPenalty;
    private JSpinner spnMultiAlnSize;
    private JSpinner spnMultiWindowSize;
    private JSpinner spnMultiMaxStructures;
    private JSpinner spnMultiMaxEnergyDiff;
    private JCheckBox chkAllowBaseInserts;
    private JPanel pnlMultiOptions;
    private JButton btnMoveUp;
    private JButton btnMoveDn;
    private JButton btnRemove;
    private JButton btnClearList;
    private JSpinner spnTurboMaxStructures;
    private JLabel lblBaseInserts;
    private JButton btnResetMultiOptions;
    private JPanel pnlSeqInfo;
    private JSpinner spnMultiMaxPairs;
    private JLabel lblProgressPercent;
    private JTextField txtSeqFile = new JTextField(), txtMSeqFile = new JTextField();

    private HashMap<JComponent, Object> defaultOptionValues = new HashMap<>();
    private ArrayList<JComponent> optionsControls = new ArrayList<>();
    private ArrayList<JComponent> settingsControls = new ArrayList<>();
    private final String settingsNodeName = "PredictMulti";

    /**
     * Constructor.
     */
    public PredictMultiStructureWindow() {
        $$$setupUI$$$();
        autoSetFieldNames();
        setCaption("Predict a Structure Common to Multiple Sequences");

        lblBaseInserts.setLabelFor(chkAllowBaseInserts);

        listenForActions(contentPanel, JButton.class, this); // add all buttons as action listeners

        addFileField("seq", txtSeqFile, FileFilters.Sequence, false, "add-seq");
        addFileField("mseq", txtMSeqFile, FileFilters.MultiSequence, false, "add-mseq");
        addFileField("out", txtOutputFile, FileFilters.CT, true, "set-output");

        pnlTaskStatus.setVisible(false);
        setContent(contentPanel);
        setPreferredSize(new Dimension(700, 600));

        lstSeqs.addListSelectionListener(this::onSequenceSelected);
        SimpleDocumentListener.listen(txtOutputFile, e -> updateSequenceInfoFromUI());

        // optionsControls is a list of all controls that we should store the default value for
        // and then check to see if the user has changed it and show a warning if they have.
        // (because options are on different tabs and the user may have forgotten they set them
        // in a previous instance of this form and their values have been restored from that previous time.)
        addOptionControls(pnlMultiOptions);
        addOptionControls(pnlOptions);
        addOptionControls(pnlTurboOptions);

        // settingsControls is a list of all controls that should have their
        // values saved and restored.
        settingsControls.addAll(optionsControls);
        settingsControls.addAll(Arrays.asList(chkShowDescription));

        loadSettingsDefaults();

        for (JComponent c : optionsControls)
            defaultOptionValues.put(c, getControlValue(c));
        handleSettings(false); // load recent settings values.

        btnResetOptions.setEnabled(false);
        btnResetTurboOptions.setEnabled(false);
        btnResetMultiOptions.setEnabled(false);
        lblOptionWarning.setVisible(false);
        btnResetAll.setVisible(false);

        updateFormUI();

        startUiTimer();
    }
    private void onSequenceSelected(final ListSelectionEvent e) {
        if (!e.getValueIsAdjusting()) loadSequenceInfo(lstSeqs.getSelectedValue());
    }
    private void loadSequenceInfo(final SequenceItem seq) {
        try {
            _seqInfoUpdating = true;
            if (seq == null) {
                txtOutputFile.setText("");
            } else {
                txtOutputFile.setText(seq.outCtFile);
            }
            pnlSeqInfo.setEnabled(seq != null);
        } finally {
            _seqInfoUpdating = false;
        }
    }
    private boolean _seqInfoUpdating;
    private void updateSequenceInfoFromUI() {
        if (_seqInfoUpdating) return;
        SequenceItem seq = getSelectedSequence();
        if (seq == null) return;
        seq.outCtFile = txtOutputFile.getText();
    }

    void addOptionControls(JComponent parent) {
        for (Component c : parent.getComponents()) {
            if (!(c instanceof JComponent)) continue;
            if (c instanceof JCheckBox || c instanceof JTextComponent || c instanceof JSpinner || c instanceof JComboBox || c instanceof JRadioButton)
                optionsControls.add((JComponent) c);
            if (((JComponent) c).getComponentCount() != 0)
                addOptionControls((JComponent) c);
        }
    }

    @Override
    public void showWindow() {
        super.showWindow();
        SwingUtilities.invokeLater(this::pack);
    }

    private void createUIComponents() {
        // TODO: place custom component creation code here
        lstSeqs = new JList<>();
    }

    private void loadSettingsDefaults() {
        setSpin(spnTemperature, TEMP_37C, 0, 600.0);
        setSpin(spnMultiIterations, 2, 1, 100);
        setSpin(spnMultiDsvChange, 1, 1, 1000);
        setSpin(spnMultiMaxEnergyDiff, 20, 1, 200);
        setSpin(spnMultiMaxStructures, 20, 1, 200);
        setSpin(spnMultiMaxPairs, -1);
        setSpin(spnMultiWindowSize, 3, 1, 200);
        setSpin(spnMultiAlnSize, 1, 1, 200);
        setSpin(spnMultiGapPenalty, 0.4, 0, 100);
        setSpin(spnTurboGamma, 0.3, 0, 100);
        setSpin(spnTurboMaxEnergyDiff, 50, 1, 500);
        setSpin(spnTurboWindowSize, 5, 1, 500);
        setSpin(spnTurboMaxStructures, 1000, 1, 5000);
        setSpin(spnTurboMEAGamma, 1.0, 0, 100);
        setSpin(spnTurboIterations, 3, 1, 100);
        setSpin(spnTurboPKIterations, 1, 0, 100);
        setSpin(spnTurboMinHelix, 3, 2, 100);
        setSpin(spnTurboThreshold, 0.0, 0, 100.0);
        chkAllowBaseInserts.setSelected(true);
        cmbTurboMode.setSelectedIndex(0);
    }

    @Override
    protected void processCommand(final CommandInfo ci) {
        switch (ci.getCommand()) {
            case "reset-options-all":
                resetOptionsInPanel(contentPanel);
                break;
            case "reset-options-turbo":
                resetOptionsInPanel(pnlTurboOptions);
                break;
            case "reset-options-multi":
                resetOptionsInPanel(pnlMultiOptions);
                break;
            case "reset-options-general":
                resetOptionsInPanel(pnlOptions);
                break;
            case "start":
                if (verifyInput()) runCalc();
                break;
            case "cancel-task":
                BackgroundTask task = getRunningTask();
                if (task != null) task.cancel();
                break;
            case "add-seq":
                addSequenceFile(txtSeqFile.getText().trim(), false);
                txtSeqFile.setText("");
                break;
            case "add-mseq":
                //addSequenceFile(txtMSeqFile.getText().trim(), true);
                txtMSeqFile.setText("");
                break;
            case "set-output":
                break;
            case "list-moveup":
                moveSequence(getSelectedSequence(), -1);
                break;
            case "list-movedn":
                moveSequence(getSelectedSequence(), 1);
                break;
            case "list-remove":
                deleteSequence(getSelectedSequence());
                break;
            case "list-clear":
                seqList.clear();
                refreshSeqList();
                break;
            default:
                super.processCommand(ci);
        }
    }

    @SuppressWarnings("unchecked")
    private void refreshSeqList(boolean reselect) {
        if (reselect)
            refreshSeqList(getSelectedSequence());
        else
            refreshSeqList();
    }
    private void refreshSeqList(SequenceItem selectAfterRefresh) {
        refreshSeqList();
        setSelectedSequence(selectAfterRefresh);
    }
    private void refreshSeqList() {
        lstSeqs.setListData(seqList.toArray(new SequenceItem[seqList.size()]));
    }

    private ArrayList<SequenceItem> seqList = new ArrayList<>();
    private SequenceItem getSelectedSequence() {
        Object sel = lstSeqs.getSelectedValue();
        return sel instanceof SequenceItem ? (SequenceItem) sel : null;
    }
    private void setSelectedSequence(SequenceItem sel) {
        lstSeqs.setSelectedValue(sel, true);
    }
    private void deleteSequence(final SequenceItem seq) {
        int pos = seq == null ? -1 : seqList.indexOf(seq);
        if (pos == -1) return;
        seqList.remove(pos);
        refreshSeqList();
        if (pos < seqList.size())
            setSelectedSequence(seqList.get(pos));
        else if (seqList.size() != 0)
            setSelectedSequence(seqList.get(seqList.size() - 1));
    }
    private void moveSequence(final SequenceItem seq, final int direction) {
        int pos = seq == null ? -1 : seqList.indexOf(seq);
        if (pos == -1) return;
        int swap = pos + direction;
        if (swap < 0 || swap >= seqList.size())
            return;
        seqList.set(pos, seqList.get(swap));
        seqList.set(swap, seq);
        refreshSeqList(seq);
    }
    private void addSequenceFile(final String path, final boolean multiSeq) {
        if (multiSeq) {
            throw new UnsupportedOperationException("Not yet implmented.");
        } else {
            RNA rna = new RNA(path, RNABackend.FILE_SEQ, "rna", false, true);
            if (rna.GetErrorCode() != 0)
                Dialogs.showError(rna.GetFullErrorMessage(), "Sequence File Error");
            else
                addSequence(rna.GetStructure().GetSequenceLabel(), path, ModuleWindow.getOutputFile(path, "ct"));
        }
    }
    private void addSequence(final String label, final String seqFile, final String ctFile) {
        SequenceItem s = new SequenceItem(label, seqFile, ctFile);
        seqList.add(s);
        refreshSeqList(s);
    }
    private void resetOptionsInPanel(JComponent parent) {
        for (JComponent c : optionsControls)
            if (parent.isAncestorOf(c))
                setControlValue(c, defaultOptionValues.get(c));
    }
    /**
     * Method generated by IntelliJ IDEA GUI Designer
     * DO NOT edit this method OR call it in your code!
     *
     * @noinspection ALL
     */
    private void $$$setupUI$$$() {
        createUIComponents();
        contentPanel = new JPanel();
        contentPanel.setLayout(new GridBagLayout());
        contentPanel.setAutoscrolls(false);
        contentPanel.setName("contentPanel");
        final JLabel label1 = new JLabel();
        label1.setFont(new Font(label1.getFont().getName(), Font.BOLD, 14));
        label1.setOpaque(false);
        label1.setText("Predict a Secondary Structure Common to Multiple Sequences");
        GridBagConstraints gbc;
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        contentPanel.add(label1, gbc);
        tabControls = new JTabbedPane();
        tabControls.setName("tabbedPane1");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.gridwidth = 2;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        gbc.insets = new Insets(5, 0, 5, 0);
        contentPanel.add(tabControls, gbc);
        pnlSequences = new JPanel();
        pnlSequences.setLayout(new GridBagLayout());
        tabControls.addTab("Sequences", pnlSequences);
        final JLabel label2 = new JLabel();
        label2.setFont(new Font(label2.getFont().getName(), label2.getFont().getStyle(), 14));
        label2.setText("Load two or more sequences to align/fold.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 3, 0, 3);
        pnlSequences.add(label2, gbc);
        final JButton button1 = new JButton();
        button1.setActionCommand("recent-seq");
        button1.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlSequences.add(button1, gbc);
        final JButton button2 = new JButton();
        button2.setActionCommand("browse-seq");
        button2.setText("Add Sequence File ...");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlSequences.add(button2, gbc);
        final JLabel label3 = new JLabel();
        label3.setText(" Add a sequence from a SEQ, FASTA, or Text file.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 10, 0, 0);
        pnlSequences.add(label3, gbc);
        final JLabel label4 = new JLabel();
        label4.setFont(new Font(label4.getFont().getName(), label4.getFont().getStyle(), 14));
        label4.setText("Loaded Sequences:");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.gridwidth = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 3, 0, 3);
        pnlSequences.add(label4, gbc);
        btnMoveUp = new JButton();
        btnMoveUp.setActionCommand("list-moveup");
        btnMoveUp.setText("Move Up");
        btnMoveUp.setToolTipText("Move the sequence down in the list.");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 3;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlSequences.add(btnMoveUp, gbc);
        btnMoveDn = new JButton();
        btnMoveDn.setActionCommand("list-movedn");
        btnMoveDn.setText("Move Down");
        btnMoveDn.setToolTipText("Move the sequence up in the list.");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 4;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlSequences.add(btnMoveDn, gbc);
        btnRemove = new JButton();
        btnRemove.setActionCommand("list-remove");
        btnRemove.setText("Remove");
        btnRemove.setToolTipText("Remvoe the selected sequence from the list.");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 5;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlSequences.add(btnRemove, gbc);
        final JScrollPane scrollPane1 = new JScrollPane();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.gridwidth = 2;
        gbc.gridheight = 4;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        pnlSequences.add(scrollPane1, gbc);
        lstSeqs.setSelectionMode(0);
        scrollPane1.setViewportView(lstSeqs);
        btnClearList = new JButton();
        btnClearList.setActionCommand("list-clear");
        btnClearList.setText("Clear");
        btnClearList.setToolTipText("Remove ALL sequences from the list (!)");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 6;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlSequences.add(btnClearList, gbc);
        pnlSeqInfo = new JPanel();
        pnlSeqInfo.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 8;
        gbc.gridwidth = 3;
        gbc.fill = GridBagConstraints.BOTH;
        pnlSequences.add(pnlSeqInfo, gbc);
        final JLabel label5 = new JLabel();
        label5.setText("Output File");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        pnlSeqInfo.add(label5, gbc);
        txtOutputFile = new JTextField();
        txtOutputFile.setName("txtOutputFile");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlSeqInfo.add(txtOutputFile, gbc);
        final JButton button3 = new JButton();
        button3.setActionCommand("browse-out");
        button3.setText("...");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlSeqInfo.add(button3, gbc);
        final JButton button4 = new JButton();
        button4.setActionCommand("recent-out");
        button4.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlSeqInfo.add(button4, gbc);
        final JLabel label6 = new JLabel();
        label6.setFont(new Font(label6.getFont().getName(), label6.getFont().getStyle(), 14));
        label6.setText("Options for the Selected Sequence:");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 7;
        gbc.gridwidth = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 3, 0, 3);
        pnlSequences.add(label6, gbc);
        pnlOptions = new JPanel();
        pnlOptions.setLayout(new GridBagLayout());
        tabControls.addTab("General Options", pnlOptions);
        final JLabel label7 = new JLabel();
        label7.setText("Temperature");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        pnlOptions.add(label7, gbc);
        spnTemperature = new JSpinner();
        spnTemperature.setName("");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOptions.add(spnTemperature, gbc);
        final JPanel spacer1 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlOptions.add(spacer1, gbc);
        btnResetOptions = new JButton();
        btnResetOptions.setActionCommand("reset-options-general");
        btnResetOptions.setText("Reset All to Defaults");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.EAST;
        pnlOptions.add(btnResetOptions, gbc);
        final JLabel label8 = new JLabel();
        label8.setText("K");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        pnlOptions.add(label8, gbc);
        final JPanel spacer2 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.VERTICAL;
        pnlOptions.add(spacer2, gbc);
        pnlTurboOptions = new JPanel();
        pnlTurboOptions.setLayout(new GridBagLayout());
        tabControls.addTab("TurboFold Options", pnlTurboOptions);
        final JLabel label9 = new JLabel();
        label9.setText("TurboFold Mode");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        pnlTurboOptions.add(label9, gbc);
        cmbTurboMode = new JComboBox();
        final DefaultComboBoxModel defaultComboBoxModel1 = new DefaultComboBoxModel();
        defaultComboBoxModel1.addElement("Maximum Expected Accuracy");
        defaultComboBoxModel1.addElement("Pseudoknots");
        defaultComboBoxModel1.addElement("Threshold");
        cmbTurboMode.setModel(defaultComboBoxModel1);
        cmbTurboMode.setName("cmbPartsMode");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.gridwidth = 2;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlTurboOptions.add(cmbTurboMode, gbc);
        final JPanel spacer3 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 5;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.VERTICAL;
        pnlTurboOptions.add(spacer3, gbc);
        final JLabel label10 = new JLabel();
        label10.setText("TurboFold Gamma");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 5, 0, 2);
        pnlTurboOptions.add(label10, gbc);
        spnTurboGamma = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlTurboOptions.add(spnTurboGamma, gbc);
        btnResetTurboOptions = new JButton();
        btnResetTurboOptions.setActionCommand("reset-options-turbo");
        btnResetTurboOptions.setText("Reset All to Defaults");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.EAST;
        pnlTurboOptions.add(btnResetTurboOptions, gbc);
        final JPanel panel1 = new JPanel();
        panel1.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.gridwidth = 2;
        gbc.gridheight = 2;
        gbc.weightx = 1.5;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(10, 0, 0, 10);
        pnlTurboOptions.add(panel1, gbc);
        panel1.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Options for Maximum Expected Accuracy Mode:", TitledBorder.DEFAULT_JUSTIFICATION, TitledBorder.DEFAULT_POSITION, new Font(panel1.getFont().getName(), Font.BOLD, panel1.getFont().getSize())));
        final JLabel label11 = new JLabel();
        label11.setText("Max % Energy Difference");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 10, 0, 2);
        panel1.add(label11, gbc);
        final JLabel label12 = new JLabel();
        label12.setText("Window Size");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 10, 0, 2);
        panel1.add(label12, gbc);
        spnTurboMaxEnergyDiff = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel1.add(spnTurboMaxEnergyDiff, gbc);
        spnTurboWindowSize = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel1.add(spnTurboWindowSize, gbc);
        final JLabel label13 = new JLabel();
        label13.setText("Max Number of Structures");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 10, 0, 2);
        panel1.add(label13, gbc);
        spnTurboMaxStructures = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel1.add(spnTurboMaxStructures, gbc);
        final JLabel label14 = new JLabel();
        label14.setText("MEA-Gamma");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 10, 0, 2);
        panel1.add(label14, gbc);
        spnTurboMEAGamma = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 3;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel1.add(spnTurboMEAGamma, gbc);
        final JLabel label15 = new JLabel();
        label15.setText("TurboFold Iterations");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 5, 0, 2);
        pnlTurboOptions.add(label15, gbc);
        spnTurboIterations = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlTurboOptions.add(spnTurboIterations, gbc);
        final JPanel panel2 = new JPanel();
        panel2.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 3;
        gbc.gridwidth = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(10, 5, 0, 0);
        pnlTurboOptions.add(panel2, gbc);
        panel2.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Options for Pseudoknot Mode:", TitledBorder.DEFAULT_JUSTIFICATION, TitledBorder.DEFAULT_POSITION, new Font(panel2.getFont().getName(), Font.BOLD, panel2.getFont().getSize())));
        final JLabel label16 = new JLabel();
        label16.setText("Pseudoknot Iterations");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 10, 0, 2);
        panel2.add(label16, gbc);
        spnTurboPKIterations = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel2.add(spnTurboPKIterations, gbc);
        final JLabel label17 = new JLabel();
        label17.setText("Minimum Helix Length");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.insets = new Insets(0, 10, 0, 2);
        panel2.add(label17, gbc);
        spnTurboMinHelix = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel2.add(spnTurboMinHelix, gbc);
        final JPanel panel3 = new JPanel();
        panel3.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 4;
        gbc.gridwidth = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(10, 5, 0, 0);
        pnlTurboOptions.add(panel3, gbc);
        panel3.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Options for Threshold Mode:", TitledBorder.DEFAULT_JUSTIFICATION, TitledBorder.DEFAULT_POSITION, new Font(panel3.getFont().getName(), Font.BOLD, panel3.getFont().getSize())));
        spnTurboThreshold = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel3.add(spnTurboThreshold, gbc);
        final JLabel label18 = new JLabel();
        label18.setText("Threshold");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 0.5;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 10, 0, 2);
        panel3.add(label18, gbc);
        pnlMultiOptions = new JPanel();
        pnlMultiOptions.setLayout(new GridBagLayout());
        tabControls.addTab("Multilign Options", pnlMultiOptions);
        final JLabel label19 = new JLabel();
        label19.setText("Multilign Iterations");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.insets = new Insets(0, 20, 0, 2);
        pnlMultiOptions.add(label19, gbc);
        spnMultiIterations = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlMultiOptions.add(spnMultiIterations, gbc);
        final JLabel label20 = new JLabel();
        label20.setText("Maximum DSV Change");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.insets = new Insets(0, 20, 0, 2);
        pnlMultiOptions.add(label20, gbc);
        spnMultiDsvChange = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlMultiOptions.add(spnMultiDsvChange, gbc);
        final JLabel label21 = new JLabel();
        label21.setText("Maximum % Energy Difference");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.insets = new Insets(0, 20, 0, 2);
        pnlMultiOptions.add(label21, gbc);
        final JLabel label22 = new JLabel();
        label22.setText("Maximum Number of Structures");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.insets = new Insets(0, 20, 0, 2);
        pnlMultiOptions.add(label22, gbc);
        final JLabel label23 = new JLabel();
        label23.setText("Structure Window Size");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 4;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.insets = new Insets(0, 20, 0, 2);
        pnlMultiOptions.add(label23, gbc);
        final JLabel label24 = new JLabel();
        label24.setText("Alignment Window Size");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 5;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.insets = new Insets(0, 20, 0, 2);
        pnlMultiOptions.add(label24, gbc);
        final JLabel label25 = new JLabel();
        label25.setText("Gap Penalty");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 6;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.insets = new Insets(0, 20, 0, 2);
        pnlMultiOptions.add(label25, gbc);
        lblBaseInserts = new JLabel();
        lblBaseInserts.setText(" Allow Single Base Pair Inserts");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 8;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.insets = new Insets(0, 20, 0, 2);
        pnlMultiOptions.add(lblBaseInserts, gbc);
        spnMultiMaxEnergyDiff = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlMultiOptions.add(spnMultiMaxEnergyDiff, gbc);
        spnMultiMaxStructures = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 3;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlMultiOptions.add(spnMultiMaxStructures, gbc);
        spnMultiWindowSize = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 4;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlMultiOptions.add(spnMultiWindowSize, gbc);
        spnMultiAlnSize = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 5;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlMultiOptions.add(spnMultiAlnSize, gbc);
        spnMultiGapPenalty = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 6;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlMultiOptions.add(spnMultiGapPenalty, gbc);
        final JPanel spacer4 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 9;
        gbc.fill = GridBagConstraints.VERTICAL;
        pnlMultiOptions.add(spacer4, gbc);
        chkAllowBaseInserts = new JCheckBox();
        chkAllowBaseInserts.setText("");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 8;
        gbc.anchor = GridBagConstraints.WEST;
        pnlMultiOptions.add(chkAllowBaseInserts, gbc);
        final JLabel label26 = new JLabel();
        label26.setText("Maximum Pairs");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 7;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.insets = new Insets(0, 20, 0, 2);
        pnlMultiOptions.add(label26, gbc);
        spnMultiMaxPairs = new JSpinner();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 7;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlMultiOptions.add(spnMultiMaxPairs, gbc);
        btnResetMultiOptions = new JButton();
        btnResetMultiOptions.setActionCommand("reset-options-multi");
        btnResetMultiOptions.setText("Reset All to Defaults");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.EAST;
        pnlMultiOptions.add(btnResetMultiOptions, gbc);
        final JLabel label27 = new JLabel();
        label27.setText("(-1 means use the default calculated value)");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 7;
        gbc.weightx = 0.5;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 3, 0, 2);
        pnlMultiOptions.add(label27, gbc);
        pnlOutput = new JPanel();
        pnlOutput.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 4;
        gbc.gridwidth = 2;
        gbc.fill = GridBagConstraints.BOTH;
        contentPanel.add(pnlOutput, gbc);
        lblOptionWarning = new JLabel();
        lblOptionWarning.setBackground(new Color(-2171392));
        lblOptionWarning.setForeground(new Color(-16777216));
        lblOptionWarning.setHorizontalAlignment(0);
        lblOptionWarning.setHorizontalTextPosition(0);
        lblOptionWarning.setOpaque(true);
        lblOptionWarning.setText("<html><center>Warning: Some options have been set to non-default values.<br>---List---");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 3);
        pnlOutput.add(lblOptionWarning, gbc);
        btnStartCalc = new JButton();
        btnStartCalc.setActionCommand("start");
        btnStartCalc.setHorizontalTextPosition(0);
        btnStartCalc.setMinimumSize(new Dimension(54, 40));
        btnStartCalc.setName("btnStartCalc");
        btnStartCalc.setOpaque(false);
        btnStartCalc.setPreferredSize(new Dimension(54, 40));
        btnStartCalc.setText("Start");
        btnStartCalc.setMnemonic('S');
        btnStartCalc.setDisplayedMnemonicIndex(0);
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlOutput.add(btnStartCalc, gbc);
        btnResetAll = new JButton();
        btnResetAll.setActionCommand("reset-options-all");
        btnResetAll.setText("Reset All Options");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.EAST;
        pnlOutput.add(btnResetAll, gbc);
        pnlTaskStatus = new JPanel();
        pnlTaskStatus.setLayout(new GridBagLayout());
        pnlTaskStatus.setVisible(true);
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 5;
        gbc.gridwidth = 2;
        gbc.fill = GridBagConstraints.BOTH;
        contentPanel.add(pnlTaskStatus, gbc);
        btnCancelTask = new JButton();
        btnCancelTask.setActionCommand("cancel-task");
        btnCancelTask.setText("Cancel");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.BOTH;
        pnlTaskStatus.add(btnCancelTask, gbc);
        prgTaskProgress = new JProgressBar();
        prgTaskProgress.setMinimumSize(new Dimension(10, 40));
        prgTaskProgress.setPreferredSize(new Dimension(146, 20));
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlTaskStatus.add(prgTaskProgress, gbc);
        lblTaskStatus = new JLabel();
        lblTaskStatus.setText("Status of background Task");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        pnlTaskStatus.add(lblTaskStatus, gbc);
        lblProgressPercent = new JLabel();
        lblProgressPercent.setHorizontalAlignment(0);
        lblProgressPercent.setPreferredSize(new Dimension(40, 16));
        lblProgressPercent.setText("0%");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        pnlTaskStatus.add(lblProgressPercent, gbc);
        chkShowDescription = new JCheckBox();
        chkShowDescription.setActionCommand("");
        chkShowDescription.setSelected(true);
        chkShowDescription.setText("Show Description");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHEAST;
        contentPanel.add(chkShowDescription, gbc);
        lblToolDescription = new JLabel();
        lblToolDescription.setName("lblToolDescription");
        lblToolDescription.setText("<html>This tool takes two or more sequences and folds them into their common lowest free energy conformations. It combines the capabilities of <b>Multilign</b> and <b>TurboFold</b> to create distinct sets of possible structures for multiple sequences.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.gridwidth = 2;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        contentPanel.add(lblToolDescription, gbc);
        final JSeparator separator1 = new JSeparator();
        separator1.setMinimumSize(new Dimension(1, 3));
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.gridwidth = 2;
        gbc.fill = GridBagConstraints.BOTH;
        contentPanel.add(separator1, gbc);
        label5.setLabelFor(txtOutputFile);
        label7.setLabelFor(spnTemperature);
        label9.setLabelFor(cmbTurboMode);
        label10.setLabelFor(spnTurboGamma);
        label11.setLabelFor(spnTurboMaxEnergyDiff);
        label12.setLabelFor(spnTurboWindowSize);
        label13.setLabelFor(spnTurboMaxStructures);
        label14.setLabelFor(spnTurboMEAGamma);
        label15.setLabelFor(spnTurboIterations);
        label16.setLabelFor(spnTurboPKIterations);
        label17.setLabelFor(spnTurboMinHelix);
        label18.setLabelFor(spnTurboThreshold);
        label19.setLabelFor(spnMultiIterations);
        label20.setLabelFor(spnMultiDsvChange);
        label21.setLabelFor(spnMultiMaxEnergyDiff);
        label22.setLabelFor(spnMultiMaxStructures);
        label23.setLabelFor(spnMultiWindowSize);
        label24.setLabelFor(spnMultiAlnSize);
        label25.setLabelFor(spnMultiGapPenalty);
        label26.setLabelFor(spnMultiMaxPairs);
    }
    /** @noinspection ALL */
    public JComponent $$$getRootComponent$$$() { return contentPanel; }

    private class SequenceItem {
        public String label;
        public String inSeqFile;
        public String outCtFile;
        @Override
        public String toString() {
            return label + " (" + inSeqFile + ")";
        }
        public SequenceItem(final String label, final String seqFile, final String ctFile) {
            this.label = label;
            inSeqFile = seqFile;
            outCtFile = ctFile;
        }
    }

    @Override
    protected void updateFormUI() {
        // If there is a running task, just update progress and return. Do not modify any form controls.
        PredictTask task = (PredictTask) getRunningTask();
        if (task != null) {
            setTaskStatus(task.getProgress(), task.getStatus());
            btnCancelTask.setEnabled(task.canCancel);
            return;
        }

        boolean modified = false;
        StringBuilder diffs = new StringBuilder();
        for (JComponent c : optionsControls) {
            Object value = getControlValue(c);
            if (!Objects.equals(value, defaultOptionValues.get(c))) {
                diffs.append(getControlDescription(c)).append(": ").append(formatValue(value)).append("    ");
                modified = true;
            }
        }

        if (modified) {
            String warn = lblOptionWarning.getText();
            int end = warn.indexOf("<br>");
            String newWarn = warn.substring(0, end + 4) + diffs;
            if (!newWarn.equals(warn))
                lblOptionWarning.setText(newWarn);
        }
        lblOptionWarning.setVisible(modified);
        btnResetAll.setVisible(modified);

        // Only enable controls when there is no running tasks.
        btnResetOptions.setEnabled(modified);
        btnResetTurboOptions.setEnabled(modified);
        btnResetMultiOptions.setEnabled(modified);

        if (lblToolDescription.isVisible() != chkShowDescription.isSelected()) {
            lblToolDescription.setVisible(chkShowDescription.isSelected());
            if (this.isVisible())
                pack();
        }
    }

    private void runCalc() {
        //TODO:  saveRecentFiles(txtConFile1, txtConFile2, txtSeqFile1, txtSeqFile2, txtOutputFile, txtOutputFile2);
        handleSettings(true);
        btnStartCalc.setVisible(false);
        pnlTaskStatus.setVisible(true);
        enableDescendants(tabControls, false);
        enableDescendants(pnlOutput, false);
        setTaskStatus(0, "Starting calculation...");

        // set variables based on form input
        PredictTask task = new PredictTask();
        task.turboMode = cmbTurboMode.getSelectedIndex();

        task.multiIterations = getSpinInt(spnMultiIterations);
        task.multiDsvChange = getSpinInt(spnMultiDsvChange);
        task.multiMaxEnergyDiff = getSpinInt(spnMultiMaxEnergyDiff);
        task.multiMaxStructures = getSpinInt(spnMultiMaxStructures);
        task.multiMaxPairs = getSpinInt(spnMultiMaxPairs);
        task.multiWindowSize = getSpinInt(spnMultiWindowSize);
        task.multiAlnSize = getSpinInt(spnMultiAlnSize);
        task.multiGapPenalty = getSpinDbl(spnMultiGapPenalty);
        task.allowBaseInserts = chkAllowBaseInserts.isSelected();
        task.temperature = getSpinDbl(spnTemperature);
        task.turboGamma = getSpinDbl(spnTurboGamma);
        task.turboMaxEnergyDiff = getSpinInt(spnTurboMaxEnergyDiff);
        task.turboWindowSize = getSpinInt(spnTurboWindowSize);
        task.turboMaxStructures = getSpinInt(spnTurboMaxStructures);
        task.turboMEAGamma = getSpinDbl(spnTurboMEAGamma);
        task.turboIterations = getSpinInt(spnTurboIterations);
        task.turboPKIterations = getSpinInt(spnTurboPKIterations);
        task.turboMinHelix = getSpinInt(spnTurboMinHelix);
        task.turboThreshold = getSpinDbl(spnTurboThreshold);

        runInBackground(task, this::calcComplete);
    }

    private void calcComplete(PredictTask task) {
        btnStartCalc.setVisible(true);
        pnlTaskStatus.setVisible(false);
        enableDescendants(tabControls, true);
        enableDescendants(pnlOutput, true);
        if (task.isCanceled)
            return; // return immediately with no message or other action, because this could be in response to the form closing.

        // Give the UI time to catch up.
        if (task.onComplete()) this.dispose();
    }

    private class PredictTask extends PredictTaskBase {
        RNAstructureBackendCalculator bc;
        int turboMode;
        int multiIterations;
        int multiDsvChange;
        int multiMaxPairs;
        int multiMaxEnergyDiff;
        int multiMaxStructures;
        int multiWindowSize;
        int multiAlnSize;
        double multiGapPenalty;
        boolean allowBaseInserts;
        double temperature;
        double turboGamma;
        int turboMaxEnergyDiff;
        int turboWindowSize;
        int turboMaxStructures;
        double turboMEAGamma;
        int turboIterations;
        int turboPKIterations;
        int turboMinHelix;
        double turboThreshold;

        // Generated output file names
        String[] turboCTs, turboPFSs, multiCTs, multiPFSs;
        String turboAln, multiAln;

        @Override
        public int getProgress() {
            if (bc == null) return workDone;
            return workDone + bc.getProgressNumber() * nextStepWork / 100;
        }
        @Override
        protected void runCalc() {
            final int count = seqList.size();
            // Generate output file names.
            turboCTs = new String[count];
            turboPFSs = new String[count];
            multiCTs = new String[count];
            multiPFSs = new String[count];
            StringBuilder combinedBase = new StringBuilder();
            combinedBase.append(PathTools.getDir(seqList.get(0).outCtFile, true)); // use output directory from first output CT file.

            bc = new RNAstructureBackendCalculator();
            bc.activateTurboFold();
            for (int i = 0; i < count; i++) {
                SequenceItem seq = seqList.get(i);
                String base = PathTools.changeExtension(seq.outCtFile, null); // remove the extension.
                turboCTs[i] = base + "_Turbo.ct";
                turboPFSs[i] = base + "_Turbo.pfs";
                multiCTs[i] = seq.outCtFile;
                bc.addTurboFoldTuple(seq.inSeqFile, turboCTs[i], turboPFSs[i]);
                combinedBase.append(PathTools.getBaseName(seq.outCtFile)).append("_");
            }
            combinedBase.setLength(combinedBase.length() - 1); // remove final slash
            turboAln = combinedBase + "_Turbo.aln";
            multiAln = combinedBase + ".aln";

            nextStep("Running TurboFold...", false, 25);
            String backendError;
            switch (turboMode) {
                case 0: // MEA
                    backendError = bc.runTurboFoldMaximumExpectedAccuracy(turboGamma, turboIterations, turboMaxEnergyDiff, turboMaxStructures, turboWindowSize, turboMEAGamma, turboAln);
                    break;
                case 1: // PK
                    backendError = bc.runTurboFoldPseudoknot(turboGamma, turboIterations, turboPKIterations, turboMinHelix, turboAln);
                    break;
                default: // Threshold
                    backendError = bc.runTurboFoldThreshold(turboGamma, turboIterations, turboThreshold, turboAln);
                    break;
            }
            if (!isEmpty(backendError)) {
                setError(backendError);
                return;
            }
            bc = new RNAstructureBackendCalculator();
            nextStep("Running Multilign...", false, 75);
            bc.activateMultilign();
            for (int i = 0; i < count; i++)
                bc.addMultilignTuple(seqList.get(i).inSeqFile, multiCTs[i]);
            int maxPairs = multiMaxPairs == -1 ? bc.getMultilignMaxPairs() : multiMaxPairs;
            backendError = bc.runMultilign(multiMaxEnergyDiff, multiMaxStructures, multiWindowSize, multiAlnSize, multiGapPenalty, allowBaseInserts, multiDsvChange, maxPairs, multiIterations, multiAln, false, true);
            if (!isEmpty(backendError)) {
                setError(backendError);
                //return;
            }
        }

        /**
         * Show the results window.
         *
         * @return True if the results window could be shown. False otherwise.
         */
        protected boolean showResults() {
            // else show results window.
            PredictionResultsWindow w = new PredictionResultsWindow("Results of Multiple Structure Prediction");
            final int count = seqList.size();
            for (int i = 0; i < count; i++)
                if (isFile(turboCTs[i])) w.addFile(turboCTs[i], "TurboFold Structures from Sequence " + (i + 1));
            for (int i = 0; i < count; i++)
                if (isFile(multiCTs[i])) w.addFile(multiCTs[i], "Multilign Structures from Sequence " + (i + 1));
            for (int i = 0; i < count; i++)
                if (isFile(turboPFSs[i]))
                    w.addFile(turboPFSs[i], "TurboFold Partition File (PFS) from Sequence " + (i + 1));

            if (isFile(turboAln)) w.addFile(turboAln, "TurboFold Alignment File");
            if (isFile(multiAln)) w.addFile(multiAln, "Multilign Alignment File");

            w.addPlotHeader("TurboFold Structures (without color annotations)", count);
            for (int i = 0; i < count; i++) {
                final String ctFile = turboCTs[i];
                if (isFile(ctFile))
                    w.addPlot("Sequence " + (i + 1), "Aligned structures from Sequence " + (i + 1),
                            () -> launchDrawing(ctFile, FileType.CT, null, null, 1));
            }
            w.addPlotHeader("TurboFold Structures, Color-Annotated by Base-Pairing Probability", count);
            for (int i = 0; i < count; i++) {
                final String ctFile = turboCTs[i], pfsFile = turboPFSs[i];
                if (isFile(ctFile) && isFile(pfsFile))
                    w.addPlot("Sequence " + (i + 1), "Aligned structures from Sequence " + (i + 1),
                            () -> launchDrawing(ctFile, FileType.CT, pfsFile, FileType.PFS, 1));
            }
            w.addPlotHeader("Multilign Structures (without color annotations)", count);
            for (int i = 0; i < count; i++) {
                final String ctFile = multiCTs[i];
                if (isFile(ctFile))
                    w.addPlot("Sequence " + (i + 1), "Aligned structures from Sequence " + (i + 1),
                            () -> launchDrawing(ctFile, FileType.CT, null, null, 1));
            }
            w.showAbbreviations(false);
            w.showWindow();
            return true;
        }
    }

    private void setTaskStatus(final int progress, final String status) {
        lblTaskStatus.setText(status);
        prgTaskProgress.setIndeterminate(progress < 0);
        lblProgressPercent.setVisible(progress >= 0);
        int p = progress > 100 ? 100 : progress;
        if (progress >= 0) {
            lblProgressPercent.setText(p + "%");
            prgTaskProgress.setValue(p);
        }
    }

    private boolean verifyInput() {
        try {
            if (seqList.size() < 2)
                validateFailed(lstSeqs, "At least two sequences are required.");
            int i = 0;
            for (SequenceItem seq : seqList) {
                ++i;
                if (!isFile(seq.inSeqFile)) {
                    setSelectedSequence(seq);
                    validateFailed(lstSeqs, "Sequence file %s does not exist.", i);
                }
            }
        } catch (FormValidationException ex) {
            Dialogs.showError(ex.getMessage(), "Error in User Input");
            focusField(ex.getControl());
            return false;
        }
        return true;
    }

    @Override
    protected void handleSettings(boolean writeValues) {
        handleRecentControlValues(writeValues, settingsNodeName, settingsControls);
    }
}
