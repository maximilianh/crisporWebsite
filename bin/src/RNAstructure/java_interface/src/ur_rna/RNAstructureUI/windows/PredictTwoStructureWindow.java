/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructure.backend.Dynalign_object;
import ur_rna.RNAstructure.backend.RNA;
import ur_rna.RNAstructure.backend.RNABackend;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.utilities.FileFilters;
import ur_rna.RNAstructureUI.utilities.FileType;
import ur_rna.Utilities.PathTools;
import ur_rna.Utilities.swing.FormValidationException;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Objects;

import static ur_rna.Utilities.PathTools.isFile;
import static ur_rna.Utilities.Strings.isWhiteSpace;

/**
 * @author Richard M. Watson
 */
public class PredictTwoStructureWindow extends PredictWindowBase {
    private static final long serialVersionUID = 20160122;

    private JPanel contentPanel;
    private JTextField txtSeqFile1;
    private JTabbedPane tabControls;
    private JSpinner spnMaxEnergyDiff;
    private JSpinner spnWindowSize;
    private JSpinner spnMaxStructs;
    private JTextField txtConFile1;
    private JButton btnStartCalc;
    private JTextField txtOutputFile1;
    private JTextField txtOutputFile2;
    private JButton btnResetOptions;
    private JLabel lblOptionWarning;
    private JButton btnCancelTask;
    private JProgressBar prgTaskProgress;
    private JLabel lblTaskStatus;
    private JPanel pnlTaskStatus;
    private JCheckBox chkShowDescription;
    private JLabel lblToolDescription;
    private JTextField txtSeqFile2;
    private JSpinner spnAlnWindowSize;
    private JSpinner spnGapPenalty;
    private JCheckBox chkBaseInserts;
    private JComboBox cmbPartsMode;
    private JLabel lblBaseInserts;
    private JSpinner spnStochasticSamples;
    private JTextField txtConFile2;
    private JPanel pnlOutput;
    private JPanel pnlSequences;
    private JPanel pnlOptions;
    private JPanel pnlPartsOptions;
    private JPanel pnlConstraints;
    private JButton btnResetAll;
    private JButton btnResetPartsOptions;
    private JLabel lblProgressPercent;

    private HashMap<JComponent, Object> defaultOptionValues = new HashMap<>();
    private ArrayList<JComponent> optionsControls = new ArrayList<>();
    private ArrayList<JComponent> settingsControls = new ArrayList<>();
    private final String settingsNodeName = "PredictDynalign";

    /**
     * Constructor.
     */
    public PredictTwoStructureWindow() {
        $$$setupUI$$$();
        autoSetFieldNames();
        setCaption("Predict a Structure Common to Two Sequences");

        lblBaseInserts.setLabelFor(chkBaseInserts);

        listenForActions(contentPanel, JButton.class, this); // add all buttons as action listeners

//        listenForActions(, btnBrowseCon2, btnBrowseOutput, btnBrowseSeq1, btnBrowseSeq2, btnBrowseShape,
//                btnStartCalc, btnResetOptions, btnCancelTask);

        addFileField("seq1", txtSeqFile1, FileFilters.Sequence, false);
        addFileField("seq2", txtSeqFile2, FileFilters.Sequence, false);
        addFileField("con1", txtConFile1, FileFilters.Constraints, false);
        addFileField("con2", txtConFile2, FileFilters.Constraints, false);
        addFileField("out1", txtOutputFile1, FileFilters.CT, true);
        addFileField("out2", txtOutputFile2, FileFilters.CT, true);

        pnlTaskStatus.setVisible(false);
        setContent(contentPanel);
        setPreferredSize(new Dimension(700, 600));

        // optionsControls is a list of all controls that we should store the default value for
        // and then check to see if the user has changed it and show a warning if they have.
        // (because options are on different tabs and the user may have forgotten they set them
        // in a previous instance of this form and their values have been restored from that previous time.)
        optionsControls.addAll(Arrays.asList(
                spnMaxEnergyDiff,
                spnMaxStructs, spnWindowSize, spnAlnWindowSize,
                spnGapPenalty, chkBaseInserts, spnStochasticSamples, cmbPartsMode
        ));

        // settingsControls is a list of all controls that should have their
        // values saved and restored.
        settingsControls.addAll(optionsControls);
        settingsControls.addAll(Arrays.asList(chkShowDescription));

        loadSettingsDefaults();

        for (JComponent c : optionsControls)
            defaultOptionValues.put(c, getControlValue(c));
        handleSettings(false); // load recent settings values.

        btnResetOptions.setEnabled(false);
        btnResetPartsOptions.setEnabled(false);
        lblOptionWarning.setVisible(false);
        btnResetAll.setVisible(false);

        updateFormUI();

        startUiTimer();
    }

    @Override
    public void showWindow() {
        super.showWindow();
        SwingUtilities.invokeLater(this::pack);
    }

    private void createUIComponents() {
        // TODO: place custom component creation code here
    }

    private void loadSettingsDefaults() {
        setSpin(spnMaxEnergyDiff, 20, 1, 200);
        setSpin(spnMaxStructs, 20, 1, 1000);
        setSpin(spnWindowSize, 3, 1, 100);
        setSpin(spnAlnWindowSize, 1, 1, 100);
        setSpin(spnGapPenalty, 0.4);
        setSpin(spnStochasticSamples, 20, 1, 1000);
        chkBaseInserts.setSelected(true);
        cmbPartsMode.setSelectedIndex(0);
    }

    @Override
    protected void processCommand(final CommandInfo ci) {
        switch (ci.getCommand()) {
            case "reset-options-all":
                for (JComponent c : optionsControls)
                    setControlValue(c, defaultOptionValues.get(c));
                break;
            case "reset-options-dynalign":
                for (JComponent c : optionsControls)
                    if (pnlOptions.isAncestorOf(c))
                        setControlValue(c, defaultOptionValues.get(c));
                break;
            case "reset-options-parts":
                for (JComponent c : optionsControls)
                    if (pnlPartsOptions.isAncestorOf(c))
                        setControlValue(c, defaultOptionValues.get(c));
                break;
            case "start":
                if (verifyInput()) runCalc();
                break;
            case "cancel-task":
                BackgroundTask task = getRunningTask();
                if (task != null) task.cancel();
                break;
            default:
                super.processCommand(ci);
        }
    }

    private String lastSeqFile1 = "", lastOutFile1 = ""; // detect changes to these fields.
    private String lastSeqFile2 = "", lastOutFile2 = ""; // detect changes to these fields.

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
        //boolean useSeqFile = optSeqFile.isSelected();
//        txtSeqFile.setEnabled(useSeqFile);
//        txtSequence.setEnabled(!useSeqFile);
//        txtSeqTitle.setEnabled(!useSeqFile);
//        scrlSequence.setEnabled(!useSeqFile);
        btnResetOptions.setEnabled(modified);
        btnResetPartsOptions.setEnabled(modified);
        btnResetOptions.setEnabled(modified);

        //if (useSeqFile) {
        if (!lastSeqFile1.equals(txtSeqFile1.getText())) {
            lastSeqFile1 = txtSeqFile1.getText();
            // only do the following if the user has NOT changed the output file name themselves
            if (txtOutputFile1.getText().isEmpty() || lastOutFile1.equals(txtOutputFile1.getText()))
                txtOutputFile1.setText(lastOutFile1 = ModuleWindow.getOutputFile(txtSeqFile1.getText(), "ct"));
        }
        if (!lastSeqFile2.equals(txtSeqFile2.getText())) {
            lastSeqFile2 = txtSeqFile2.getText();
            // only do the following if the user has NOT changed the output file name themselves
            if (txtOutputFile2.getText().isEmpty() || lastOutFile2.equals(txtOutputFile2.getText()))
                txtOutputFile2.setText(lastOutFile2 = ModuleWindow.getOutputFile(txtSeqFile2.getText(), "ct"));
        }
        //}

        if (lblToolDescription.isVisible() != chkShowDescription.isSelected()) {
            lblToolDescription.setVisible(chkShowDescription.isSelected());
            if (this.isVisible())
                pack();
        }
    }

    private void runCalc() {
        saveRecentFiles(txtConFile1, txtConFile2, txtSeqFile1, txtSeqFile2, txtOutputFile1, txtOutputFile2);
        handleSettings(true);
        btnStartCalc.setVisible(false);
        pnlTaskStatus.setVisible(true);
        enableDescendants(tabControls, false);
        enableDescendants(pnlOutput, false);
        setTaskStatus(0, "Starting calculation...");

        // set variables based on form input
        PredictTask task = new PredictTask();
//        if (optEnterSeq.isSelected()) {
//            task.sequence = txtSequence.getText();
//            task.sequenceTitle = txtSeqTitle.getText();
//            task.seqType = RNABackend.SEQUENCE_STRING;
//        } else {
        task.sequence1 = txtSeqFile1.getText();
        task.sequence2 = txtSeqFile1.getText();
        task.seqType1 = RNABackend.FILE_SEQ;
        task.seqType2 = RNABackend.FILE_SEQ;
//        }
        //task.alphabet = "rna"; //optDNA.isSelected() ? "dna" : "rna";
        task.baseOutputName1 = PathTools.changeExtension(txtOutputFile1.getText().trim(), "");
        task.baseOutputName2 = PathTools.changeExtension(txtOutputFile2.getText().trim(), "");
        task.baseOutputNameShared = PathTools.getDir(task.baseOutputName1, true) + PathTools.getBaseName(task.baseOutputName1) + '_' + PathTools.getBaseName(task.baseOutputName2);

        // task.temperature = 310.15; //(Double) spnTemperature.getValue();
        // task.maxLoopSize = (Integer) spnMaxLoopSize.getValue();
        task.maxEnergyDiff = (Integer) spnMaxEnergyDiff.getValue();
        task.maxStructures = (Integer) spnMaxStructs.getValue();
        task.windowSize = (Integer) spnWindowSize.getValue();
        task.alnWindowSize = (Integer) spnWindowSize.getValue();
        task.gapPenalty = (Double) spnGapPenalty.getValue();
        task.allowBaseInserts = chkBaseInserts.isSelected();
        //task.shapeSlope = (Double) spnShapeSlope.getValue();
        //task.shapeIntercept = (Double) spnShapeIntercept.getValue();
        //task.shapeFile = txtShapeFile.getText();
        task.conFile1 = txtConFile1.getText();
        task.conFile2 = txtConFile2.getText();

        task.stochasticSamples = (Integer) spnStochasticSamples.getValue();
        task.partsMode = cmbPartsMode.getSelectedIndex();

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
    /**
     * Method generated by IntelliJ IDEA GUI Designer
     * DO NOT edit this method OR call it in your code!
     *
     * @noinspection ALL
     */
    private void $$$setupUI$$$() {
        contentPanel = new JPanel();
        contentPanel.setLayout(new GridBagLayout());
        contentPanel.setAutoscrolls(false);
        contentPanel.setName("contentPanel");
        final JLabel label1 = new JLabel();
        label1.setFont(new Font(label1.getFont().getName(), Font.BOLD, 14));
        label1.setOpaque(false);
        label1.setText("Predict a Secondary Structure Common to Two Sequences");
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
        txtSeqFile1 = new JTextField();
        txtSeqFile1.setName("txtSeqFile1");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlSequences.add(txtSeqFile1, gbc);
        final JButton button1 = new JButton();
        button1.setActionCommand("browse-seq1");
        button1.setText("...");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlSequences.add(button1, gbc);
        final JButton button2 = new JButton();
        button2.setActionCommand("recent-seq1");
        button2.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlSequences.add(button2, gbc);
        final JLabel label2 = new JLabel();
        label2.setText("Sequence File 1");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 3, 0, 3);
        pnlSequences.add(label2, gbc);
        final JLabel label3 = new JLabel();
        label3.setText("Sequence File 2");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 3, 0, 3);
        pnlSequences.add(label3, gbc);
        txtSeqFile2 = new JTextField();
        txtSeqFile2.setName("txtSeqFile2");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlSequences.add(txtSeqFile2, gbc);
        final JButton button3 = new JButton();
        button3.setActionCommand("browse-seq2");
        button3.setText("...");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlSequences.add(button3, gbc);
        final JButton button4 = new JButton();
        button4.setActionCommand("recent-seq2");
        button4.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlSequences.add(button4, gbc);
        final JPanel spacer1 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.VERTICAL;
        pnlSequences.add(spacer1, gbc);
        pnlOptions = new JPanel();
        pnlOptions.setLayout(new GridBagLayout());
        tabControls.addTab("Dynalign Options", pnlOptions);
        final JLabel label4 = new JLabel();
        label4.setText("Maximum % Energy Difference");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        pnlOptions.add(label4, gbc);
        final JLabel label5 = new JLabel();
        label5.setText("Maximum Number of Structures");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        pnlOptions.add(label5, gbc);
        final JLabel label6 = new JLabel();
        label6.setText("Structure Window Size");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        pnlOptions.add(label6, gbc);
        spnMaxEnergyDiff = new JSpinner();
        spnMaxEnergyDiff.setName("spnMaxEnergyDiff");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOptions.add(spnMaxEnergyDiff, gbc);
        spnMaxStructs = new JSpinner();
        spnMaxStructs.setName("spnMaxStructs");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOptions.add(spnMaxStructs, gbc);
        spnWindowSize = new JSpinner();
        spnWindowSize.setName("spnWindowSize");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOptions.add(spnWindowSize, gbc);
        spnGapPenalty = new JSpinner();
        spnGapPenalty.setName("spnGapPenalty");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 4;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOptions.add(spnGapPenalty, gbc);
        final JLabel label7 = new JLabel();
        label7.setText("Gap Penalty");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 4;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        pnlOptions.add(label7, gbc);
        final JLabel label8 = new JLabel();
        label8.setText("Alignment Window Size");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        pnlOptions.add(label8, gbc);
        spnAlnWindowSize = new JSpinner();
        spnAlnWindowSize.setName("spnAlnWindowSize");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 3;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOptions.add(spnAlnWindowSize, gbc);
        final JPanel spacer2 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 4;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlOptions.add(spacer2, gbc);
        btnResetOptions = new JButton();
        btnResetOptions.setActionCommand("reset-options-dynalign");
        btnResetOptions.setText("Reset All to Defaults");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.EAST;
        pnlOptions.add(btnResetOptions, gbc);
        pnlPartsOptions = new JPanel();
        pnlPartsOptions.setLayout(new GridBagLayout());
        tabControls.addTab("PARTS Options", pnlPartsOptions);
        final JLabel label9 = new JLabel();
        label9.setText("PARTS Mode");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        pnlPartsOptions.add(label9, gbc);
        cmbPartsMode = new JComboBox();
        final DefaultComboBoxModel defaultComboBoxModel1 = new DefaultComboBoxModel();
        defaultComboBoxModel1.addElement("MAP (Maximum A Posteriori) Prediction");
        defaultComboBoxModel1.addElement("Pairing Probability Prediction");
        defaultComboBoxModel1.addElement("Stochastic Sampling");
        cmbPartsMode.setModel(defaultComboBoxModel1);
        cmbPartsMode.setName("cmbPartsMode");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlPartsOptions.add(cmbPartsMode, gbc);
        final JLabel label10 = new JLabel();
        label10.setText("Number of Stochastic Samples");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        pnlPartsOptions.add(label10, gbc);
        spnStochasticSamples = new JSpinner();
        spnStochasticSamples.setName("spnStochasticSamples");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlPartsOptions.add(spnStochasticSamples, gbc);
        final JPanel spacer3 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.VERTICAL;
        pnlPartsOptions.add(spacer3, gbc);
        final JPanel spacer4 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlPartsOptions.add(spacer4, gbc);
        btnResetPartsOptions = new JButton();
        btnResetPartsOptions.setActionCommand("reset-options-parts");
        btnResetPartsOptions.setText("Reset All to Defaults");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.EAST;
        pnlPartsOptions.add(btnResetPartsOptions, gbc);
        pnlConstraints = new JPanel();
        pnlConstraints.setLayout(new GridBagLayout());
        tabControls.addTab("Constraints", pnlConstraints);
        final JLabel label11 = new JLabel();
        label11.setText("Folding Constraints (Seq 1)");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        pnlConstraints.add(label11, gbc);
        txtConFile1 = new JTextField();
        txtConFile1.setName("txtConFile1");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 3);
        pnlConstraints.add(txtConFile1, gbc);
        final JButton button5 = new JButton();
        button5.setActionCommand("browse-con1");
        button5.setText("...");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlConstraints.add(button5, gbc);
        final JLabel label12 = new JLabel();
        label12.setText("SHAPE Slope");
        label12.setVisible(false);
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 4;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        pnlConstraints.add(label12, gbc);
        final JButton button6 = new JButton();
        button6.setActionCommand("recent-con1");
        button6.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlConstraints.add(button6, gbc);
        final JLabel label13 = new JLabel();
        label13.setText("Folding Constraints (Seq 2)");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        pnlConstraints.add(label13, gbc);
        txtConFile2 = new JTextField();
        txtConFile2.setName("txtConFile2");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 3);
        pnlConstraints.add(txtConFile2, gbc);
        final JButton button7 = new JButton();
        button7.setActionCommand("browse-con2");
        button7.setText("...");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 2;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlConstraints.add(button7, gbc);
        final JButton button8 = new JButton();
        button8.setActionCommand("recent-con2");
        button8.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 2;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlConstraints.add(button8, gbc);
        lblBaseInserts = new JLabel();
        lblBaseInserts.setText("Allow Single Base Inserts");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        pnlConstraints.add(lblBaseInserts, gbc);
        chkBaseInserts = new JCheckBox();
        chkBaseInserts.setActionCommand("");
        chkBaseInserts.setHorizontalTextPosition(2);
        chkBaseInserts.setName("chkBaseInserts");
        chkBaseInserts.setSelected(true);
        chkBaseInserts.setText("");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        pnlConstraints.add(chkBaseInserts, gbc);
        final JPanel spacer5 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.VERTICAL;
        pnlConstraints.add(spacer5, gbc);
        final JLabel label14 = new JLabel();
        label14.setText("<html>This calculation will produce multiple files in the same directory as the output file you choose below.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 4;
        gbc.gridwidth = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        contentPanel.add(label14, gbc);
        pnlOutput = new JPanel();
        pnlOutput.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 5;
        gbc.gridwidth = 2;
        gbc.fill = GridBagConstraints.BOTH;
        contentPanel.add(pnlOutput, gbc);
        final JLabel label15 = new JLabel();
        label15.setText("Output File 1");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        pnlOutput.add(label15, gbc);
        txtOutputFile1 = new JTextField();
        txtOutputFile1.setName("txtOutputFile");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.gridwidth = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlOutput.add(txtOutputFile1, gbc);
        lblOptionWarning = new JLabel();
        lblOptionWarning.setBackground(new Color(-2171392));
        lblOptionWarning.setForeground(new Color(-16777216));
        lblOptionWarning.setHorizontalAlignment(0);
        lblOptionWarning.setHorizontalTextPosition(0);
        lblOptionWarning.setOpaque(true);
        lblOptionWarning.setText("<html><center>Warning: Some options have been set to non-default values.<br>---List---");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.gridwidth = 4;
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
        gbc.gridx = 1;
        gbc.gridy = 3;
        gbc.gridwidth = 3;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlOutput.add(btnStartCalc, gbc);
        final JLabel label16 = new JLabel();
        label16.setText("Output File 2");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        pnlOutput.add(label16, gbc);
        txtOutputFile2 = new JTextField();
        txtOutputFile2.setName("txtOutputFile2");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.gridwidth = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlOutput.add(txtOutputFile2, gbc);
        btnResetAll = new JButton();
        btnResetAll.setActionCommand("reset-options-all");
        btnResetAll.setText("Reset All Options");
        gbc = new GridBagConstraints();
        gbc.gridx = 4;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.EAST;
        pnlOutput.add(btnResetAll, gbc);
        final JButton button9 = new JButton();
        button9.setActionCommand("recent-out2");
        button9.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 4;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOutput.add(button9, gbc);
        final JButton button10 = new JButton();
        button10.setActionCommand("recent-out1");
        button10.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 4;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOutput.add(button10, gbc);
        final JButton button11 = new JButton();
        button11.setActionCommand("browse-out2");
        button11.setText("...");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOutput.add(button11, gbc);
        final JButton button12 = new JButton();
        button12.setActionCommand("browse-out1");
        button12.setText("...");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOutput.add(button12, gbc);
        pnlTaskStatus = new JPanel();
        pnlTaskStatus.setLayout(new GridBagLayout());
        pnlTaskStatus.setVisible(true);
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 6;
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
        lblToolDescription.setText("<html>This tool takes <b>two sequences</b>, determines common secondary structures, and calculates base pair probabilities. <br>Additionally, it creates two <b>alignment</b> files. <br>Both sequences, as well as the resulting alignments, can be restrained using <b>constraint</b> files. <br>Note that this tool assumes the sequences are <b>RNA</b>.");
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
        label2.setLabelFor(txtSeqFile1);
        label3.setLabelFor(txtSeqFile2);
        label4.setLabelFor(spnMaxEnergyDiff);
        label5.setLabelFor(spnMaxStructs);
        label6.setLabelFor(spnWindowSize);
        label7.setLabelFor(spnGapPenalty);
        label8.setLabelFor(spnAlnWindowSize);
        label9.setLabelFor(cmbPartsMode);
        label10.setLabelFor(spnStochasticSamples);
        label11.setLabelFor(txtConFile1);
        label13.setLabelFor(txtConFile2);
        label15.setLabelFor(txtOutputFile1);
        label16.setLabelFor(txtOutputFile2);
    }
    /** @noinspection ALL */
    public JComponent $$$getRootComponent$$$() { return contentPanel; }

    private class PredictTask extends PredictTaskBase {
        public String sequence1, sequence2;
        public int seqType1, seqType2;

        //public double temperature;
        //public int maxLoopSize;
        public int maxEnergyDiff;
        public int maxStructures;
        public int windowSize;
        public int alnWindowSize;
        public double gapPenalty;
        //public int iterations;
        //public int minHelixLen;
        public boolean allowBaseInserts;
//        public double shapeSlope;
//        public double shapeIntercept;

        public String shapeFile;
        public String conFile1, conFile2;
        //public String sequenceTitle;

        // Output files
        public String outDynCt1, outDynCt2, outPartsCt1, outPartsCt2, outDsv, outDynAln;
        public int stochasticSamples;
        public int partsMode;
        public String baseOutputName1;
        public String baseOutputName2;
        public String baseOutputNameShared;

        @Override
        protected void runCalc() {
            nextStep("Loading sequences...", true, 2);

            Dynalign_object dyn = new Dynalign_object(sequence1, seqType1, sequence2, seqType2, true);
            if (checkError(dyn) || isCanceled) return;
            dyn.SetProgress(stepProgress);
            RNA rna1 = dyn.GetRNA1(), rna2 = dyn.GetRNA2();

            if (!isWhiteSpace(conFile1)) if (checkError(rna1.ReadConstraints(conFile1))) return;
            if (!isWhiteSpace(conFile2)) if (checkError(rna2.ReadConstraints(conFile2))) return;
            if (checkError(dyn) || isCanceled) return;

            outDsv = baseOutputNameShared + ".dsv";
            outDynAln = baseOutputNameShared + ".aln";
            outDynCt1 = baseOutputName1 + ".ct";
            outDynCt2 = baseOutputName2 + ".ct";

            nextStep("Running Dynalign...", true, 100);
            int ret = dyn.Dynalign((short) maxStructures, (short) windowSize, (short) alnWindowSize, (short) maxEnergyDiff, (short) -99, (float) gapPenalty, allowBaseInserts, outDsv);
            if (checkError(ret)) return;

            if (checkError(rna1.WriteCt(outDynCt1))) return;
            if (checkError(rna2.WriteCt(outDynCt2))) return;
            if (checkError(rna1, rna2)) return;

            dyn.WriteAlignment(outDynAln);
            if (checkError(dyn) || checkError(rna1, rna2)) return;
            dyn.StopProgress();
        }

        /**
         * Show the results window.
         *
         * @return True if the results window could be shown. False otherwise.
         */
        protected boolean showResults() {
            // else show results window.
            PredictionResultsWindow w = new PredictionResultsWindow("Results of Aligned Structure Prediction");
            if (isFile(outDynCt1)) w.addFile(outDynCt1, "Structures from Sequence 1");
            if (isFile(outDynCt2)) w.addFile(outDynCt2, "Structures from Sequence 2");
            if (isFile(outDsv)) w.addFile(outDsv, "Dynalign Save File");
            if (isFile(outDynAln)) w.addFile(outDynAln, "Alignment File");

//            if (isFile(outDsv)) {
//                w.addPlotHeader("Base-pairing Probabilities (from Partition Function)", 1);
//                if (isFile(outDsv))
//                    w.addPlot("Basepair Dot Plot", "Pairwise base-pairing probabilities in dot-plot form.",
//                            () -> launchDrawing(outDsv, FileType.PFS, null, null, 1));
//            }

            w.addPlotHeader("Predicted Structures", 3);
            if (isFile(outDynCt1))
                w.addPlot("Structures from Sequence 1", "Aligned structures from Sequence 1.",
                        () -> launchDrawing(outDynCt1, FileType.CT, null, null, 1));
            if (isFile(outDynCt2))
                w.addPlot("Structures from Sequence 2", "Aligned structures from Sequence 2.",
                        () -> launchDrawing(outDynCt2, FileType.CT, null, null, 1));

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
            validateFileField(txtOutputFile1, false, true);
            validateFileField(txtOutputFile2, false, true);
            validateFileField(txtSeqFile1, false);
            validateFileField(txtSeqFile2, false);
            validateFileField(txtConFile1, true);
            validateFileField(txtConFile2, true);
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
