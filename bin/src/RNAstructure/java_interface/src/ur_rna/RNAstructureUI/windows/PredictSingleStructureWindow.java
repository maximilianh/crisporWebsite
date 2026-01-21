/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructure.backend.RNA;
import ur_rna.RNAstructure.backend.RNABackend;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.utilities.FileFilters;
import ur_rna.RNAstructureUI.utilities.FileType;
import ur_rna.Utilities.PathTools;
import ur_rna.Utilities.Strings;
import ur_rna.Utilities.swing.FormValidationException;

import javax.swing.*;
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
public class PredictSingleStructureWindow extends PredictWindowBase {
    private static final long serialVersionUID = 20160122;

    private JPanel contentPanel;
    private JTextField txtSeqFile;
    private JRadioButton optSeqFile;
    private JButton btnBrowseSeq;
    private JTextArea txtSequence;
    private JTextField txtSeqTitle;
    private JRadioButton optRNA;
    private JSpinner spnTemperature;
    private JTabbedPane tabControls;
    private JSpinner spnMaxLoopSize;
    private JSpinner spnMaxEnergyDiff;
    private JSpinner spnMinHelixLen;
    private JSpinner spnIterations;
    private JSpinner spnGamma;
    private JSpinner spnWindowSize;
    private JSpinner spnMaxStructs;
    private JRadioButton optDNA;
    private JRadioButton optEnterSeq;
    private JTextField txtConFile;
    private JButton btnBrowseCon;
    private JTextField txtShapeFile;
    private JButton btnBrowseShape;
    private JSpinner spnShapeIntercept;
    private JSpinner spnShapeSlope;
    private JButton btnStartCalc;
    private JTextField txtOutputFile;
    private JButton btnBrowseOutput;
    private JButton btnRecentSeq;
    private JButton btnRecentOutput;
    private JButton btnRecentCon;
    private JButton btnRecentShape;
    private JButton btnResetOptions;
    private JLabel lblOptionWarning;
    private JScrollPane scrlSequence;
    private JButton btnCancelTask;
    private JProgressBar prgTaskProgress;
    private JLabel lblTaskStatus;
    private JPanel pnlTaskStatus;
    private JCheckBox chkShowDescription;
    private JLabel lblToolDescription;
    private JPanel pnlOutput;
    private JLabel lblProgressPercent;

    private HashMap<JComponent, Object> defaultOptionValues = new HashMap<>();
    private ArrayList<JComponent> optionsControls = new ArrayList<>();
    private ArrayList<JComponent> settingsControls = new ArrayList<>();
    private final String settingsNodeName = "PredictSingle";

    /**
     * Constructor.
     */
    public PredictSingleStructureWindow() {
        $$$setupUI$$$();
        autoSetFieldNames();
        setCaption("Predict a Secondary Structure");

//        listenForActions(btnBrowseCon, btnBrowseOutput, btnBrowseSeq, btnBrowseShape,
//                btnStartCalc, btnResetOptions, btnCancelTask);

        listenForActions(contentPanel, JButton.class, this);

        addFileField("seq", txtSeqFile, FileFilters.Sequence, false, "select-seqfile");
        addFileField("con", txtConFile, FileFilters.Constraints, false);
        addFileField("shape", txtShapeFile, FileFilters.SHAPE, false);
        addFileField("output", txtOutputFile, FileFilters.CT, true);

        pnlTaskStatus.setVisible(false);
        setContent(contentPanel);
        setPreferredSize(new Dimension(700, 600));

        // optionsControls is a list of all controls that we should store the default value for
        // and then check to see if the user has changed it and show a warning if they have.
        // (because options are on different tabs and the user may have forgotten they set them
        // in a previous instance of this form and their values have been restored from that previous time.)
        optionsControls.addAll(Arrays.asList(
                spnTemperature, spnMaxLoopSize, spnMaxEnergyDiff,
                spnMaxStructs, spnWindowSize, spnGamma,
                spnIterations, spnMinHelixLen,
                spnShapeIntercept, spnShapeSlope));

        // settingsControls is a list of all controls that should have their
        // values saved and restored.
        settingsControls.addAll(optionsControls);
        settingsControls.addAll(Arrays.asList(optDNA, optRNA, optEnterSeq, optSeqFile, chkShowDescription));

        loadSettingsDefaults();

        for (JComponent c : optionsControls)
            defaultOptionValues.put(c, getControlValue(c));
        handleSettings(false); // load recent settings values.

        btnResetOptions.setEnabled(false);
        lblOptionWarning.setVisible(false);

        updateFormUI();

        startUiTimer();
    }

    @Override
    public void showWindow() {
        super.showWindow();
        SwingUtilities.invokeLater(this::pack);
    }

    @SuppressWarnings("MagicNumber")
    private void loadSettingsDefaults() {
        setSpin(spnTemperature, 310.15, 0, 600, 1);
        setSpin(spnMaxLoopSize, 30, 3, 100);
        setSpin(spnMaxEnergyDiff, 10, 1, 200);
        setSpin(spnMaxStructs, 20, 1, 1000);
        setSpin(spnWindowSize, 3, 1, 100);
        setSpin(spnGamma, 1.0);
        setSpin(spnIterations, 1, 0, 1000);
        setSpin(spnMinHelixLen, 3, 1, 20);
        setSpin(spnShapeIntercept, -0.6);
        setSpin(spnShapeSlope, 1.8);
    }

    @Override
    protected void processCommand(final CommandInfo ci) {
        switch (ci.getCommand()) {
            case "select-seqfile": // called after browse-seq or recent-seq sets a new file path
                optSeqFile.setSelected(true);
                break;
            case "reset-options":
                for (JComponent c : optionsControls)
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

    private String lastSeqFile = "", lastOutFile = ""; // detect changes to these fields.

    @Override
    protected void updateFormUI() {
        // If there is a running task, just update progress and return. Do not modify any form controls.
        PredictSingleStructureTask task = (PredictSingleStructureTask) getRunningTask();
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

        boolean useSeqFile = optSeqFile.isSelected();
        // Only enable controls when there is no running tasks.
        txtSeqFile.setEnabled(useSeqFile);
        txtSequence.setEnabled(!useSeqFile);
        txtSeqTitle.setEnabled(!useSeqFile);
        scrlSequence.setEnabled(!useSeqFile);
        btnResetOptions.setEnabled(modified);

        if (useSeqFile) {
            if (!lastSeqFile.equals(txtSeqFile.getText())) {
                lastSeqFile = txtSeqFile.getText();
                // only do the following if the user has NOT changed the output file name themselves
                if (txtOutputFile.getText().isEmpty() || lastOutFile.equals(txtOutputFile.getText()))
                    txtOutputFile.setText(lastOutFile = ModuleWindow.getOutputFile(txtSeqFile.getText(), "ct"));
            }
        }

        if (lblToolDescription.isVisible() != chkShowDescription.isSelected()) {
            lblToolDescription.setVisible(chkShowDescription.isSelected());
            if (this.isVisible())
                pack();
        }
    }

    private void runCalc() {
        saveRecentFiles(txtShapeFile, txtConFile, txtSeqFile, txtOutputFile);
        handleSettings(true);
        btnStartCalc.setVisible(false);
        pnlTaskStatus.setVisible(true);
        enableDescendants(tabControls, false);
        enableDescendants(pnlOutput, false);

        setTaskStatus(0, "Starting calculation...");

        // set variables based on form input
        PredictSingleStructureTask task = new PredictSingleStructureTask();
        if (optEnterSeq.isSelected()) {
            task.sequence = txtSequence.getText();
            task.sequenceTitle = txtSeqTitle.getText();
            task.seqType = RNABackend.SEQUENCE_STRING;
        } else {
            task.sequence = txtSeqFile.getText();
            task.seqType = RNABackend.FILE_SEQ;
        }
        task.alphabet = optDNA.isSelected() ? "dna" : "rna";
        task.baseOutputName = PathTools.changeExtension(txtOutputFile.getText().trim(), "");
        task.temperature = (Double) spnTemperature.getValue();
        task.maxLoopSize = (Integer) spnMaxLoopSize.getValue();
        task.maxEnergyDiff = (Integer) spnMaxEnergyDiff.getValue();
        task.maxStructures = (Integer) spnMaxStructs.getValue();
        task.windowSize = (Integer) spnWindowSize.getValue();
        task.gamma = (Double) spnGamma.getValue();
        task.iterations = (Integer) spnIterations.getValue();
        task.minHelixLen = (Integer) spnMinHelixLen.getValue();
        task.shapeSlope = (Double) spnShapeSlope.getValue();
        task.shapeIntercept = (Double) spnShapeIntercept.getValue();
        task.shapeFile = txtShapeFile.getText();
        task.conFile = txtConFile.getText();

        runInBackground(task, this::calcComplete);
    }

    private void calcComplete(PredictSingleStructureTask task) {
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
        contentPanel.setAutoscrolls(true);
        contentPanel.setName("contentPanel");
        final JLabel label1 = new JLabel();
        label1.setFont(new Font(label1.getFont().getName(), Font.BOLD, 14));
        label1.setOpaque(false);
        label1.setText("Predict a Secondary Structure");
        GridBagConstraints gbc;
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.NORTHWEST;
        contentPanel.add(label1, gbc);
        tabControls = new JTabbedPane();
        tabControls.setName("tabbedPane1");
        tabControls.setPreferredSize(new Dimension(512, 260));
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.gridwidth = 2;
        gbc.fill = GridBagConstraints.BOTH;
        gbc.insets = new Insets(5, 0, 5, 0);
        contentPanel.add(tabControls, gbc);
        final JPanel panel1 = new JPanel();
        panel1.setLayout(new GridBagLayout());
        tabControls.addTab("Sequence", panel1);
        final JLabel label2 = new JLabel();
        label2.setText("Nucleic Acid Type");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        panel1.add(label2, gbc);
        optRNA = new JRadioButton();
        optRNA.setName("optRNA");
        optRNA.setSelected(true);
        optRNA.setText("RNA");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        panel1.add(optRNA, gbc);
        txtSeqFile = new JTextField();
        txtSeqFile.setName("txtSeqFile");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 2;
        gbc.gridwidth = 3;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        panel1.add(txtSeqFile, gbc);
        btnBrowseSeq = new JButton();
        btnBrowseSeq.setActionCommand("browse-seq");
        btnBrowseSeq.setName("btnBrowseSeq");
        btnBrowseSeq.setText("...");
        gbc = new GridBagConstraints();
        gbc.gridx = 4;
        gbc.gridy = 2;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel1.add(btnBrowseSeq, gbc);
        txtSeqTitle = new JTextField();
        txtSeqTitle.setName("txtSeqTitle");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 3;
        gbc.gridwidth = 3;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        panel1.add(txtSeqTitle, gbc);
        scrlSequence = new JScrollPane();
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 4;
        gbc.gridwidth = 6;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        gbc.insets = new Insets(3, 0, 0, 0);
        panel1.add(scrlSequence, gbc);
        txtSequence = new JTextArea();
        txtSequence.setName("txtSequence");
        txtSequence.setText("(Enter Sequence)");
        scrlSequence.setViewportView(txtSequence);
        btnRecentSeq = new JButton();
        btnRecentSeq.setActionCommand("recent-seq");
        btnRecentSeq.setName("btnRecentSeq");
        btnRecentSeq.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 5;
        gbc.gridy = 2;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel1.add(btnRecentSeq, gbc);
        final JLabel label3 = new JLabel();
        label3.setText("Sequence Source");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        panel1.add(label3, gbc);
        optSeqFile = new JRadioButton();
        optSeqFile.setName("optSeqFile");
        optSeqFile.setSelected(true);
        optSeqFile.setText("Load from File");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.gridwidth = 2;
        gbc.anchor = GridBagConstraints.WEST;
        panel1.add(optSeqFile, gbc);
        optEnterSeq = new JRadioButton();
        optEnterSeq.setName("optEnterSeq");
        optEnterSeq.setText("Enter Sequence Below");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        panel1.add(optEnterSeq, gbc);
        final JLabel label4 = new JLabel();
        label4.setText("Sequence File");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 3, 0, 3);
        panel1.add(label4, gbc);
        final JLabel label5 = new JLabel();
        label5.setText("Sequence Title");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 3, 0, 3);
        panel1.add(label5, gbc);
        optDNA = new JRadioButton();
        optDNA.setName("optDNA");
        optDNA.setText("DNA");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        panel1.add(optDNA, gbc);
        final JPanel panel2 = new JPanel();
        panel2.setLayout(new GridBagLayout());
        tabControls.addTab("Options", panel2);
        final JLabel label6 = new JLabel();
        label6.setText("Temperature");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        panel2.add(label6, gbc);
        spnTemperature = new JSpinner();
        spnTemperature.setName("spnTemperature");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel2.add(spnTemperature, gbc);
        final JLabel label7 = new JLabel();
        label7.setText("K");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        panel2.add(label7, gbc);
        final JLabel label8 = new JLabel();
        label8.setText("(all)");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        panel2.add(label8, gbc);
        final JLabel label9 = new JLabel();
        label9.setText("(used for)");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        panel2.add(label9, gbc);
        final JLabel label10 = new JLabel();
        label10.setText("Maximum Loop Size");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        panel2.add(label10, gbc);
        spnMaxLoopSize = new JSpinner();
        spnMaxLoopSize.setName("spnMaxLoopSize");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel2.add(spnMaxLoopSize, gbc);
        final JLabel label11 = new JLabel();
        label11.setText("bases");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        panel2.add(label11, gbc);
        final JLabel label12 = new JLabel();
        label12.setText("Maximum % Energy Difference");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        panel2.add(label12, gbc);
        final JLabel label13 = new JLabel();
        label13.setText("Maximum Number of Structures");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 4;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        panel2.add(label13, gbc);
        final JLabel label14 = new JLabel();
        label14.setText("Window Size");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 5;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        panel2.add(label14, gbc);
        spnMaxEnergyDiff = new JSpinner();
        spnMaxEnergyDiff.setName("spnMaxEnergyDiff");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 3;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel2.add(spnMaxEnergyDiff, gbc);
        spnMaxStructs = new JSpinner();
        spnMaxStructs.setName("spnMaxStructs");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 4;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel2.add(spnMaxStructs, gbc);
        spnWindowSize = new JSpinner();
        spnWindowSize.setName("spnWindowSize");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 5;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel2.add(spnWindowSize, gbc);
        final JLabel label15 = new JLabel();
        label15.setText("(all)");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        panel2.add(label15, gbc);
        final JLabel label16 = new JLabel();
        label16.setText("(MFE, MEA)");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 3;
        gbc.anchor = GridBagConstraints.WEST;
        panel2.add(label16, gbc);
        final JLabel label17 = new JLabel();
        label17.setText("(MFE, MEA)");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 4;
        gbc.anchor = GridBagConstraints.WEST;
        panel2.add(label17, gbc);
        final JLabel label18 = new JLabel();
        label18.setText("(MFE, MEA)");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 5;
        gbc.anchor = GridBagConstraints.WEST;
        panel2.add(label18, gbc);
        spnGamma = new JSpinner();
        spnGamma.setName("spnGamma");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 6;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel2.add(spnGamma, gbc);
        final JLabel label19 = new JLabel();
        label19.setText("Gamma");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 6;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        panel2.add(label19, gbc);
        final JLabel label20 = new JLabel();
        label20.setText("(MEA)");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 6;
        gbc.anchor = GridBagConstraints.WEST;
        panel2.add(label20, gbc);
        final JLabel label21 = new JLabel();
        label21.setText("Iterations");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 7;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        panel2.add(label21, gbc);
        final JLabel label22 = new JLabel();
        label22.setText("Minimum Helix Length");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 8;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 0, 0, 2);
        panel2.add(label22, gbc);
        spnIterations = new JSpinner();
        spnIterations.setName("spnIterations");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 7;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel2.add(spnIterations, gbc);
        spnMinHelixLen = new JSpinner();
        spnMinHelixLen.setName("spnMinHelixLen");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 8;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel2.add(spnMinHelixLen, gbc);
        final JLabel label23 = new JLabel();
        label23.setText("(Pseudoknot Prediction)");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 7;
        gbc.anchor = GridBagConstraints.WEST;
        panel2.add(label23, gbc);
        final JLabel label24 = new JLabel();
        label24.setText("(Pseudoknot Prediction)");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 8;
        gbc.anchor = GridBagConstraints.WEST;
        panel2.add(label24, gbc);
        final JLabel label25 = new JLabel();
        label25.setText("bases");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 8;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        panel2.add(label25, gbc);
        final JPanel spacer1 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 4;
        gbc.gridy = 8;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        panel2.add(spacer1, gbc);
        btnResetOptions = new JButton();
        btnResetOptions.setActionCommand("reset-options");
        btnResetOptions.setText("Reset All to Defaults");
        gbc = new GridBagConstraints();
        gbc.gridx = 4;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        panel2.add(btnResetOptions, gbc);
        final JPanel panel3 = new JPanel();
        panel3.setLayout(new GridBagLayout());
        tabControls.addTab("Constraints", panel3);
        final JLabel label26 = new JLabel();
        label26.setText("Folding Constraints File");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        panel3.add(label26, gbc);
        txtConFile = new JTextField();
        txtConFile.setName("txtConFile");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.gridwidth = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 3);
        panel3.add(txtConFile, gbc);
        btnBrowseCon = new JButton();
        btnBrowseCon.setActionCommand("browse-con");
        btnBrowseCon.setName("btnBrowseCon");
        btnBrowseCon.setText("...");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel3.add(btnBrowseCon, gbc);
        final JLabel label27 = new JLabel();
        label27.setText("SHAPE File");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        panel3.add(label27, gbc);
        txtShapeFile = new JTextField();
        txtShapeFile.setName("txtShapeFile");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.gridwidth = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 3);
        panel3.add(txtShapeFile, gbc);
        btnBrowseShape = new JButton();
        btnBrowseShape.setActionCommand("browse-shape");
        btnBrowseShape.setName("btnBrowseShape");
        btnBrowseShape.setText("...");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel3.add(btnBrowseShape, gbc);
        final JLabel label28 = new JLabel();
        label28.setText("SHAPE Intercept");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        panel3.add(label28, gbc);
        final JLabel label29 = new JLabel();
        label29.setText("SHAPE Slope");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        panel3.add(label29, gbc);
        spnShapeIntercept = new JSpinner();
        spnShapeIntercept.setName("spnShapeIntercept");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel3.add(spnShapeIntercept, gbc);
        spnShapeSlope = new JSpinner();
        spnShapeSlope.setName("spnShapeSlope");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 3;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel3.add(spnShapeSlope, gbc);
        final JPanel spacer2 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 2;
        gbc.weightx = 2.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        panel3.add(spacer2, gbc);
        final JPanel spacer3 = new JPanel();
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 4;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.VERTICAL;
        panel3.add(spacer3, gbc);
        btnRecentCon = new JButton();
        btnRecentCon.setActionCommand("recent-con");
        btnRecentCon.setName("btnRecentCon");
        btnRecentCon.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 4;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel3.add(btnRecentCon, gbc);
        btnRecentShape = new JButton();
        btnRecentShape.setActionCommand("recent-shape");
        btnRecentShape.setName("btnRecentShape");
        btnRecentShape.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 4;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        panel3.add(btnRecentShape, gbc);
        final JLabel label30 = new JLabel();
        label30.setText("<html>This calculation will produce multiple files in the same directory as the output file you choose below.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.gridwidth = 2;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        contentPanel.add(label30, gbc);
        pnlOutput = new JPanel();
        pnlOutput.setLayout(new GridBagLayout());
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 4;
        gbc.gridwidth = 2;
        gbc.fill = GridBagConstraints.BOTH;
        contentPanel.add(pnlOutput, gbc);
        final JLabel label31 = new JLabel();
        label31.setText("Output File");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(0, 2, 0, 3);
        pnlOutput.add(label31, gbc);
        txtOutputFile = new JTextField();
        txtOutputFile.setName("txtOutputFile");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.weightx = 1.0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlOutput.add(txtOutputFile, gbc);
        btnBrowseOutput = new JButton();
        btnBrowseOutput.setActionCommand("browse-output");
        btnBrowseOutput.setName("btnBrowseOutput");
        btnBrowseOutput.setText("...");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOutput.add(btnBrowseOutput, gbc);
        btnRecentOutput = new JButton();
        btnRecentOutput.setActionCommand("recent-output");
        btnRecentOutput.setName("btnRecentOutput");
        btnRecentOutput.setText("Recent");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 2, 0, 2);
        pnlOutput.add(btnRecentOutput, gbc);
        lblOptionWarning = new JLabel();
        lblOptionWarning.setBackground(new Color(-2171392));
        lblOptionWarning.setForeground(new Color(-16777216));
        lblOptionWarning.setHorizontalAlignment(0);
        lblOptionWarning.setHorizontalTextPosition(0);
        lblOptionWarning.setOpaque(true);
        lblOptionWarning.setText("<html><center>Warning: Some options have been set to non-default values.<br>---List---");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
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
        gbc.gridy = 2;
        gbc.gridwidth = 2;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlOutput.add(btnStartCalc, gbc);
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
        lblToolDescription.setText("<html>This tool combines four separate prediction and analysis routines: <ol><li>Calculate a partition function (<b>PF</b>).</li><li>Predict a maximum free energy (<b>MFE</b>) structure.</li><li>Find structures with maximum expected accuracy (<b>MEA</b>).</li><li>Predict pseudoknots (<b>PK</b>).</li></ol>This generates a highly probable, probability-annotated list of secondary structures, starting with the lowest free energy structure and including others with high probabilities of correctness. SHAPE constraints may be specified and are applied to the probability-annotated structures.<br>If SHAPE contraints are specified, a second group of SHAPE constrained, SHAPE annotated structures are be generated. This SHAPE structure group is distinct from the probability-annotated structure group, and is not probability-annotated itself. </html>");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.gridwidth = 2;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        contentPanel.add(lblToolDescription, gbc);
        label4.setLabelFor(txtSeqFile);
        label5.setLabelFor(txtSeqTitle);
        label6.setLabelFor(spnTemperature);
        label10.setLabelFor(spnMaxLoopSize);
        label12.setLabelFor(spnMaxEnergyDiff);
        label13.setLabelFor(spnMaxStructs);
        label14.setLabelFor(spnWindowSize);
        label19.setLabelFor(spnGamma);
        label21.setLabelFor(spnIterations);
        label22.setLabelFor(spnMinHelixLen);
        label26.setLabelFor(txtConFile);
        label27.setLabelFor(txtShapeFile);
        label28.setLabelFor(spnShapeIntercept);
        label29.setLabelFor(spnShapeSlope);
        label31.setLabelFor(txtOutputFile);
        ButtonGroup buttonGroup;
        buttonGroup = new ButtonGroup();
        buttonGroup.add(optRNA);
        buttonGroup.add(optDNA);
        buttonGroup = new ButtonGroup();
        buttonGroup.add(optSeqFile);
        buttonGroup.add(optEnterSeq);
    }
    /** @noinspection ALL */
    public JComponent $$$getRootComponent$$$() { return contentPanel; }

    private class PredictSingleStructureTask extends PredictTaskBase {
        public String sequence;
        public int seqType;
        public String alphabet;
        public String baseOutputName;

        public double temperature;
        public int maxLoopSize;
        public int maxEnergyDiff;
        public int maxStructures;
        public int windowSize;
        public double gamma;
        public int iterations;
        public int minHelixLen;
        public double shapeSlope;
        public double shapeIntercept;

        public String shapeFile;
        public String conFile;
        public String sequenceTitle;

        // Output files
        public String outFoldCt, outMaxExpectCt, outPfs, outPkCt, outSav;

        @Override
        protected void runCalc() {
            nextStep("Loading sequence...", true, 2);
            RNA rna = new RNA(sequence, seqType, alphabet, false, false);
            if (!isEmpty(sequenceTitle)) rna.GetStructure().SetSequenceLabel(sequenceTitle);
            if (temperature != TEMP_37C) rna.SetTemperature(temperature);
            if (checkError(rna) || isCanceled) return;
            if (!Strings.isWhiteSpace(conFile)) rna.ReadConstraints(conFile);
            if (checkError(rna) || isCanceled) return;
            if (!Strings.isWhiteSpace(shapeFile)) rna.ReadSHAPE(shapeFile, shapeSlope, shapeIntercept);
            if (checkError(rna) || isCanceled) return;

            rna.SetProgress(stepProgress);

            nextStep("Running partition...", true, 50);
            rna.PartitionFunction(outPfs = baseOutputName + ".pfs", temperature);
            if (checkError(rna) || isCanceled) return;

            nextStep("Running MaxExpect...", true, 5);
            rna.MaximizeExpectedAccuracy(maxEnergyDiff, maxStructures, windowSize, gamma);
            rna.WriteCt(outMaxExpectCt = baseOutputName + "_MaxExpect.ct");
            if (checkError(rna) || isCanceled) return;
            rna.GetStructure().RemoveAllStructures(); //remove structures already written to file so they don't accumulate.

            nextStep("Running ProbKnot...", false, 3);
            rna.ProbKnot(iterations, minHelixLen);
            rna.WriteCt(outPkCt = baseOutputName + "_ProbKnot.ct");
            if (checkError(rna) || isCanceled) return;
            rna.GetStructure().RemoveAllStructures(); //remove structures already written to file so they don't accumulate.

            nextStep("Running Fold...", true, 40);
            rna.FoldSingleStrand(maxEnergyDiff, maxStructures, windowSize,
                    outSav = baseOutputName + ".sav", maxLoopSize);
            rna.WriteCt(outFoldCt = baseOutputName + "_Fold.ct");
            if (checkError(rna) || isCanceled) return;
            rna.GetStructure().RemoveAllStructures(); //remove structures already written to file so they don't accumulate.

            nextStep("Done.", false, 0);
        }
        /**
         * Show the results window.
         *
         * @return True if the results window could be shown. False otherwise.
         */
        protected boolean showResults() {
            // else show results window.
            PredictionResultsWindow w = new PredictionResultsWindow("Results of Single Structure Prediction");
            if (isFile(outFoldCt)) w.addFile(outFoldCt, "MFE (Fold) Structure Prediction");
            if (isFile(outMaxExpectCt)) w.addFile(outMaxExpectCt, "MEA (MaxExpect) Structure Prediction");
            if (isFile(outPkCt)) w.addFile(outPkCt, "ProbKnot Pseudoknot Structure Prediction");
            if (isFile(outPfs)) w.addFile(outPfs, "Partition Function File (Base-pair Probabilities)");
            if (isFile(outPfs)) {
                w.addPlotHeader("Base-pairing Probabilities (from Partition Function)", 1);
                if (isFile(outPfs))
                    w.addPlot("Basepair Dot Plot", "Pairwise base-pairing probabilities in dot-plot form.",
                            () -> launchDrawing(outPfs, FileType.PFS, null, null, 1));
            }

            w.addPlotHeader("Predicted Structures (with no color annotations).", 3);
            if (isFile(outFoldCt))
                w.addPlot("MFE Structure", "Minimum Free Energy structure drawing (Without color annotation).",
                        () -> launchDrawing(outFoldCt, FileType.CT, null, null, 1));
            if (isFile(outMaxExpectCt))
                w.addPlot("MEA Structure", "Maximum Expected Accuracy structure drawing (Without color annotation).",
                        () -> launchDrawing(outMaxExpectCt, FileType.CT, null, null, 1));
            if (isFile(outPkCt))
                w.addPlot("Pseudoknot Structure", "Structure drawing with possible pseudoknots (Without color annotation).",
                        () -> launchDrawing(outPkCt, FileType.CT, null, null, 1));

            w.addPlotHeader("Predicted Structures color-annotated with base-pairing probabilities.", 3);
            if (isFile(outFoldCt))
                w.addPlot("MFE Structure + Probabilities", "Minimum Free Energy structure drawing color-annotated with base-pairing probabilities.",
                        () -> launchDrawing(outFoldCt, FileType.CT, outPfs, FileType.PFS, 1));
            if (isFile(outMaxExpectCt))
                w.addPlot("MEA Structure + Probabilities", "Maximum Expected Accuracy structure drawing color-annotated with base-pairing probabilities.",
                        () -> launchDrawing(outMaxExpectCt, FileType.CT, outPfs, FileType.PFS, 1));
            if (isFile(outPkCt))
                w.addPlot("Pseudoknot Structure + Probabilities", "Structure drawing with possible pseudoknots color-annotated with base-pairing probabilities.",
                        () -> launchDrawing(outPkCt, FileType.CT, outPfs, FileType.PFS, 1));

            if (isFile(shapeFile)) {
                w.addPlotHeader("Predicted Structures color-annotated by SHAPE reactivity.", 3);
                if (isFile(outFoldCt))
                    w.addPlot("MFE Structure + SHAPE", "Minimum Free Energy structure drawing color-annotated by SHAPE reactivity.",
                            () -> launchDrawing(outFoldCt, FileType.CT, shapeFile, FileType.SHAPE, 1));
                if (isFile(outMaxExpectCt))
                    w.addPlot("MEA Structure + SHAPE", "Maximum Expected Accuracy structure drawing color-annotated by SHAPE reactivity.",
                            () -> launchDrawing(outMaxExpectCt, FileType.CT, shapeFile, FileType.SHAPE, 1));
                if (isFile(outPkCt))
                    w.addPlot("Pseudoknot Structure + SHAPE", "Structure drawing with possible pseudoknots color-annotated by SHAPE reactivity.",
                            () -> launchDrawing(outPkCt, FileType.CT, shapeFile, FileType.SHAPE, 1));
            }
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
            validateFileField(txtOutputFile, false, true);
            if (optSeqFile.isSelected()) {
                validateFileField(txtSeqFile, false);
            } else {
                validateRequiredField(txtSeqTitle);
                validateRequiredField(txtSequence);
            }
            validateFileField(txtConFile, true);
            validateFileField(txtShapeFile, true);
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
