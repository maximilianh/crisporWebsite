package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructure.backend.ProgressHandler;
import ur_rna.RNAstructure.backend.RNA;
import ur_rna.RNAstructure.backend.SimpleProgressHandler;
import ur_rna.RNAstructure.backend.TwoRNA;
import ur_rna.RNAstructureUI.RNAstructure;
import ur_rna.RNAstructureUI.RNAstructureBackendCalculator;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.ui.RecentFileButton;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.FileType;
import ur_rna.RNAstructureUI.utilities.MRUFileStorage;
import ur_rna.Utilities.Convert;
import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.swing.FormValidationException;

import javax.swing.*;
import javax.swing.event.InternalFrameEvent;
import javax.swing.text.JTextComponent;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.lang.reflect.Field;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.function.Consumer;
import java.util.prefs.Preferences;

import static ur_rna.Utilities.PathTools.isFile;
import static ur_rna.Utilities.Strings.isEmpty;
import static ur_rna.Utilities.Strings.isWhiteSpace;

/**
 * A base class for the Structure Prediction Windows
 * @author Richard M. Watson
 */
public abstract class PredictWindowBase  extends InternalWindow {
    public static final double TEMP_37C  = 310.15;

    private Timer uiUpdateTimer;
    protected JPanel pnlClient = new JPanel();
    protected JPanel pnlContent = new JPanel();

    class FileFieldInfo {
        boolean isOutput;
        String filters;
        JTextComponent txt;
        String actionWhenSet;
        FileFieldInfo(final JTextComponent txt, final String filters, final boolean output, final String actionWhenSet) {
            isOutput = output;
            this.txt=txt;
            this.filters=filters;
            this.actionWhenSet=actionWhenSet;
        }
    }
    class FieldInfo {
        public String name;
        public String desc;
        public JComponent field;
        public FieldInfo(final JComponent field, final String desc, final String name) {
            this.name = name;
            this.desc = desc;
            this.field = field;
        }
    }
    private HashMap<String,FileFieldInfo> fileFields = new HashMap<>();
    private HashMap<JComponent,FieldInfo> fields = new HashMap<>();
    protected void addFileField(String actionLookupName, JTextComponent txt, String filters) { addFileField(actionLookupName, txt, filters, false, null); }
    protected void addFileField(String actionLookupName, JTextComponent txt, String filters, boolean isOutput) { addFileField(actionLookupName, txt, filters, isOutput, null); }
    protected void addFileField(String actionLookupName, JTextComponent txt, String filters, boolean isOutput, String action) {
        fileFields.put(actionLookupName, new FileFieldInfo(txt, filters, isOutput, action));
    }
    protected void addField(JComponent c, String description, String settingsName) {
        fields.put(c, new FieldInfo(c, description, settingsName));
    }

    protected FileFieldInfo getFileField(String fieldName, boolean showErrorIfMissing) {
        FileFieldInfo ff = fileFields.get(fieldName);
        if (ff==null&&showErrorIfMissing)
            Dialogs.showWarning("Programming error: No file field info for \"" + fieldName + "\"");
        return ff;
    }
    protected FileFieldInfo getFileField(JTextComponent txt, boolean showErrorIfMissing) {
        for(FileFieldInfo ff : fileFields.values())
            if (txt==ff.txt)
                return ff;
        if (showErrorIfMissing)
            Dialogs.showWarning("Programming error: No file field info for component \"" + getControlDescription(txt) + "\"");
        return null;
    }


    public PredictWindowBase() {
        setLayout(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.BOTH;
        gbc.weightx = gbc.weighty = 1;
        add(pnlClient, gbc);
        pnlClient.setLayout(new GridBagLayout());
    }
    void setContent(JPanel content) {
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets.set(5, 5, 5, 5);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.weightx = 1;
        pnlClient.add(content, gbc);
        pnlContent = content;
    }

    protected void startUiTimer() {
        if (uiUpdateTimer==null)
            uiUpdateTimer = new Timer(200, e -> updateFormUI());
        uiUpdateTimer.start();
    }
    protected void stopUiTimer() { if (uiUpdateTimer!=null) uiUpdateTimer.stop(); }

    /**
     * Called periodically by a timer to update the user interface.
     */
    protected abstract void updateFormUI();

    protected void listenForActions(JComponent parent) { listenForActions(parent, AbstractButton.class, this);}
    protected void listenForActions(JComponent parent, Class<? extends AbstractButton> type, ActionListener listener) {
        if (listener == null) listener = this;
        for (Component c : parent.getComponents()) {
            if (type.isAssignableFrom(c.getClass()))
                ((AbstractButton)c).addActionListener(listener);
            if (c instanceof JComponent && ((JComponent)c).getComponentCount()!=0)
                listenForActions((JComponent)c, type, listener);
        }
    }

    protected void listenForActions(AbstractButton... buttons) {
        for (AbstractButton b : buttons) b.addActionListener(this);
    }

    private PredictSingleStructureWindow.BackgroundTask runningTask = null;
    protected BackgroundTask getRunningTask() { return runningTask; }

    @Override
    protected boolean frameClosing(final InternalFrameEvent e) {
        handleSettings(true);
        return cancelRunningTask() && super.frameClosing(e);
    }

    @Override
    public void actionPerformed(final ActionEvent e) {
        String action = e.getActionCommand();
        Object source = e.getSource();
        if (action.startsWith("browse-"))
            browseFile(getFileField(action.substring(7), true));
        else if (action.startsWith("recent-"))
            browseRecentFile(source instanceof JComponent ? (JComponent)source : null, getFileField(action.substring(7), true));
        else
            super.actionPerformed(e);
    }
    @Override
    public void pack() {
        super.pack();
        if (pnlContent != null)
            setSize(getWidth(), pnlContent.getHeight() + 10 + getHeight() - pnlClient.getHeight());
    }
    /**
     * Cancel the active background task (if there is one) if it is cancelable.
     * Show a message to the user if the task cannot be canceled at the current time.
     * @return Returns true if there is no task or if it is cancelable.
     *         Returns false if the active task cannot be canceled at the current time.
     */
    protected boolean cancelRunningTask() {
        if (runningTask==null) return true;
        if (runningTask.isCancelable()) {
            runningTask.cancel();
            return true;
        }
        Dialogs.showWarning("The current running calculation cannot be canceled at this time.\nPlease wait for it to complete.\n"+
                "This calculation may have periods during which it can be canceled and others when it cannot.");
        return false;
    }

    @Override
    public void showWindow() {
        //already packed
        setLocation(0, 0);
        pack();
        setVisible(true);
    }

    protected DecimalFormat numberFormatter = new DecimalFormat("0.#######");
    protected String formatValue(Object value) {
        if (value instanceof Number)
            return numberFormatter.format(value);
        return value.toString();
    }

    protected String getControlDescription(final JComponent c) {
        FieldInfo fi = fields.get(c);
        if (fi!=null&&fi.desc!=null) return fi.desc;

        if (c instanceof AbstractButton) {
            String text = ((AbstractButton) c).getText();
            if (!isWhiteSpace(text)) return text;
        }
        for (Component other : c.getParent().getComponents())
            if (other instanceof JLabel && ((JLabel) other).getLabelFor() == c) {
                String text = ((JLabel) other).getText().trim();
                // remove colon at the end.
                if (text.endsWith(":"))
                    text = text.substring(0, text.length()-1);
                return text;
            }
        return stripControlPrefix(getSettingsName(c));
    }
    private static String[] controlPrefixes = "txt lbl chk opt rdo spn pnl btn cmd scrl fld cmp cmb lst cbo uic ui".split(" ");
    protected String stripControlPrefix(String name) {
        for(String pfx : controlPrefixes)
            if (name.startsWith(pfx))
                return name.substring(pfx.length());
        return name;
    }

    protected interface BackgroundTask extends Runnable {
        int getProgress();
        String getStatus();
        void run();
        boolean isCancelable();
        void cancel();
    }

    protected abstract void handleSettings(boolean writeValues);

    protected Object getControlValue(final JComponent c) {
        if (c instanceof JSpinner)
            return ((JSpinner) c).getValue();
        else if (c instanceof AbstractButton)
            return ((AbstractButton) c).isSelected();
        else if (c instanceof JTextField)
            return ((JTextField) c).getText();
        else if (c instanceof JComboBox)
            return ((JComboBox) c).getSelectedIndex();
        else
            Dialogs.showWarning("Unknown field type in getControlValue: " + c.getClass().getTypeName());
        return null;
    }
    protected void setControlValue(final JComponent c, Object value) {
        if (c instanceof JSpinner)
            ((JSpinner) c).setValue(Convert.toType(value, ((JSpinner) c).getValue()));
        else if (c instanceof AbstractButton)
            ((AbstractButton) c).setSelected(Convert.toBool(value));
        else if (c instanceof JTextField)
            ((JTextField) c).setText(value.toString());
        else if (c instanceof JComboBox)
            ((JComboBox) c).setSelectedIndex(Convert.toInt(value));
        else
            Dialogs.showWarning("Unknown field type in setControlValue: " + c.getClass().getTypeName());
    }

    protected void handleRecentControlValues(boolean writeValues, String nodeName, Iterable<JComponent> fields) {
        Preferences p = RNAstructure.getPrefs().node(nodeName);
        for (JComponent c : fields) {
            String key = getSettingsName(c);
            Object value = getControlValue(c);
            if (writeValues) {
                if (value == null)
                    p.remove(key);
                else
                    p.put(key, formatValue(value));
            } else {
                String defval = value == null ? null : value.toString();
                setControlValue(c, p.get(key, defval));
            }
        }
    }
    protected void loadRecent(String nodeName, JComponent... fields) {
        handleRecentControlValues(false, nodeName, Arrays.asList(fields));
    }
    protected void saveRecent(String nodeName, JComponent... fields) {
        handleRecentControlValues(true, nodeName, Arrays.asList(fields));
    }

    /**
     * Automatically set the names of controls (JComponents) that are
     * direct Fields of this class.
     */
    protected void autoSetFieldNames() {
        Field[] fields = this.getClass().getDeclaredFields();
        for (Field f : fields)
            try {
                f.setAccessible(true);
                if (JComponent.class.isAssignableFrom(f.getType())) {
                    JComponent c = (JComponent) f.get(this);
                    if (c != null && isEmpty(c.getName()))
                        c.setName(f.getName());
                    //For debugging:  if (c != null && !isEmpty(c.getName())) c.setToolTipText(c.getName());
                }
            } catch (Exception ex) {
                System.out.println("Error in autoSetFieldNames: " + ex.toString());
            }
    }

    protected String getSettingsName(JComponent c) {
        FieldInfo fi = fields.get(c);
        if (fi!=null&&fi.name!=null) return fi.name;

        String name = c.getName();
        if (name != null) return name;
        if (c instanceof JRadioButton) return "opt" + ((JRadioButton) c).getText();
        if (c instanceof JCheckBox) return "chk" + ((JCheckBox) c).getText();
        Dialogs.showWarning("Unknown field type in getSettingsName: " + c.getClass().getTypeName());
        return c.getClass().getTypeName();
    }

    protected void saveRecentFiles(JTextField... fields) {
        MRUFileStorage recent = RNAstructure.MRUFiles;
        for (JTextField f : fields) {
            String file = f.getText();
            if (!isEmpty(file))
                recent.add(file);
        }
        recent.saveToStorage();
    }

    protected void browseRecentFile(JComponent displayTarget, FileFieldInfo field) {
        if (field != null)
            RecentFileButton.showRecentFileMenu(displayTarget, field.filters, f -> {
                field.txt.setText(f.toString());
                if (!isEmpty(field.actionWhenSet))
                    actionPerformed(new ActionEvent(field.txt, ActionEvent.ACTION_PERFORMED, field.actionWhenSet));
            });
    }
    protected boolean browseFile(FileFieldInfo field) {
        if (field!=null && browseFile(field.txt, field.filters, field.isOutput)) {
            if (!isEmpty(field.actionWhenSet))
                actionPerformed(new ActionEvent(field.txt, ActionEvent.ACTION_PERFORMED, field.actionWhenSet));
            return true;
        }
        return false;
    }
    protected boolean browseFile(JTextComponent field, String filter, boolean forSave) {
        String file;
        if (forSave)
            file = StandardFileChooser.getSaveName(filter, field.getText());
        else
            file = StandardFileChooser.getOpenName(filter, field.getText());
        if (file != null)
            field.setText(file);
        return file != null;
    }

    protected <TTask extends BackgroundTask> void runInBackground(final TTask task, final Consumer<TTask> done) {
        if (runningTask!=null) {
            Dialogs.showWarning("Another calculation is currently running. Please wait for that one to finish before starting a new one.");
            return;
        }
        SwingWorker<Void,Void>  worker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                runningTask = task;
                task.run();
                return null;
            }
            @Override
            protected void done() {
                done.accept(task);
                runningTask = null;
            }
        };
        worker.execute();
    }

    protected void enableDescendants(Container parent, boolean enabled, Component...exclude) {
        for(Component c : parent.getComponents()) {
            if (ObjTools.contains(exclude, c)) continue;
            c.setEnabled(enabled);
            if (c instanceof Container && ((Container) c).getComponentCount()!=0)
                enableDescendants(((Container) c), enabled, exclude);
        }
    }

    protected static DrawingWindow launchDrawing(String file, FileType type, String addFile, FileType addType, int strand) {
        final String MISSING_FILE_MESSAGE = "The drawing cannot be created because a necessary file was not found.\n" +
                "(Possibly because it was moved or deleted after the calculation.)\nFile: ";
        if (!isFile(file)) {
            Dialogs.showError(MISSING_FILE_MESSAGE + file);
            return null;
        }
        if (addFile != null && !isFile(addFile)) {
            Dialogs.showError(MISSING_FILE_MESSAGE + addFile);
            return null;
        }
        DrawingWindow dw = new DrawingWindow(file, strand, type);
        if (dw.isError()) {
            Dialogs.showWarning("Failed to draw: " + file);
            dw.dispose();
            return null;
        } else {
            dw.showWindow();
            if (addFile != null) {
                dw.getStructure().setAnnotation(addFile, addType);
                if (dw.isError())
                    Dialogs.showWarning("Failed to add color annotations from file: " + addFile);
            }
        }
        return dw;
    }

    /**
     * Initialize a spinner with a SpinnerNumberModel.
     * @param spn The spinner
     * @param value The initial value.
     * @param min The minimum value.
     * @param max The maximim value.
     * @param step Teh step size.
     */
    protected void setSpin(JSpinner spn, double value, double min, double max, double step) {        spn.setModel(new SpinnerNumberModel(value, min, max, step));    }
    protected void setSpin(JSpinner spn, int value, int min, int max, int step) {        spn.setModel(new SpinnerNumberModel(value, min, max, step));    }
    protected void setSpin(JSpinner spn, double value, double min, double max) {        spn.setModel(new SpinnerNumberModel(value, min, max, 0.1));    }
    protected void setSpin(JSpinner spn, int value, int min, int max) {        spn.setModel(new SpinnerNumberModel(value, min, max, 1));    }
    protected void setSpin(JSpinner spn, double value) {        spn.setModel(new SpinnerNumberModel(value, null, null, 1));    }
    protected void setSpin(JSpinner spn, int value) {        spn.setModel(new SpinnerNumberModel(value, null, null, 1));    }
    protected int getSpinInt(JSpinner spn) { return  ((Number)spn.getValue()).intValue(); }
    protected double getSpinDbl(JSpinner spn) { return  ((Number)spn.getValue()).doubleValue(); }
    protected Number getSpin(JSpinner spn) { return  (Number)spn.getValue(); }

    public static abstract class PredictTaskBase implements BackgroundTask {
        public boolean hadError;
        public String errorMessage;
        protected ProgressHandler stepProgress = new SimpleProgressHandler();
        protected String status = "Starting calculation...";
        protected boolean canCancel = false;
        protected boolean isCanceled = false;

        protected void setError(String message) {
            hadError = true;
            errorMessage = message;
        }

        protected int workDone, nextStepWork;
        protected void nextStep(String stepStatus, boolean stepCancel, int stepWork) {
            canCancel = stepCancel;
            status = stepStatus;
            workDone += nextStepWork;
            stepProgress.update(0);
            nextStepWork = stepWork;
        }
        @Override
        public int getProgress() {
            return workDone + stepProgress.progress() * nextStepWork / 100;
        }
        @Override
        public String getStatus() {
            return status;
        }
        @Override
        public boolean isCancelable() {
            return canCancel;
        }
        @Override
        public void cancel() {
            //overallProgress.cancel();
            stepProgress.cancel();
            isCanceled = true;
            status = "Canceling calculation...";
        }

        @Override
        public void run() {
            try {
                runCalc();
                if (!hadError) nextStep("Done.", false, 0);
            } catch (Exception ex) {
                if (!hadError) setError("Uncaught exception in runCalc: " + ex.toString());
            }
        }

        /**
         * Show the results window.
         *
         * @return True if the results window could be shown. False otherwise.
         */
        public boolean onComplete() {
            try {
                if (isCanceled) return false;
                Thread.yield(); // allow the UI thread some time to refresh
                if (hadError) {
                    Dialogs.showWarning(errorMessage);
                    return false;
                }
                return showResults();
            } catch (Exception ex) {
                Dialogs.showError("An unexpected error occurred while post-processing the calculations results:\n" + ex.toString());
            }
            return false;
        }

        protected abstract void runCalc();
        protected abstract boolean showResults();

        protected boolean checkError(String... list) {
            for (String s : list)
                if (!isWhiteSpace(s)) {
                    setError(s);
                    return true;
                }
            return false;
        }
        protected boolean checkError(RNAstructureBackendCalculator... list) {
            for (RNAstructureBackendCalculator rna : list)
                if (rna.GetErrorCode() != 0) {
                    setError(rna.GetFullErrorMessage());
                    return true;
                }
            return false;
        }
        protected boolean checkError(RNA... rnaList) {
            for (RNA rna : rnaList)
                if (rna.GetErrorCode() != 0) {
                    setError(rna.GetFullErrorMessage());
                    return true;
                }
            return false;
        }
        protected boolean checkError(int... list) {
            for (int code : list)
                if (code != 0) {
                    setError(RNA.GetErrorMessage(code));
                    return true;
                }
            return false;
        }
        protected boolean checkError(TwoRNA... rnaList) {
            for (TwoRNA rna : rnaList)
                if (rna.GetErrorCode() != 0) {
                    setError(rna.GetErrorMessage(rna.GetErrorCode())+rna.GetErrorDetails());
                    return true;
                }
            return false;
        }
    }


    protected void validateFailed(JComponent control, String message, Object...args) throws FormValidationException {
        if (args.length!=0) message = String.format(message, args);
        if (control!=null)
            message = message.replace("<FIELD>", getControlDescription(control));
        throw new FormValidationException(message, control);
    }

    /**
     * Verifies the existence of a file in a JTextComponent that has been added as a named file-field.
     * @param txt The JTextComponent used to accept file entry.
     * @param allowEmpty Whether the field is optional (and so can be empty and still be valid input).
     * @param forOutput Whether the field is the name of an output file. (If so, it does not need to exist, but must be in a valid directory.
     * @throws FormValidationException if txt fails validation.
     */
    protected void validateFileField(JTextComponent txt, boolean allowEmpty, boolean forOutput)
            throws FormValidationException {
        String file = txt.getText().trim();
        if (file.isEmpty()) {
            if (!allowEmpty)
                validateFailed(txt,"The <FIELD> entry cannot be blank. Please enter a valid file path.");
        } else if (forOutput) {
            File dir = new File(file).getParentFile();
            if (!dir.isDirectory())
                validateFailed(txt,"The <FIELD> must be saved in an existing directory.\n(%s does not exist)", dir.toString());
        } else {
            if (!isFile(file))
                validateFailed(txt, "The <FIELD> was not found. Please enter a valid file path.");
        }
    }
    protected void validateFileField(JTextComponent txt, boolean allowEmpty) throws FormValidationException { validateFileField(txt, allowEmpty, false); }
    protected void validateRequiredField(JTextComponent txt)
            throws FormValidationException {
        String text = txt.getText().trim();
        if (text.isEmpty())
            validateFailed(txt,"The <FIELD> entry cannot be blank.");
    }
    /**
     * Attempts to show the field and set focus to it, possibly by activating the tab on which the control is placed.
     * @param c The JComponent to show and set focus to.
     */
    protected void focusField(JComponent c) {
        if (c == null) return;
        //Container parent = c.getParent();
        //while(parent != null)
        // Just try the simple method first.
        c.grabFocus();
    }
}
