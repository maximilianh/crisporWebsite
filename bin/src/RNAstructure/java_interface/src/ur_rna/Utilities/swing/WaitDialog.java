package ur_rna.Utilities.swing;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * A simple wait dialog
 */
public class WaitDialog extends JDialog {
    private static final String DEFAULT_TITLE = "Please wait...";
    private static final String DEFAULT_MESSAGE = "Operation in progress. Please wait...";
    private String message;
    private boolean showProgress;
    private boolean showCancel;
    private float progressValue;
    private float progressTotal = 100;
    private ActionListener cancelListener;
    private boolean canceled;

    private JLabel lblMessage;
    private JProgressBar prgProgress;
    private JButton btnCancel;
    private boolean closeOnCancel = true;
    public static WaitDialog show(final Frame owner, String message) {
        WaitDialog w = new WaitDialog(owner, message);
        w.setVisible(true);
        return w;
    }
    public WaitDialog() {
        this(null, DEFAULT_TITLE);
    }
    public WaitDialog(final Frame owner) {
        this(owner, DEFAULT_TITLE);
    }
    public WaitDialog(final String message) {
        this(null, null, message);
    }
    public WaitDialog(final Frame owner, final String title) {
        this(owner, title, DEFAULT_MESSAGE);
    }
    public WaitDialog(final Frame owner, final String title, final String message) {
        this(owner, title, message, false, false);
    }
    public WaitDialog(final Frame owner, final String title, final String message, boolean showProgress,  ActionListener cancelListener) {
        this(owner, title, message, showProgress, cancelListener!=null);
        setCancelListener(cancelListener);
    }
    public WaitDialog(final Frame owner, final String title, final String message, boolean showProgress, boolean showCancel) {
        super(owner, title==null?DEFAULT_TITLE:title);
        this.message = message==null?DEFAULT_MESSAGE:message;
        this.showProgress = showProgress;
        this.showCancel = showCancel;
        buildUI();
        updateUI();
    }
    private void updateUI() {
        lblMessage.setText(message);
        prgProgress.setVisible(showProgress);
        if (showProgress)
            prgProgress.setValue(Math.round(progressValue/progressTotal));
        btnCancel.setVisible(showCancel);
    }
    private void buildUI() {
        setResizable(false);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setSize(200, 100);
        btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(this::canceled);
        prgProgress = new JProgressBar();
        lblMessage = new JLabel();
        //setLayout(new GridBagLayout());
        //GridBagConstraints c = new GridBagConstraints();
        add(lblMessage, BorderLayout.NORTH);
        add(prgProgress, BorderLayout.CENTER);
        add(btnCancel, BorderLayout.SOUTH);
    }
    private void canceled(final ActionEvent event) {
        canceled = true;
        if (cancelListener!=null)
            cancelListener.actionPerformed(event);
        if (closeOnCancel)
            this.dispose();
    }

    public void close() { this.dispose(); }

    public String getMessage() {
        return message;
    }
    public void setMessage(final String message) {
        this.message = message;
        updateUI();
    }
    public boolean isProgressVisible() {
        return showProgress;
    }
    public WaitDialog setProgressVisible(final boolean showProgress) {
        this.showProgress = showProgress;
        updateUI();
        return this;
    }
    public boolean isCancelVisible() {
        return showCancel;
    }
    public WaitDialog setCancelVisible(final boolean showCancel) {
        this.showCancel = showCancel;
        updateUI();
        return this;
    }
    public float getProgressValue() {
        return progressValue;
    }
    public WaitDialog setProgressValue(final float progressValue) {
        this.progressValue = progressValue;
        updateUI();
        return this;
    }
    public float getProgressTotal() {
        return progressTotal;
    }
    public WaitDialog setProgressTotal(final float progressTotal) {
        this.progressTotal = progressTotal;
        updateUI();
        return this;
    }
    public ActionListener getCancelListener() {
        return cancelListener;
    }
    public WaitDialog setCancelListener(final ActionListener cancelListener) {
        this.cancelListener = cancelListener;
        return this;
    }
    public boolean isCloseOnCancel() {
        return closeOnCancel;
    }
    public WaitDialog setCloseOnCancel(final boolean closeOnCancel) {
        this.closeOnCancel = closeOnCancel;
        return this;
    }
    public boolean isCanceled() { return canceled; }
    public void cancel() {
        btnCancel.doClick();
    }
}
