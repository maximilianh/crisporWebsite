package ur_rna.RNAstructureUI.ui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;

/**
 * @author Richard M. Watson
 */
public class ProgressPanel extends JPanel {
    private final JProgressBar progressBar;
    private final JLabel label;
    private final JButton button;
    private final Timer timer;
    private boolean started;

    public ProgressPanel() { this("Progress: ", null, null); }
    public ProgressPanel(String message, ActionListener cancelListener, ActionListener updateListener) {
        super(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();

        label = new JLabel(message);
        gbc.insets = new Insets(2, 2, 2, 2);
        gbc.gridx = 1;
        gbc.gridy = 0;
        label.setFont(label.getFont().deriveFont(Font.BOLD));
        this.add(label, gbc);

        //Dimension dim = UIManager.getDimension("ProgressBar.horizontalSize");
        //Dialogs.showMessage(dim.toString());

        // Build the progress bar.
        progressBar = new JProgressBar();
        progressBar.setPreferredSize(new Dimension(20, 25));
        //progressBar.setMinimumSize(new Dimension(25,50));
        progressBar.setValue(0);
        gbc.gridx++;
        gbc.weightx = 1;
        gbc.weighty = 1;
        gbc.fill = GridBagConstraints.BOTH;
        this.add(progressBar, gbc);

        button = new JButton("Cancel");
        gbc.fill = GridBagConstraints.NONE;
        gbc.gridx++;
        gbc.weightx = gbc.weighty = 0;
        button.setVisible(cancelListener != null);
        this.add(button, gbc);

        timer = new Timer(250, updateListener);
        timer.setRepeats(true);

        if (cancelListener != null)
            addCancelListener(cancelListener);
    }
    public void start() {
        started = true;
        if (!isIndeterminate())
            timer.start();
    }
    public void stop() {
        started = false;
        timer.stop();
    }

    public void addCancelListener(ActionListener listener) { button.addActionListener(listener); }
    public void removeCancelListener(ActionListener listener) { button.removeActionListener(listener); }
    public ActionListener[] getCancelListeners() { return button.getActionListeners(); }

    public void addUpdateListener(ActionListener listener) { timer.addActionListener(listener); }
    public void removeUpdateListener(ActionListener listener) { timer.removeActionListener(listener); }
    public ActionListener[] getUpdateListeners() { return timer.getActionListeners(); }

    public String getMessage() {
        return label.getText();
    }
    public void setMessage(final String message) {
        label.setText(message);
    }
    public boolean isIndeterminate() {
        return progressBar.isIndeterminate();
    }
    public void setIndeterminate(final boolean indeterminate) {
        progressBar.setIndeterminate(indeterminate);
        if (indeterminate)
            timer.stop();
        else if (started) timer.start();
    }
    public int getProgress() {
        return progressBar.getValue();
    }
    public void setProgress(final int progress) {
        progressBar.setValue(progress);
    }
    public boolean allowCancel() {
        return button.isVisible();
    }
    public void setAllowCancel(final boolean canCancel) {
        this.button.setVisible(canCancel);
    }
//
//        // If the progress bar is indeterminate, explicitly set it as such.
//        // This creates a simple entirely GUI-based dialog that just shows
//        // that a calculation is happening. Otherwise, build a connection to
//        // the back end to monitor and show explicit progress.
//        if( !determinate ) { bar.setIndeterminate( true ); }
//        else {
//
//            // Build the SwingWorker that handles the background monitoring of
//            // the progress bar. The SwingWorker runs a loop that, every half
//            // second, checks the progress of the back end calculation and
//            // updates the progress bar. The "done" method is called once the
//            // calculation is done to dispose of the progress bar.
//            SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {
//                public Void doInBackground() {
//                    try {
//                        setProgress( 0 );
//                        int progress = 0;
//
//                        while( progress < 100) {
//                            try { Thread.sleep( 250 ); }
//                            catch( InterruptedException ex ) {
//                                break;
//                            }
//                            if (progressIndicator == null || self.isClosed())
//                                break;
//                            progress = backend.getProgressNumber();
//                            setProgress( Math.min( progress, 100 ) );
//                        }
//                    } catch( Exception e ) { e.printStackTrace(); }
//
//                    return null;
//                }
//                public void done() {
//                    bar.setIndeterminate(true);
//                }
//            };
//
//            // Add a property change listener to the SwingWorker so the value
//            // of the progress bar updates properly.
//            worker.addPropertyChangeListener( new PropertyChangeListener() {
//                public void propertyChange( PropertyChangeEvent e ) {
//                    if( "progress".equals(e.getPropertyName()) ) {
//                        int progress = (Integer)e.getNewValue();
//                        bar.setValue( progress );
//                    }
//                }
//            });
//
//            // Execute the progress bar through the SwingWorker.
//            worker.execute();
//        }
}
