package ur_rna.RNAstructureUI.utilities;

import ur_rna.Utilities.EventArgs;
import ur_rna.Utilities.EventSource;
import ur_rna.Utilities.Strings;

import javax.swing.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Richard M. Watson
 */
public abstract class BackgroundWorker extends SwingWorker<Void, Void>  implements PropertyChangeListener {
    public final static int DEFAULT_UPDATE_INTERVAL = 250;
    protected String status = "";
    protected boolean progressIndeterminate;

    //    private Runnable work;
    private final ArrayList<ErrorInfo> errors = new ArrayList<>();
    public final EventSource.TwoArgs<BackgroundWorker, Integer> ProgressChange = new EventSource.TwoArgs<>();
    public final EventSource.TwoArgs<BackgroundWorker, StateChangeEventArgs> StateChange = new EventSource.TwoArgs<>();
    public final EventSource.TwoArgs<BackgroundWorker,String> StatusChange = new EventSource.TwoArgs<>();
    public final EventSource.OneArg<BackgroundWorker> WorkDone = new EventSource.OneArg<>();

    public static class StateChangeEventArgs extends EventArgs {
        public final StateValue old;
        public final StateValue state;
        public StateChangeEventArgs(final StateValue old, final StateValue state) {
            this.old = old;
            this.state = state;
        }
    }

    public BackgroundWorker() {
        addPropertyChangeListener(this);
    }

    public void setStatus(String status) {
        String old = this.status;
        this.status = status;
        firePropertyChange("status", old, status);
    }

    public void setProgressIndeterminate(boolean isIndeterminate) {
        progressIndeterminate = isIndeterminate;
        int progress = getProgress();
        firePropertyChange("progress", progress, -1);
    }
    public boolean isProgressIndeterminate() {
        return progressIndeterminate;
    }

    public String getStatus() {
        return status;
    }

    @Override // Implements PropertyChangeListener
    public void propertyChange(final PropertyChangeEvent evt) {
        switch(evt.getPropertyName()) {
            case "state":
                StateChangeEventArgs ea = new StateChangeEventArgs(
                        (StateValue)evt.getOldValue(),
                        (StateValue)evt.getNewValue()
                );
                stateChanged(ea.old, ea.state);
                StateChange.invoke(this, ea);
                if (StateValue.DONE.equals(ea.state))
                    WorkDone.invoke(this);
                break;
            case "progress":
                progressChanged((Integer) evt.getNewValue());
                break;
            case "status":
                statusChanged((String)evt.getNewValue());
                break;
        }
    }

    protected void progressChanged(final int percent) {
        ProgressChange.invoke(this, percent);
    }
    protected void stateChanged(final StateValue oldState, final StateValue state) {
        StateChange.invoke(this, new StateChangeEventArgs(oldState, state));
    }
    protected void statusChanged(final String status) {
        StatusChange.invoke(this, status);
    }

    /**
     * Receives data chunks from the {@code publish} method asynchronously on the
     * <i>Event Dispatch Thread</i>.
     * <p>
     * <p>
     * Please refer to the {@link #publish} method for more details.
     *
     * @param chunks intermediate results to process
     * @see #publish
     */
    @Override
    protected void process(final List<Void> chunks) {
        super.process(chunks);
        for(Object o : chunks)
            processOutput(o);
    }
    protected void processOutput(final Object o) {   }

    public static class ErrorInfo {
        public String message;
        public Exception ex;
    }

    public void setError(Exception ex) { setError(null, ex); }
    public void setError(String message) { setError(message, null); }
    public void setError(String message, Exception ex) {
        ErrorInfo err = new ErrorInfo();
        err.message = message;
        err.ex = ex;
        errors.add(err);
    }

    public boolean hadErrors() { return errors.size() != 0; }
    public List<ErrorInfo> getErrors() { return errors; }
    public String getErrorMessages() {
        StringBuilder sb = new StringBuilder();
        int i = 0;
        for (ErrorInfo ei : this.errors) {
            sb.append(String.format("%4d.", ++i));
            sb.append("  ");
            if (!Strings.isEmpty(ei.message))
                sb.append(ei.message);
            else if (ei.ex != null)
                sb.append(ei.ex.getMessage());
            sb.append("\n");
        }
        return sb.toString();
    }

    @Override
    protected Void doInBackground() {
        try {
            work();
        }catch (Exception ex) {
            setError("Uncaught error while running background calculation.", ex);
        }
        return null;
    }

    public abstract void work();
}
