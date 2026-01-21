package ur_rna.StructureEditor.models;

/**
 *
 */
public class HistoryUpdateEvent {
    public final Object source;
    public final boolean storing;
    public final HistoryState state;
    public HistoryUpdateEvent(final Object source, final boolean storing, final HistoryState state) {
        this.source = source;
        this.storing = storing;
        this.state = state;
    }
}
