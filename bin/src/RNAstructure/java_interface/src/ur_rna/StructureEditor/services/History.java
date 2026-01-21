package ur_rna.StructureEditor.services;


import java.util.ArrayList;
import java.util.List;

/**
 * Stores the history of an arbitrary model for undo/redo
 */
public class History<T> {
    private ArrayList<T> stack = new ArrayList<>();
    private int pos = -1;
    /**
     * Store the current state of the model.
     * Typically this clears all history in front of the current entry.
     * @param state The current state of the model.
     */
    public void store(T state) {
        if (canRedo())
            clearRedo();
        stack.add(state);
        pos++;
    }

    public List<T> entries() { return stack; }

    public int size() { return stack.size();}
    public int position() { return pos;}
    public T get(int position) { return stack.get(position); }
    public T current() {
        if (pos == -1) //should only occur if size() == 0
            throw new IndexOutOfBoundsException();
        return stack.get(pos);
    }

    /**
     * Clear all entries in front of the current position.
     */
    private void clearRedo() {
        for (int i = stack.size() - 1; i > pos; i--)
            stack.remove(i);
    }

    public void clear() {
        stack.clear();
        pos = -1;
    }

    public T undo() {
        if (!canUndo())
            throw new IndexOutOfBoundsException();
        return stack.get(--pos);
    }
    public T redo() {
        if (!canRedo())
            throw new IndexOutOfBoundsException();
        return stack.get(++pos);
    }

    public boolean canUndo() {
        return pos > 0;
    }
    public boolean canRedo() {
        return pos < stack.size() - 1;
    }
}
