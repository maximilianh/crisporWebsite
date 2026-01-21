package ur_rna.StructureEditor.services;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author Richard M. Watson
 */
public class HistoryTest {
    History<String> history;

    @Before
    public void setUp() throws Exception {
        history = new History<>();
    }
    @Test
    public void testStore() throws Exception {
        history.store("first");
        assertEquals(1, history.size());
        assertEquals(0, history.position());
        assertEquals(false, history.canUndo());
        assertEquals(false, history.canRedo());
        history.store("second");
        assertEquals(2, history.size());
        assertEquals(1, history.position());
        assertEquals(true, history.canUndo());
        assertEquals(false, history.canRedo());
    }
    @Test
    public void testUndoSingleItem() throws Exception {
        history.store("first");
        assertEquals(0, history.position());
        assertEquals(1, history.size());
        assertEquals(false, history.canUndo());
        assertEquals(false, history.canRedo());
        history.store("second");
        assertEquals("first", history.undo());
        assertEquals(0, history.position());
        assertEquals(2, history.size());
        assertEquals(false, history.canUndo());
        assertEquals(true, history.canRedo());
        assertEquals("second", history.redo());
        assertEquals(true, history.canUndo());
        assertEquals(false, history.canRedo());
    }
    @Test
    public void testUndoSecondItem() throws Exception {
        history.store("first");
        history.store("second");
        assertEquals(1, history.position());
        assertEquals(2, history.size());
        assertEquals("first", history.undo());
        assertEquals("first", history.current());
        assertEquals(0, history.position());
        assertEquals(2, history.size());
        assertEquals(false, history.canUndo());
        assertEquals(true, history.canRedo());
        assertEquals("second", history.redo());
        assertEquals("second", history.current());
        assertEquals(1, history.position());
        assertEquals(2, history.size());
        assertEquals(true, history.canUndo());
        assertEquals(false, history.canRedo());
    }
    @Test(expected = IndexOutOfBoundsException.class)
    public void testUndoFirstItem() throws Exception {
        history.store("first");
        history.undo();
    }
    @Test(expected = IndexOutOfBoundsException.class)
    public void testRedoFinalItem() throws Exception {
        history.store("first");
        history.redo();
    }
    @Test
    public void testRedoSecondItem() throws Exception {
        history.store("first");
        history.store("second");
        assertEquals("first", history.undo());
        assertEquals("second", history.redo());
    }

    @Test(expected = IndexOutOfBoundsException.class)
    public void testRedoFinalItem2() throws Exception {
        history.store("first");
        history.store("second");
        history.undo();
        history.redo();
        history.redo();
    }
}