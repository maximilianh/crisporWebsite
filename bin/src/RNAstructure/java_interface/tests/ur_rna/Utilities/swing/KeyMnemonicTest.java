package ur_rna.Utilities.swing;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static ur_rna.Utilities.swing.KeyMnemonic.*;

/**
 * @author Richard M. Watson
 */
public class KeyMnemonicTest {

    @Test
    public void testGetMnemonic() throws Exception {
        assertEquals(0, getMnemonic("File"));
        assertEquals('F', getMnemonic("&File"));
        assertEquals('I', getMnemonic("F&ile"));
        assertEquals('E', getMnemonic("Fil&e"));
        assertEquals(0, getMnemonic("File&"));
        assertEquals('J', getMnemonic("PB && &Jelly"));
        assertEquals('J', getMnemonic("PB&&&Jelly"));
        assertEquals(0, getMnemonic("PB&&&&Jelly"));
        assertEquals('J', getMnemonic("P&&B&J&T"));
        assertEquals('T', getMnemonic("P&&B&&J&T"));
        assertEquals('B', getMnemonic("P&B&J&T"));
        assertEquals('B', getMnemonic("&&P&B&J&T"));
        assertEquals('P', getMnemonic("&&&P&B&J&T"));
    }
    @Test
    public void testStripMnemonics() throws Exception {
        assertEquals("File", stripMnemonics("File"));
        assertEquals("File", stripMnemonics("&File"));
        assertEquals("File", stripMnemonics("F&ile"));
        assertEquals("File", stripMnemonics("Fil&e"));
        assertEquals("File", stripMnemonics("File&"));
        assertEquals("PB & Jelly", stripMnemonics("PB && &Jelly"));
        assertEquals("PB&Jelly", stripMnemonics("PB&&&Jelly"));
        assertEquals("PB&&Jelly", stripMnemonics("PB&&&&Jelly"));
        assertEquals("P&BJ&T", stripMnemonics("P&&B&J&T"));
        assertEquals("P&B&JT", stripMnemonics("P&&B&&J&T"));
        assertEquals("PB&J&T", stripMnemonics("P&B&J&T"));
        assertEquals("PB&J&T&", stripMnemonics("P&B&J&T&"));
        assertEquals("P&B&J&T&", stripMnemonics("&P&B&J&T&"));
        assertEquals("&PB&J&T", stripMnemonics("&&P&B&J&T"));
        assertEquals("&P&B&J&T", stripMnemonics("&&&P&B&J&T"));
        assertEquals("", stripMnemonics("&"));
        assertEquals("&", stripMnemonics("&&"));
        assertEquals("&", stripMnemonics("&&&"));
        assertEquals("&T", stripMnemonics("&&&T"));
        assertEquals("&&", stripMnemonics("&&&&"));
    }
    @Test
    public void testGetMnemonicIndex() throws Exception {
        assertEquals(-1, getMnemonicIndex("File"));
        assertEquals(0, getMnemonicIndex("&File"));
        assertEquals(1, getMnemonicIndex("F&ile"));
        assertEquals(3, getMnemonicIndex("Fil&e"));
        assertEquals(-1, getMnemonicIndex("File&"));
        assertEquals(5, getMnemonicIndex("PB && &Jelly"));
        assertEquals(3, getMnemonicIndex("PB&&&Jelly"));
        assertEquals(-1, getMnemonicIndex("PB&&&&Jelly"));
        assertEquals(3, getMnemonicIndex("P&&B&J&T"));
        assertEquals(5, getMnemonicIndex("P&&B&&J&T"));
        assertEquals(1, getMnemonicIndex("P&B&J&T"));
        assertEquals(1, getMnemonicIndex("P&B&J&T&"));
        assertEquals(0, getMnemonicIndex("&P&B&J&T&"));
        assertEquals(2, getMnemonicIndex("&&P&B&J&T"));
        assertEquals(1, getMnemonicIndex("&&&P&B&J&T"));
        assertEquals(-1, getMnemonicIndex("&"));
        assertEquals(-1, getMnemonicIndex("&&"));
        assertEquals(-1, getMnemonicIndex("&&&"));
        assertEquals(1, getMnemonicIndex("&&&T"));
        assertEquals(-1, getMnemonicIndex("&&&&"));
        assertEquals(12, getMnemonicIndex("&& P && B && J & "));
        assertEquals(3, getMnemonicIndex("&&&&&&&T"));
    }
}