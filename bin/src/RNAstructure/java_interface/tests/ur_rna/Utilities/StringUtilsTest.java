package ur_rna.Utilities;

import org.junit.Test;

import java.text.ParseException;

import static org.junit.Assert.assertEquals;
import static ur_rna.Utilities.Strings.*;
import static ur_rna.Utilities.testing.TestUtils.assertError;

/**
 * Test the methods of the utility class ur_rna.Utilities.Strings
 */
public class StringUtilsTest {
    @Test
    public void toFriendlyNameTest() throws Exception {
        assertEquals("John", toFriendlyName("john"));
        assertEquals("First Name", toFriendlyName("firstName"));
        assertEquals("Get ID", toFriendlyName("getID"));
        assertEquals("Ssn Number", toFriendlyName("SsnNumber"));
        assertEquals("SSN Number", toFriendlyName("SSNNumber"));
        assertEquals("ID Card", toFriendlyName("IDCard"));
        assertEquals("ID Card", toFriendlyName("ID_Card"));
        assertEquals("X Card", toFriendlyName("XCard"));
        assertEquals("X Card", toFriendlyName("xCard"));
        assertEquals("I Am Nice", toFriendlyName("iAmNice"));
        assertEquals("I Am Nice", toFriendlyName("i am  nice"));
        assertEquals("My Name", toFriendlyName("_myName"));
        assertEquals("My Name", toFriendlyName("__my__Name"));
        assertEquals("Name 1", toFriendlyName("name1"));
        assertEquals("First Name", toFriendlyName("first_name"));
        assertEquals("First Name", toFriendlyName("first__name"));
        assertEquals("First Name", toFriendlyName("First__Name"));
        assertEquals("Game 22 Rules", toFriendlyName("game22rules"));
        assertEquals("Game 22 Rules", toFriendlyName("game_22__Rules"));
        assertEquals("", toFriendlyName("__"));
        assertEquals("", toFriendlyName(" _ _"));
        assertEquals("X", toFriendlyName(" _ x _"));
        assertEquals("USA", toFriendlyName("USA"));
        assertEquals("USA 1", toFriendlyName("USA1"));
        assertEquals("US Am 1", toFriendlyName("USAm1"));
    }

    @Test
    public void unescapeStringLiteralTest() throws Exception {
        // simple sanity-checks
        assertEquals("", unescapeStringLiteral(""));
        assertEquals("a", unescapeStringLiteral("a"));
        assertEquals(" ", unescapeStringLiteral(" "));
        assertEquals("\t", unescapeStringLiteral("\t"));
        assertEquals("\0", unescapeStringLiteral("\0"));
        assertEquals("hello!", unescapeStringLiteral("hello!"));

        // standard single-chars
        assertEquals("\\", unescapeStringLiteral("\\\\"));
        assertEquals("\7", unescapeStringLiteral("\\a"));
        assertEquals("\b", unescapeStringLiteral("\\b"));
        assertEquals("\33", unescapeStringLiteral("\\e"));
        assertEquals("\f", unescapeStringLiteral("\\f"));
        assertEquals("\n", unescapeStringLiteral("\\n"));
        assertEquals("\r", unescapeStringLiteral("\\r"));
        assertEquals("\t", unescapeStringLiteral("\\t"));
        assertEquals("\0", unescapeStringLiteral("\\0"));

        assertEquals("\"", unescapeStringLiteral("\\\""));
        assertEquals("\'", unescapeStringLiteral("\\\'"));

        // unknown sequence
        assertEquals("\\k", unescapeStringLiteral("\\k"));
        // trailing backslash
        assertEquals("\\", unescapeStringLiteral("\\"));
        assertEquals("hi\\", unescapeStringLiteral("hi\\"));

        // printable ASCII
        String ascii = " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
        assertEquals(ascii, unescapeStringLiteral(ascii));
        StringBuilder hiAscii = new StringBuilder(); for(int i = 128; i < 256; i++) hiAscii.append((char)i);
        assertEquals(hiAscii.toString(), unescapeStringLiteral(hiAscii.toString()));

        // adjacent escapes
        assertEquals("\\\n", unescapeStringLiteral("\\\\\n"));
        assertEquals("\r\n", unescapeStringLiteral("\\r\\n"));
        assertEquals("\r\n\t", unescapeStringLiteral("\\r\\n\\t"));
        assertEquals("r\rn\nt\tx", unescapeStringLiteral("r\\rn\\nt\\tx"));

        // control
        assertEquals("\r", unescapeStringLiteral("\\cM"));
        assertEquals("\r\n", unescapeStringLiteral("\\cM\\cJ"));
        assertEquals(";=", unescapeStringLiteral("\\c{\\c}"));
        assertError(()->unescapeStringLiteral("\\c"), ParseException.class);

        // octal
        assertEquals("\3", unescapeStringLiteral("\\3"));
        assertEquals("\37", unescapeStringLiteral("\\37"));
        assertEquals("\377", unescapeStringLiteral("\\377"));
        assertEquals("\378", unescapeStringLiteral("\\378"));
        assertEquals("\37s", unescapeStringLiteral("\\37s"));
        assertEquals("A\37Z", unescapeStringLiteral("A\\37Z"));
        assertEquals("\37\7", unescapeStringLiteral("\\37\\7"));

        // hex
        assertEquals("\r", unescapeStringLiteral("\\x0D"));
        assertEquals("\rE", unescapeStringLiteral("\\x0DE"));
        assertEquals("\r\n", unescapeStringLiteral("\\x0D\\x0A"));
        assertEquals("A\r\nZ", unescapeStringLiteral("A\\x0D\\x0AZ"));
        assertError(()->unescapeStringLiteral("\\x"), ParseException.class);
        assertError(()->unescapeStringLiteral("\\x0"), ParseException.class);
        assertError(()->unescapeStringLiteral("\\x0G"), ParseException.class);
        assertError(()->unescapeStringLiteral("\\xs1"), ParseException.class);

        // unicode
        assertEquals("\r", unescapeStringLiteral("\\u000D"));
        assertEquals("\rD", unescapeStringLiteral("\\u000DD"));
        assertEquals("\r\n", unescapeStringLiteral("\\u000D\\u000A"));
        assertEquals("\u8201", unescapeStringLiteral("\\u8201"));
        assertEquals("\u8Fed", unescapeStringLiteral("\\u8Fed"));
        assertError(()->unescapeStringLiteral("\\u"), ParseException.class);
        assertError(()->unescapeStringLiteral("\\u034"), ParseException.class);
        assertError(()->unescapeStringLiteral("\\u034G"), ParseException.class);
        assertError(()->unescapeStringLiteral("\\us031"), ParseException.class);
    }

    @Test
    public void padTest() throws Exception {
        assertEquals("hi", Strings.padr(1, "hi"));
        assertEquals("hi", padr(2, "hi"));
        assertEquals("hi ", padr(3, "hi"));
        assertEquals("hi--", pad(4, "hi", '-', 0));
        assertEquals("hi", pad(2, "hi there", ' ', 1));
        assertEquals("re", pad(2, "hi there", ' ', -1));

        assertEquals("hi", padl(1, "hi"));
        assertEquals("hi", padl(2, "hi"));
        assertEquals(" hi", padl(3, "hi"));
        assertEquals("--hi", pad(-4, "hi", '-', 0));
        assertEquals("i", pad(-1, "hi", ' ', -1));
        assertEquals("h", pad(-1, "hi", ' ', 1));
    }
}