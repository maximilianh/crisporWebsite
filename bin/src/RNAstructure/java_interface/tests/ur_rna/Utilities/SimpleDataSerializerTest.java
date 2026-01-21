package ur_rna.Utilities;

import org.junit.Before;
import org.junit.Test;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.junit.Assert.*;
import static ur_rna.Utilities.SimpleDataSerializer.*;
import static ur_rna.Utilities.Strings.escapeStringLiteral;
import static ur_rna.Utilities.testing.TestUtils.assertInstanceOf;

public class SimpleDataSerializerTest {
    private SimpleDataSerializer p;

    @Before
    public void setUp() {
        p = new SimpleDataSerializer();
    }

    @Test
    public void setUpTest() throws Exception {
        assertNotNull(p);
    }

    @Test
    public void regexText() throws Exception {
        Pattern p = Pattern.compile("/-(?:[^-]|-(?!/))++-/");
        Matcher m = p.matcher("000/-aaa/bbb/ccc-/ddd/eee-/");
        assertTrue(m.find());
        assertEquals("/-aaa/bbb/ccc-/", m.group(0));
    }

    @Test
    public void tokenizeTest() throws Exception {
        String textMLC = "/* multi-line comment  \" {} [] !@#$%^&*()_/*\r\n/?<>\n\n / */";
        String textSLC = "// single-line comment \" {} [] !@#$%^&*()_";
        String textTypes = "true false null 5 ID1 \"fred\" \"joe is good\" -6.25E-21 [ ] { } , : " + textMLC + " " + textSLC + "\n THE_END";

        Token[] parsed = p.tokenize(textTypes, null, true, false);
        Object[] expected = {
                TokenPattern.BOOL, "true",
                TokenPattern.BOOL, "false",
                TokenPattern.NULL, "null",
                TokenPattern.NUMBER, "5",
                TokenPattern.IDENTIFIER, "ID1",
                TokenPattern.STRING,"\"fred\"",
                TokenPattern.STRING,"\"joe is good\"",
                TokenPattern.NUMBER,"-6.25E-21",
                TokenPattern.OPEN_BRACKET, "[",
                TokenPattern.CLOSE_BRACKET, "]",
                TokenPattern.OPEN_BRACE, "{",
                TokenPattern.CLOSE_BRACE, "}",
                TokenPattern.COMMA,  ",",
                TokenPattern.COLON, ":",
                TokenPattern.COMMENT, textMLC,
                TokenPattern.COMMENT, textSLC,
                TokenPattern.IDENTIFIER, "THE_END"
        };
        for (int i = 0; i < parsed.length; i++) {
            assertEquals(expected[i * 2], parsed[i].type());
            assertEquals(expected[i * 2 + 1], parsed[i].text());
        }
    }

    @Test
    public void parseObjTest() throws Exception {
        String textObj = "{ name: \"Apples are nice.\", name2: 3.45, \"3-name\": true }";
        Object o = p.parse(textObj, null);
        assertInstanceOf(Map.class, o);
        Map obj = (Map)o;
        assertEquals("Apples are nice.", obj.get("name"));
        assertEquals(3.45, obj.get("name2"));
        assertEquals(true, obj.get("3-name"));
        assertNull(obj.get("noname"));
    }

    @Test
    public void parseArrTest() throws Exception {
        String textArr = "[ \"fred\", \"joe is good\", -6.25E-21, false ]";
        Object o = p.parse(textArr, null);
        assertInstanceOf(List.class, o);
        List obj = (List)o;
        assertEquals("fred", obj.get(0));
        assertEquals("joe is good", obj.get(1));
        assertEquals(-6.25E-21, obj.get(2));
        assertEquals(false, obj.get(3));
        assertEquals(4, obj.size());
    }

    @Test
    public void parseArrNestedTest() throws Exception {
        String textArr = "[ \"fred\", { tulips: false, oranges: \"yes\", arr: [ 1, \"s\"] }, [ \"Yum! Banana\", null, +6.25E+21], false ]";

        Object parsed = p.parse(textArr, null);
        assertInstanceOf(List.class, parsed);

        List a1 = (List)parsed;
        assertEquals("fred", a1.get(0));
        //assertEquals("joe is good", arr.getValue(1));
        //assertEquals(-6.25E-21, obj.getValue(2));
        assertEquals(false, a1.get(3));
        assertEquals(4, a1.size());

        assertInstanceOf(Map.class, a1.get(1));
        Map o2 = (Map)a1.get(1);
        assertEquals(false, o2.get("tulips"));
        assertEquals("yes", o2.get("oranges"));
        assertEquals(3, o2.size());

        assertInstanceOf(List.class, a1.get(2));
        List a2 = (List)a1.get(2);
        assertEquals("Yum! Banana", a2.get(0));
        assertEquals(null, a2.get(1));
        assertEquals(+6.25e+21, a2.get(2)); //+6.25e+21
        assertEquals(3, a2.size());

        assertInstanceOf(List.class, o2.get("arr"));
        List a3 = (List)o2.get("arr");
        assertEquals(1d, a3.get(0));
        assertEquals("s", a3.get(1));
        assertEquals(2, a3.size());
    }

    @Test
    public void tokenizeErrorTests() throws Exception {
        String[] errorStrs = {
                "/* un-terminated multiline-line comment",
                "// single-line comment after newline \n $%^&*()_",
                "hello \"unterminated String",
                "-word",
                "\0",
                "\u001F", // ASCII 31
                "6frog",
                "6+5",
                "frog+6",
        };
        int line = 0;
        for(String s : errorStrs) {
            line++;
            try {
                p.tokenize(s,null, false, false);
                throw new AssertionError("tokenize failed to throw the expected exception. text: " + escapeStringLiteral(s) + " on test " + line);
            } catch (SyntaxErrorException ex) {
                // success
            }
            //assertError(()->p.tokenize(s,null, false, false), SyntaxErrorException.class);
        }
    }

    @Test
    public void parseErrorTests() throws Exception {
        String[] errorStrs = {
                "adjacent words",
                "\"adjacent\" \"strings\"",
                "5 -6",
                "[ unterminated array ",
                " ]",
                "{ unterminated: obj ",
                " }",
                "{ name1: \"obj\"  no-colon }",
                "{ name1: \"obj\" , , two-commas: 3 }",
                "[ a , , b ]",
                "6 , 3",
                ",",
                "name: \"prop\"",
        };
        int line = 0;
        for(String s : errorStrs) {
            line++;
            try {
                p.parse(s,null);
                throw new AssertionError("parse failed to throw the expected exception. text: " + escapeStringLiteral(s) + " on test " + line);
            } catch (SyntaxErrorException ex) {
                System.out.println(ex.getMessage());
            }
            //assertError(()->p.tokenize(s,null, false, false), SyntaxErrorException.class);
        }
    }

    @Test
    public void writeTests() throws Exception {
        OutputFormatSettings s = new OutputFormatSettings();
        // Start with the most "terse" settings
        s.indent = 0;
        s.quoteAlways = false;
        s.quoteProps = false;
        s.inlineLists = true;
        s.inlineMaps = true;
        s.commasInLists = false;
        s.commasInMaps = false;
        s.spaceAfterProp = false;

        assertEquals(OutputFormatSettings.SSON, s);
        // Define and set up some OutputFormatSettings variations.
        OutputFormatSettings sUseQuotes = s.clone(), sUseCommas = s.clone(), sUseIndent = s.clone(),
                sUseSpacesInMaps = s.clone(), sTmp, sNoInline = s.clone();


        sUseQuotes.quoteAlways = sUseQuotes.quoteProps = true;
        sUseIndent.indent = 2;
        sUseCommas.commasInMaps = sUseCommas.commasInLists = true;
        sUseSpacesInMaps.spaceAfterProp = true;
        sNoInline.inlineMaps = sNoInline.inlineLists = false; sNoInline.indent = 2; // No inline lists. Do use an indent.

        assertEquals("6", p.toString(6, s));
        assertEquals("6.5", p.toString(6.5, s));
        assertEquals("0.5", p.toString(1f/2, s));
        //assertEquals("0.57", p.toString(.566, sNum));

        assertEquals("true", p.toString(Boolean.TRUE, s));

        assertEquals("null", p.toString(null, s));

        assertEquals("fred", p.toString("fred", s));
        assertEquals("\"fred\"", p.toString("fred", sUseQuotes));
        assertEquals("\"Hello World!\"", p.toString("Hello World!", s));
        assertEquals("\"Hello\\nWorld!\"", p.toString("Hello\nWorld!", s));

        List<Object> list = Arrays.asList(1.20E+7, "fred", "hi\tyou", false, null, 3);
        assertEquals("[ 1.2E7 fred \"hi\\tyou\" false null 3 ]", p.toString(list, s));
        assertEquals("[ 1.2E7, fred, \"hi\\tyou\", false, null, 3 ]", p.toString(list, sUseCommas));

        Map<Object, Object> m = new LinkedHashMap<>();
        m.put("name", "fred");
        m.put("say", "Hello World!");
        m.put("age", 65.2);
        m.put(25, true);

        assertEquals("{ name:fred say:\"Hello World!\" age:65.2 \"25\":true }", p.toString(m, s));
        assertEquals("{ name: fred say: \"Hello World!\" age: 65.2 \"25\": true }", p.toString(m, sUseSpacesInMaps));
        assertEquals("{ name:fred, say:\"Hello World!\", age:65.2, \"25\":true }", p.toString(m, sUseCommas));
        sTmp = sUseCommas.clone(); sTmp.spaceAfterProp = true;
        assertEquals("{ name: fred, say: \"Hello World!\", age: 65.2, \"25\": true }", p.toString(m, sTmp));
        assertEquals("{ name:fred say:\"Hello World!\" age:65.2 \"25\":true }", p.toString(m, sUseIndent));
        assertEquals(LF("{\n"+
                        "  name:fred\n"+
                        "  say:\"Hello World!\"\n"+
                        "  age:65.2\n"+
                        "  \"25\":true\n"+
                        "}"),
                p.toString(m, sNoInline));

        list = new ArrayList<>(list); // allow changes to list.
        list.add(3, new LinkedHashMap<>(m));
        m.put("list", list);
        assertEquals(LF("{\n"+
                        "  name:fred\n"+
                        "  say:\"Hello World!\"\n"+
                        "  age:65.2\n"+
                        "  \"25\":true\n"+
                        "  list:[\n"+
                        "    1.2E7\n"+
                        "    fred\n"+
                        "    \"hi\\tyou\"\n"+
                        "    {\n"+
                        "      name:fred\n"+
                        "      say:\"Hello World!\"\n"+
                        "      age:65.2\n"+
                        "      \"25\":true\n"+
                        "    }\n"+
                        "    false\n"+
                        "    null\n"+
                        "    3\n"+
                        "  ]\n"+
                        "}"),
                p.toString(m, sNoInline));
        list.remove(3);
        assertEquals(LF("{\n"+
                        "  name:fred\n"+
                        "  say:\"Hello World!\"\n"+
                        "  age:65.2\n"+
                        "  \"25\":true\n"+
                        "  list:[ 1.2E7 fred \"hi\\tyou\" false null 3 ]\n"+
                        "}"),
                                p.toString(m, sUseIndent));
    }
    private static String LF(String s) {
        return s.replace("\n", System.lineSeparator());
    }
}