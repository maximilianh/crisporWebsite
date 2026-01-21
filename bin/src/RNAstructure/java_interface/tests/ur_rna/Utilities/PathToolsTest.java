package ur_rna.Utilities;

import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static ur_rna.Utilities.PathTools.*;

/**
 * Test  ur_rna.Utilities.PathTools including  ur_rna.Utilities.PathTools.FileName
 */
public class PathToolsTest {
    @Test
    public void getBaseName_test1() throws Exception {
        assertEquals("hello", getBaseName("hello.txt"));
        assertEquals("hello", getBaseName("hello."));
        assertEquals("hello", getBaseName("hello"));
        assertEquals("", getBaseName(".txt"));
        assertEquals("hello", getBaseName("/path/to/the.dir/hello.txt"));
        assertEquals("hello", getBaseName("C:\\Users.dir\\hello.txt"));
        assertEquals("hello.txt", getBaseName("hello.txt.bak"));
        assertEquals("hello.txt", getBaseName("/hi.dir/hello.txt.bak"));
    }

    @Test
    public void getBaseName_test2() throws Exception {
        assertEquals("hello", getBaseName("hello.txt", true));
        assertEquals("hello", getBaseName("hello.", true));
        assertEquals("hello", getBaseName("hello", true));
        assertEquals("", getBaseName(".txt", true));
        assertEquals("hello", getBaseName("/path/to/the.dir/hello.txt", true));
        assertEquals("hello", getBaseName("C:\\Users.dir\\hello.txt", true));
        assertEquals("hello", getBaseName("hello.txt.bak", true));
        assertEquals("hello", getBaseName("/hi.dir/hello.txt.bak", true));
    }
    @Test
    public void getExt_test1() throws Exception {
        assertEquals(".txt", getExt("hello.txt"));
        assertEquals(".", getExt("hello."));
        assertEquals("", getExt("hello"));
        assertEquals(".txt", getExt(".txt"));
        assertEquals(".txt", getExt("/path/to/the.dir/hello.txt"));
        assertEquals(".txt", getExt("C:\\Users.dir\\hello.txt"));
        assertEquals(".bak", getExt("hello.txt.bak"));
        assertEquals(".bak", getExt("/hi.dir/hello.txt.bak"));
        assertEquals("", getExt("/hi.dir/hello-txt-bak"));
    }

    @Test
    public void getExt1_test2() throws Exception {
        assertEquals(".txt", getExt("hello.txt", true, true));
        assertEquals(".", getExt("hello.", true, true));
        assertEquals("", getExt("hello", true, true));
        assertEquals(".txt", getExt(".txt", true, true));
        assertEquals(".txt", getExt("/path/to/the.dir/hello.txt", true, true));
        assertEquals(".txt", getExt("C:\\Users.dir\\hello.txt", true, true));
        assertEquals(".txt.bak", getExt("hello.txt.bak", true, true));
        assertEquals(".txt.bak", getExt("/hi.dir/hello.txt.bak", true, true));
        assertEquals("", getExt("/hi.dir/hello-txt-bak", true, true));
    }

    @Test
    public void getDir_test1() throws Exception {
        assertEquals("", getDir(".txt"));
        assertEquals("/path/to/the.dir/", getDir("/path/to/the.dir/hello.txt"));
        assertEquals("C:\\Users.dir\\", getDir("C:\\Users.dir\\hello.txt"));
        assertEquals("", getDir("hello.txt.bak"));
        assertEquals("/hi.dir/", getDir("/hi.dir/hello.txt.bak"));
        assertEquals("/", getDir("/hello-txt-bak"));
        assertEquals("C:\\", getDir("C:\\hello-txt-bak"));
        assertEquals("C:\\Users\\", getDir("C:\\Users\\hello-txt-bak"));
        assertEquals("C:\\", getDir("C:\\"));
    }

    @Test
    public void getDir_test2() throws Exception {
        assertEquals("", getDir(".txt", false));
        assertEquals("/path/to/the.dir", getDir("/path/to/the.dir/hello.txt", false));
        assertEquals("C:\\Users.dir", getDir("C:\\Users.dir\\hello.txt", false));
        assertEquals("", getDir("hello.txt.bak", false));
        assertEquals("/hi.dir", getDir("/hi.dir/hello.txt.bak", false));
        assertEquals("", getDir("/hello-txt-bak", false));
        assertEquals("C:", getDir("C:\\hello-txt-bak", false));
        assertEquals("C:\\Users", getDir("C:\\Users\\hello-txt-bak", false));
        assertEquals("C:", getDir("C:\\", false));
    }


    @Test
    public void parse_test1() throws Exception {
        FileName fn = PathTools.parse("hello.txt");
        assertEquals("", fn.dir());
        assertEquals("hello.txt", fn.name());
        assertEquals(".txt", fn.ext());
        assertEquals("hello", fn.baseName());
        assertArrayEquals(new String[0], fn.splitDirs());

        fn = PathTools.parse("/root/hello.txt");
        assertEquals("/root/", fn.dir());
        assertEquals("hello.txt", fn.name());
        assertEquals(".txt", fn.ext());
        assertEquals("hello", fn.baseName());
        assertArrayEquals(new String[] {"","root"}, fn.splitDirs());

        fn = PathTools.parse("/root.dir/hello.bak.txt");
        assertEquals("/root.dir/", fn.dir());
        assertEquals("hello.bak.txt", fn.name());
        assertEquals(".txt", fn.ext());
        assertEquals("hello.bak", fn.baseName());
        assertArrayEquals(new String[] {"","root.dir"}, fn.splitDirs());
    }

    @Test
    public void FileName_splitDirs_test() throws Exception {
        assertArrayEquals(new String[] {"","root.dir", "hello", "bak", "txt" }, PathTools.parse("/root.dir/hello/bak/txt/").splitDirs());
        assertArrayEquals(new String[] {"","root.dir", "hello", "bak" }, PathTools.parse("/root.dir/hello/bak/txt").splitDirs());
        assertArrayEquals(new String[] { "" }, PathTools.parse("/root").splitDirs());
        assertArrayEquals(new String[] { "" }, PathTools.parse("/").splitDirs());
        assertArrayEquals(new String[] { }, PathTools.parse("file.txt").splitDirs());
    }

    @Test
    public void FileName_changeExt_test() throws Exception {
        assertEquals(".txt", PathTools.parse("/path-to/hello.txt").ext());
        assertEquals("hello", PathTools.parse("/path-to/hello.txt").baseName());
        assertEquals("/path-to/hello", PathTools.parse("/path-to/hello.txt").removeExt());
        assertEquals("/path-to/hello", PathTools.parse("/path-to/hello.txt").changeExt(null));

        assertEquals("hello.zip", PathTools.parse("hello.txt").changeExt("zip"));
        assertEquals("/path-to/hello.zip", PathTools.parse("/path-to/hello.txt").changeExt("zip"));
        assertEquals("hello.zip", PathTools.parse("hello.txt").changeExt(".zip"));
        assertEquals("/path-to/hello.zip", PathTools.parse("/path-to/hello.txt").changeExt(".zip"));
        assertEquals("hello.zip", PathTools.parse("hello").changeExt("zip"));
        assertEquals("/path-to/hello.zip", PathTools.parse("/path-to/hello").changeExt("zip"));
        assertEquals("hello.zip", PathTools.parse("hello").changeExt(".zip"));
        assertEquals("/path-to/hello.zip", PathTools.parse("/path-to/hello").changeExt(".zip"));
        assertEquals("/path-to/hello", PathTools.parse("/path-to/hello.txt").changeExt(""));
    }

    @Test
    public void parse_test2() throws Exception {
        FileName fn = PathTools.parse("/root.dir/hello.bak.txt", true);
        assertEquals("/root.dir/", fn.dir());
        assertEquals("hello.bak.txt", fn.name());
        assertEquals(".bak.txt", fn.ext());
        assertEquals("hello", fn.baseName());
        assertArrayEquals(new String[] {"","root.dir"}, fn.splitDirs());
    }

    @Test
    public void addSlash_test() throws Exception {
        assertEquals("/root.dir\\", addSlash("/root.dir"));
        assertEquals("/root.dir/file.txt\\", addSlash("/root.dir/file.txt"));
        assertEquals("\\", addSlash(""));
    }
}