package ur_rna.Utilities.prefs;

import org.junit.Test;

import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

/**
 * @author Richard M. Watson
 */
public class FilePreferencesTest {
    @Test
    public void testPrefs() throws Exception {
        // Runtime.getRuntime().addShutdownHook(hook);

        FilePreferences.Factory.register();
        System.setProperty("program.title", "PrefsTest");

        Preferences root = Preferences.userRoot();
        Preferences p = root;

        System.out.println("Preferences File:" + ((FilePreferences)root).getPreferencesFile());
        System.out.println("### Initial settings ###");
        printAll("", root);

        p.putBoolean("testBool", true);
        p.put("testNumber", String.valueOf(System.currentTimeMillis()));

        p = p.node("settings");

        p.putByteArray("bananas", new byte[] { 1, 2, 4, 7, 21, 0xF, Byte.MAX_VALUE, (byte)255 });
        p.put("fred", "nancy");

        System.out.println("### Before remove settings ###");
        printAll("", root);
        p.removeNode();

        p = Preferences.userRoot().node("settings");

        p.put("hi there", "tango");

        p = p.node("level2");
        p.put("fred2", "nancy2");
        p.put("testDate", new java.util.Date().toString());

        System.out.println("### Before clear level2 ###");
        printAll("", root);
        p.clear();

        System.out.println("### Final ###");
        printAll("", root);
    }

    static void printAll(String level, Preferences p) throws BackingStoreException {
        System.out.println(level + p.name() + "/");

        for (String s : p.keys())
            System.out.println(level + "•" + s + "=" + p.get(s, null));

        final String nextLevel = level + "  ";
        for (String s : p.childrenNames())
            printAll(nextLevel, p.node(s));
    }
}