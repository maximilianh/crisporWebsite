package ur_rna.RNAstructureUI.utilities;

import ur_rna.RNAstructureUI.RNAstructure;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import static ur_rna.Utilities.Strings.isEmpty;

/**
 * @author Richard M. Watson
 */
public class MRUFileStorage {
    private ArrayList<String> list = new ArrayList<>();

    public java.util.List<String> getList() { return Collections.unmodifiableList(list); }

    public void add(String path) {
        list.remove(path);
        list.add(path);
    }
    public void remove(String path) { list.remove(path); }

    public void saveToStorage() {
        Preferences p = RNAstructure.getPrefs().node("MRUFiles");
        int pos = 1;
        for (String s : list)
            if (new File(s).exists())
                p.put("" + (pos++), s);
        try {
            p.flush();
        } catch (BackingStoreException ex) {
            //
        }
    }

    public void loadFromStorage() {
        try {
            list.clear();
            Preferences p = RNAstructure.getPrefs().node("MRUFiles");
            String[] all = p.keys();
            for (String s : all) {
                String path = p.get(s, null);
                if (!isEmpty(path))
                    list.add(path);
            }
        } catch (BackingStoreException ex) {
            //
        }
    }

}
