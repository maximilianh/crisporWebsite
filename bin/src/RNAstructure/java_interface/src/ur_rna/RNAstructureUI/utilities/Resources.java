package ur_rna.RNAstructureUI.utilities;

import ur_rna.RNAstructureUI.RNAstructure;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.Utilities.ResourceLoader;

import java.io.IOException;
import java.io.InputStream;

/**
 * @author Richard M. Watson
 */
public final class Resources {
    private Resources(){} //cannot instantiate
    private static ResourceLoader _res;
    public static ResourceLoader loader() {
        if (_res == null) {
            _res = new ResourceLoader();
            searchForResources(_res);
        }
        return _res;
    }

    public static InputStream get(String name) throws IOException { return loader().getStream(name); }
    public static InputStream tryGet(String name) { return loader().tryGetStream(name); }

    private static void searchForResources(ResourceLoader r) {
        final String resDir = "resources";
        final String resFile = "res.dir";

        ResourceLoader res = loader();
        res.setBaseDir(RNAstructure.class, resDir);
        String[] subDirs = { "images", "sounds" };
        String[] altDirs = {
                ResourceLoader.getDirFromPackage(Resources.class, resDir),
                "ur_rna/" + resDir,
                resDir
        };
        res.addSearchDirs(subDirs);

        for (String alt : altDirs) {
            res.addSearchDir("/" + alt);
            for (String sub : subDirs)
                res.addSearchDir("/" + alt + "/" + sub);
        }
    }

    public static boolean verify() {
        if (loader().hasResource("res.dir"))
            return true;
        Dialogs.showMessage( "Could not locate program resources. The following locations were searched:" +
                String.join("\n\t", loader().getAllSearchDirs()));
        return false;
    }
}
