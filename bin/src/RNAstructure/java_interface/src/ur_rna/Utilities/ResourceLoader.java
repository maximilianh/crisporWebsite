package ur_rna.Utilities;

import ur_rna.Utilities.annotation.NotNull;
import ur_rna.Utilities.annotation.Nullable;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;
import java.util.List;

import static ur_rna.Utilities.PathTools.getCanonicalPath;

/**
 * @author Richard M. Watson
 *
 * A simple class used to load resources contained within a jar file
 * in the class-path of the application (or alternatively from the file system).
 *
 * The purpose of this class is to make it so that client code can request
 * a resource by name, without needing to specify the path
 * of the actual resource file. I.e. the client code does not need to know the full
 * path to the resources or the details of retrieving them.
 *
 * To achieve this, the ResourceLoader is initialized with a base directory and
 * optionally one or more sub-directories to searched for each requested resource.
 *
 * The base directory can be specified using either an absolute path, or a Class.
 * If the latter approach is used, the package name of the Class is used to
 * determine the base directory by converting each dot (.) in the package name to a slash (/).
 * In addition, a sub-directory can be appended after the package-derived path.
 * E.g. the call {@code setBaseDir(my.package.MyClass.class, "resources") } would set the
 * resource base directory to {@code my/package/resources}.
 *
 * TODO: Make it so that resources can be requested WITHOUT specifying the file extension. E.g. by associating one or more extensions with each sub-directory.
 *
 */
public class ResourceLoader {
    protected ClassLoader loader;
    protected String baseDir;
    protected List<String> searchDirs = new ArrayList<>();
    protected static final char pathSep = '/'; //File.separatorChar
    private List<String> allSearchDirs;

    public ResourceLoader() { this("resources/"); }
    public ResourceLoader(String resourceBaseDir) {
        setBaseDir(resourceBaseDir);
        setLoader(null);
    }
    public ResourceLoader(String resourceBaseDir, ClassLoader loader) {
        this(resourceBaseDir);
        setLoader(loader);
    }
    /**
     * Create a new ResourceLoader and set its resource directory based on the
     * package name of the specified class by converting each dot (.) in the package
     * name to a slash (/). If sub-directory is also specified, it will be appended to
     * the package-derived directory.
     * E.g. {@code ResourceLoader(my.package.MyClass.class, "resources")} would set the
     * resource base directory to {@code my/package/resources}.
     *
     * @param resourcePeerClass The Class from which to obtain the package name.
     * @param pathRelativeToClassPackage This specifies a sub-directory relative to the
     * package-derived directory from which resources can be obtained. (This can be
     * null or empty, in which case resources will be served directly from the
     * package-derived directory).
     */
    public ResourceLoader(Class resourcePeerClass, String pathRelativeToClassPackage) {
        this(getDirFromPackage(resourcePeerClass, pathRelativeToClassPackage));
    }
    public void setLoader(@Nullable ClassLoader c) {
        if (c == null)
            loader = ClassLoader.getSystemClassLoader();
        else
            loader = c;
    }
//    public void setLoaderFrom(@NotNull Class<?> c) {
//        if (c == null) throw new NullPointerException("Class cannot be null.");
//        setLoader(c.getClassLoader());
//    }
//    public void setLoaderFrom(@NotNull Object o) {
//        if (o == null) throw new NullPointerException("Object cannot be null.");
//        setLoaderFrom(o.getClass());
//    }

    public void setBaseDir(String pathRelativeToClassPath) {
        baseDir = pathRelativeToClassPath;
        if (baseDir == null)
            baseDir = "";
        else {
            baseDir = getCanonicalPath(pathRelativeToClassPath, true, false, true, true);
            if (baseDir.equals("/")) baseDir = "";
        }
    }


    /**
     * Set the resource directory based on the package name of the specified class
     * by converting each dot (.) in the package name to a slash (/). If a non-empty
     * sub-directory is also specified, it will be appended to the package-derived
     * directory.
     * E.g. {@code ResourceLoader(my.package.MyClass.class, "resources")} would set the
     * resource base directory to {@code my/package/resources}.
     *
     * @param resourcePeerClass The Class from which to obtain the package name.
     * @param pathRelativeToClassPackage This specifies a sub-directory relative to the
     * package-derived directory from which resources can be obtained. (This can be
     * null or empty, in which case resources will be served directly from the
     * package-derived directory).
     */
    public void setBaseDir(@NotNull Class<?> resourcePeerClass, String pathRelativeToClassPackage) {
        setBaseDir(getDirFromPackage(resourcePeerClass, pathRelativeToClassPackage));
    }
//    public void setBaseDir(@NotNull Object obj, String pathRelativeToClassPackage) {
//        setBaseDir(obj.getClass(), pathRelativeToClassPackage);
//    }
    public static String getDirFromPackage(@NotNull Class<?> c, String relativeSubDir) {
        return getDirFromPackage(c.getPackage(), relativeSubDir);
    }
    public static String getDirFromPackage(@NotNull Package p, String relativeSubDir) {
        String path = p.getName().replace('.', pathSep);
        if (relativeSubDir != null && relativeSubDir.length() != 0)
            path += pathSep + relativeSubDir;
        return path;
    }

    public String getBaseDir() {
        return baseDir;
    }

    /**
     * Add a sub-directory (relative to the resource base directory) that will be
     * searched to find requested resources.
     * @param dir The subdirectory to add. If the directory starts with a slash (/)
     *            it indicates that it is an absolute directory instead of being
     *            relative to the resource base directory.
     */
    public void addSearchDir(String dir) {
        searchDirs.add(getCanonicalPath(dir, true, false, false, true));
    }

    /**
     * Add a list of sub-directories (relative to the resource base directory) that will be
     * searched to find requested resources.
     * @param dirs A list of subdirectories to add. If a directory starts with a slash (/)
     *            it indicates that it is an absolute directory instead of being
     *            relative to the resource base directory.
     */
    public void addSearchDirs(String... dirs) {
        for(String s : dirs)
            addSearchDir(s);
    }

    public List<String> getSearchDirs() { return searchDirs; }

    /**
     * Searches for the specified resource in all search directories and returns either
     * an InputStream (if the resource is found) or null (if it cannot be found).
     * @param resourceName The name of the resource to find.
     * @return An InputStream representing the resource if it is found, or null otherwise.
     */
    public InputStream tryGetStream(String resourceName) {
        String path = getResourceLocation(resourceName);
        if (path == null)
            return null;
        return loader.getResourceAsStream(path);
    }

    /**
     * Searches for the specified resource in all search directories and returns an
     * an InputStream representing the resource.
     * @param resourceName The name of the resource to find.
     * @return An InputStream representing the resource if it is found.
     * @throws IOException - The specified resource was not found or could not be loaded.
     */
    public InputStream getStream(String resourceName) throws IOException {
        return getStream(resourceName, false);
    }

    public InputStream getStream(String resourceName, boolean ignoreMissing)
            throws IOException {
        InputStream rs = tryGetStream(resourceName);
        if (rs == null) {
            if (ignoreMissing)
                return null;
            throw createNotFoundException(resourceName);
        }
        return rs;
    }

//    public boolean testResourceDir() {
//        return loader.getResource(baseDir) != null;
//    }
//
//    public boolean testResourceDir(String testDir) {
//        return loader.getResource(testDir) != null;
//    }

//    public void verifyResourceDir() throws IOException {
//        folderVerified = testResourceDir();
//        if (!folderVerified) {
//
//
//            StringBuilder sb = new StringBuilder();
//
//
//            URL[] urls = ((URLClassLoader)loader).getURLs();
//
//            String base = urls[0].toString();
//            String[] paths = new String[] {
//                    base,
//                    base + "/",
//                    base + "/resources",
//                    base + "/ur_rna",
//                    "/resources",
//                    "/ur_rna",
//                    "/",
//                    "resources",
//                    "ur_rna",
//                    "",
//                    "ur_rna/RNAstructureUI/resources/images/Draw.gif",
//                    "/ur_rna/RNAstructureUI/resources/images/Draw.gif",
//                    "RNAstructureUI/resources/images/Draw.gif",
//                    "/RNAstructureUI/resources/images/Draw.gif",
//                    "resources/images/Draw.gif",
//                    "/resources/images/Draw.gif",
//                    "ur_rna/RNAstructureUI/resources/images/*.gif",
//                    "ur_rna/RNAstructureUI/resources/images/*.",
//                    "ur_rna/RNAstructureUI/resources/images/*.*",
//                    "ur_rna/RNAstructureUI/resources/images/*",
//                    "ur_rna/RNAstructureUI/resources/images/.",
//                    "ur_rna/RNAstructureUI/resources/images/./",
//                    "ur_rna/RNAstructureUI/resources/images/./Draw.gif"
//            };
//            for (String s : paths) {
//                Object dir = loader.getResource(s);
//                if (dir == null)
//                    System.out.println("invalid resource dir: " + s);
//                else
//                    System.out.println("!! RESOURCE: " + s + " ==> " + ObjTools.toDisplayString(dir));
//            }
//
//            sb.append("Class paths searched for resources: \n");
//            for (URL u : urls)
//                sb.append('\t').append(u.toString()).append('\n');
//            AppLog.getDefault().warn(sb.toString());
//            throw new IOException("The application resources folder is missing. Expected location (relative to any path listed in the classPath): " + baseDir);
//        }
//    }

//    public InputStream getFirstStream(String... alternatePaths) throws IOException {
//        if (!folderVerified) verifyResourceDir();
//        for (String s : alternatePaths) {
//            InputStream ss = getStream(s, true);
//            if (ss != null) return ss;
//        }
//        throw createNotFoundException(String.join(" or ", alternatePaths));
//    }

    /**
     * Returns the full path to a resource if it is found, or null otherwise.
     */
    public String getResourceLocation(String resourceName) {
        String path = baseDir + resourceName;
        if (loader.getResource(path) != null) return path;
        for (String subDir : searchDirs) {
            if (subDir.startsWith("/"))
                path = subDir + resourceName;
            else
                path = baseDir + subDir + resourceName;
            if (loader.getResource(path) != null) return path;
        }
        return null; //not found
    }

    private static IOException createNotFoundException(String path) {
        return new IOException(String.format("Resource not found: %s.", path));
    }

    public List<String> getAllSearchDirs() {
        ArrayList<String> list = new ArrayList<>();
        URL[] urls = ((URLClassLoader)loader).getURLs();
        for (URL u : urls) {
            String base = u.toString();
            if (base.endsWith(".jar"))
                base += "!";
            else if (!(base.endsWith("/")||base.endsWith(".jar!")))
                base += "/";

            list.add(base + baseDir);
            for(String subDir : searchDirs) {
                if (subDir.startsWith("/"))
                    list.add(base + subDir);
                else
                    list.add(base + baseDir + subDir);
            }
        }
        return list;
    }
    /**
     * Determine whether the specified resource exists.
     * @param resourceName The resource to locate.
     * @return True if the resource could be found in one of the search directories or false otherwise.
     */
    public boolean hasResource(final String resourceName) {
        return getResourceLocation(resourceName) != null;
    }

//    private static String getSearchList(String path) {
//        StringBuilder sb = new StringBuilder();
//        for (String sub : _searchPaths) {
//            Path full = Paths.get(sub, path);
//            if (sb.length() != 0)
//                sb.append(", ");
//            sb.append(full.toString());
//        }
//        return sb.toString();
//    }

//    private static String findResource(String path) {
//        for (String sub : _searchPaths) {
//            String full = Paths.get(sub, path).toString();
//            if (loader.getResource(full) != null)
//                return full;
//        }
//        return null;
//    }

//    public InputStream getImage(String name, boolean ignoreMissing) throws IOException {
//        try {
//            return getFirstStream(imagesSubFolder + pathSep + name, name); //First try the images folder, then the name itself.
//        } catch (IOException ex) {
//            if (!ignoreMissing)
//                throw ex;
//        }
//        return null;
//    }

//    public static java.net.URL getResPath(String name) {
//        java.net.URL s = loader.getResource(prependPath + name);
//        return s;
//    }
}
