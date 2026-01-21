package ur_rna.Utilities;

import ur_rna.Utilities.annotation.NotNull;
import ur_rna.Utilities.annotation.Nullable;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.file.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

public final class PathTools {
    private PathTools() {} //cannot instantiate.

    /**
     * Gets the base name of the file -- i.e. the name without extension.
     * If the file name contains multiple dot (.) characters, the extension is defined as the segment starting with
     * the <em>final</em> dot.
     *
     * @param fullName The full file name or path.
     * @return The base name of the file -- i.e. the name without extension.
     * @see #getBaseName(String, boolean)
     */
    public static String getBaseName(String fullName) { return getBaseName(fullName, false); }

    /**
     * Gets the base name of the file -- i.e. the name without extension.
     *
     * @param fullName                The full file name or path.
     * @param allowMultiDotExtensions If true, the extension starts at the <em>first</em> dot (.) encountered, which
     *                                means the resulting base name will never contain any dots.
     *                                If false, the extension starts at the <em>last</em> dot, so the base name will
     *                                contain all dots except the final one.
     * @return The base name of the file -- i.e. the name without any extension.
     */
    public static String getBaseName(String fullName, boolean allowMultiDotExtensions) { return parse(fullName, allowMultiDotExtensions).baseName(); }

    /**
     * Gets the extension of the file <em>including</em> the preceding dot (.)
     * If the file name contains multiple dot (.) characters, the extension is defined as the segment starting with
     * the <em>final</em> dot.
     *
     * @param fullName The full file name or path.
     * @return The extension of the file <em>including</em> the preceding dot (.)
     */
    public static String getExt(String fullName) { return getExt(fullName, true, false); }

    /**
     * Gets the extension of the file <em>including</em> the preceding dot (.)
     *
     * @param fullName                The full file name or path.
     * @param includeDot
     *@param allowMultiDotExtensions If true, the extension starts at the <em>first</em> dot (.) encountered in
     *                                the file name, which means the resulting extension may contain more than one dot.
     *                                If false, the extension starts at the <em>last</em> dot, so the extension will
     *                                either be empty (if the file has no extension) or it will contain exactly one
     *                                dot as its first character.  @return The extension of the file <em>including</em> the preceding dot (.)
     */
    public static String getExt(String fullName, final boolean includeDot, boolean allowMultiDotExtensions) { return parse(fullName, allowMultiDotExtensions).ext(includeDot); }

    /**
     * Get the full parent directory path, including the final slash.
     *  <table>
     *  <caption>Examples:</caption>
     *  <tr><td>  /path/to/fileName.ext </td><td>{@code -->}</td><td>  /path/to/     </td></tr>
     *  <tr><td>  /fileName.ext         </td><td>{@code -->}</td><td>  /             </td></tr>
     *  <tr><td>  C:\Windows\System32   </td><td>{@code -->}</td><td>  C:\Windows\   </td></tr>
     *  <tr><td>  images\banana.jpg     </td><td>{@code -->}</td><td>  images\       </td></tr>
     *  <tr><td>  HelloWorld.txt        </td><td>{@code -->}</td><td>  (empty)       </td></tr>
     *  </table>
     * @param fullName The full path to the file.
     * @return The full parent directory path, including the final slash, or an empty string if the full
     * path does not contain a directory separator character {@code \} or {@code /}.
     */
    public static String getDir(String fullName) { return parse(fullName).dir(); }
    public static String getDir(String fullName, boolean includeFinalSlash) { return parse(fullName).dir(includeFinalSlash); }

    public static FileName parse(String fullFilePath) { return parse(fullFilePath, false); }
    public static FileName parse(String fullFilePath, boolean allowMultiDotExtensions) { return new FileName(fullFilePath, allowMultiDotExtensions); }

    /**
     * Gets the absolute path to the directory containing the class file for referenceClass.
     * For a jar'd application, this is the directory containing the Jar file in which referenceClass resides.
     * @return The absolute path to the directory containing the referenceClass. The function returns {@code null}
     * if security restrictions do not permit resolution of paths.
     */
    public static @Nullable String getAppPath(Class referenceClass) {
        try {
            URL jarPath = referenceClass.getProtectionDomain().getCodeSource().getLocation();
            File jarFile = new File(jarPath.toURI());
            return jarFile.getParent();
        } catch (Exception ex) {
            return null;
        }
    }
    public static File getHomeDir() {
        File f = fileFromPath(System.getProperty("user.home"), true);
        if (f == null) f = fileFromPath(System.getProperty("user.dir"), false);
        return f;
    }
    private final static String[] documentsDirs = "Documents documents".split(" ");
    public static File getDocumentsDir() {
        File f = getHomeDir();
        try {
            if (!f.exists()) return f;
            Path p = f.toPath();
            File tmp;
            for (String s : documentsDirs)
              if ((tmp = p.resolve(s).toFile()).isDirectory()) return tmp;
        } catch (Exception ex) {
            // Path resolution failed for some reason.
            ex.printStackTrace();
        }
        return f;
    }
    public static File fileFromPath(String path) { return fileFromPath(path, false); }
    public static File fileFromPath(String path, boolean mustExist) {
        if (path == null) return null;
        File f = new File(path);
        return mustExist && !f.exists() ? null : f;
    }
    /**
     * Verify that a directory is writeable by writing a new file (and then deleting it).
     * If the file represents an existing file, the return value will be file.canWrite()
     * If the file does not exist, the return value is false.
     * @param file the file or directory to test for writing
     */
    public static boolean verifyWritable(File file) {
        if (file == null)
            throw new NullPointerException("In verifyWritable: file cannot be null.");
        if (file.isDirectory()) {
            try {
                Path p = Files.createTempFile(file.toPath(), "~", ".tmp");
                Files.delete(p);
                return true;
            } catch (Exception ex) {
                return false;
            }
        }
        return file.canWrite();
    }
    public static String changeExtension(final String file, final String extension) {
        return new FileName(file).changeExt(extension);
    }

    public static boolean exists(String fileOrDirectory) {
        try {
            if (fileOrDirectory != null) return new File(fileOrDirectory).exists();
        } catch (Exception ex) {
            // Possible RuntimeException, e.g. AccessControlException
        }
        return false;
    }
    public static boolean isFile(String file) {
        try { if (file != null) return new File(file).isFile(); }
        catch (Exception ex) { /* Possible RuntimeException, e.g. AccessControlException */ }
        return false;
    }
    public static boolean isDir(String directory) {
        try { if (directory != null) return new File(directory).isDirectory(); }
        catch (Exception ex) { /* Possible RuntimeException, e.g. AccessControlException */ }
        return false;
    }

    @Nullable
    public static File getFirstExisting(@Nullable File parentDir, boolean includeDirs, boolean includeFiles, String ... paths) {
        if (!includeDirs && !includeFiles) return null;

        for (String path : paths) {
            File test = new File(parentDir, path);
            if (test.exists()) {
                if (includeDirs && includeFiles)
                    return test;
                // return this file if:
                //   includeDirs==true and test.isDirectory()==true
                //     OR
                //   includeFiles==true and test.isFile()==true
                //     ..at this point this is the same as:
                //   includeDirs=false and test.isDirectory() == false
                //  i.e. return test if they are both true or both false.
                if (includeDirs == test.isDirectory())
                    return test;
            }
        }
        return null; //not found
    }

    public static class FileName {
        private static String[] EMPTY_ARRAY = new String[0];
        //private static String[] SINGLE_EMPTY_ITEM_ARRAY = new String[]{""};
        public final String fullPath;
        private final int sepPos, extPos;

        public boolean hasDir() { return sepPos!=-1;}
        public boolean hasExt() { return extPos!=-1;}

        /**
         * @return the file extension, preserving the initial dot.
         */
        public String ext() {
            return extPos==-1?"":fullPath.substring(extPos);
        }
        public String ext(boolean includeDot) {
            return extPos==-1?"":fullPath.substring(extPos+(includeDot?0:1));
        }

        /**
         * Returns the directory portion of the path, INCLUDING the final slash (!!)
         * (The final slash may be necessary to distinguish eg.  "/file.txt" from "file.txt" which would
         * both have a directory of "" if the slash were discarded.)
         */
        public String dir() {
            //    /path/to/fileName.ext -->     path=/path/to/      name=fileName.ext
            //    /fileName.ext         -->     path=/              name=fileName.ext
            return sepPos==-1?"":fullPath.substring(0, sepPos + 1); //keep ending slash if present, to distinguish "" from "/"
        }
        public String dir(boolean includeFinalSlash) {
            //    /path/to/fileName.ext -->     path=/path/to/      name=fileName.ext
            //    /fileName.ext         -->     path=/              name=fileName.ext
            return sepPos==-1?"":fullPath.substring(0, sepPos+(includeFinalSlash?1:0)); //keep ending slash if present, to distinguish "" from "/"
        }
        public String name() {
            return sepPos==-1?fullPath:fullPath.substring(sepPos+1);
        }
        public String baseName() {
            return extPos==-1?name():fullPath.substring(sepPos+1,extPos);
        }


        public FileName(String fullPath) { this(fullPath, false); }
        public FileName(String fullPath, boolean allowMultiDotExtensions) {
            this.fullPath = fullPath;
            sepPos = Math.max(fullPath.lastIndexOf('/'), fullPath.lastIndexOf('\\'));

            // allow for dots in the directory e.g. /conf.d/hello.txt
            // if allowMultiDotExtensions, use first dot.
            // otherwise, use last dot.
            //  /conf.d/hello.txt.bak
            //       |       |   |
            //     ignore    |   |
            //       multi-dot   standard
            if (allowMultiDotExtensions)
                extPos = fullPath.indexOf('.', sepPos+1);
            else if (sepPos==-1)
                extPos = fullPath.lastIndexOf('.');
            else {
                int pos = fullPath.substring(sepPos + 1).lastIndexOf('.');
                extPos = pos==-1?-1:pos+sepPos+1;
            }

            //Note: the index argument to lastIndexOf is the starting position for a reverse search,
            // NOT the end position. So it cannot be used to prevent the search from entering the directory.
            // So we must use substring to prevent matches in the directory part.
        }
        /**
         * Return an array of all ancestor directories of the file.
         * If the path begins with a slash, an empty string is the first element in the array.
         *
         * Examples:
         *      "/path/to/the/file.txt"    would yield  { "", "path", "to", "the }
         *      "relpath/to/the/file.txt"  would yield  { "relpath", "to", "the" }
         *      "file.txt"  would yield  { } (an empty array)
         *      "/file.txt" would yield  { "" } (an array containing a single empty string)
         */
        public String[] splitDirs() {
            //   "/path/to/"  --> [ "", "path", "to" ]
            //   "/"          --> [ "" ]
            if (sepPos==-1)
                return EMPTY_ARRAY;
            if (sepPos==0)
                return new String[]{""}; // e.g.  /file.txt would return [""]
            return reSlash.split(fullPath.substring(0, sepPos));
        }
        private final static Pattern reSlash = Pattern.compile("[/\\\\]");
        /**
         * Returns the original path, including the directory, but without any extension.
         * e.g. for  "C:\Users\John\Docs\hello.txt" this returns "C:\Users\John\Docs\hello"
         * The final slash will
         * @return
         */
        public String removeExt() {
            return extPos==-1?fullPath:fullPath.substring(0,extPos);
        }
        public String changeExt(final String extension) {
            if (extension == null || extension.isEmpty())
                return removeExt();
            if (extension.charAt(0)=='.')
                return removeExt() + extension;
            return removeExt() + "." + extension;
        }
    }

    public static String addSlash(String path) {
        // add a slash if it is not already there.
        if ((path.endsWith("/") || path.endsWith("\\") || path.endsWith(Character.toString(File.separatorChar))))
            return path;
        return path + File.separatorChar;
    }

    /**
     * Get the canonical form of a file-system path by:
     * <ol>
     *     <li>Converting it to unix format (if Windows)</li>
     *     <li>Replacing {@code //} with {@code /}</li>
     *     <li>Replacing {@code /./} with {@code /}</li>
     *     <li>Replacing {@code dirname/../} with {@code /}</li>
     * </ol>
     * @param path The path to canonicalize
     * @return A canonical unix-style form of the path
     * @see #getCanonicalPath(String, boolean, boolean, boolean, boolean)
     */
    public static String getCanonicalPath(String path) {
        return getCanonicalPath(path, false, false, false, false);
    }

    /**
     * Get the canonical form of a file-system path by:
     * <ol>
     *     <li>Converting it to unix format (if Windows)</li>
     *     <li>Removing starting slashes (if {@code removeStartSlash} is true)</li>
     *     <li>Adding an ending slash (if {@code addEndSlash} is true}</li>
     *     <li>Removing an ending slash (if {@code removeEndSlash} is true</li>
     *     <li>Replacing {@code //} with {@code /}</li>
     *     <li>Replacing {@code /./} with {@code /}</li>
     *     <li>Replacing {@code dirname/../} with {@code /}</li>
     *     <li>Removing {@code ./} at the start of a path (if {@code removeDotSlashAtStart} is true)</li>
     * </ol>
     * @param path The path to canonicalize
     * @param addEndSlash If true, any slashes at the end of the path will be removed. i.e.  {@code path/}  becomes  {@code path}
     * @param removeEndSlash If true any path not ending with a slash will have one appended. i.e.  {@code path}  becomes  {@code path/}
     * @param removeStartSlash If true, any slashes at the beginning of the path will be removed. i.e.  {@code /path}  becomes  {@code path}
     * @param removeDotSlashAtStart If this is true and the path starts with {@code ./}, it will be removed.  i.e.  {@code ./path}  becomes  {@code path}
     * @return A canonical unix-style form of the path
     */
    public static String getCanonicalPath(String path, boolean addEndSlash, boolean removeEndSlash, boolean removeStartSlash, boolean removeDotSlashAtStart) {
        if (path == null) return null;
        if (path.length() == 0) return path;

        path = path.replace("\\", "/");

        // TODO: Is "/" considered an starting slash or ending slash? I.e. should it be removed if removeStartSlash is true? Should it be removed if removeEndSlash is true?

        if (removeStartSlash) {
            // Remove starting slash
            while (path.startsWith("/"))
                path = path.substring(1);
        }

        if (addEndSlash && !path.endsWith("/"))
            path += "/";

        if (removeEndSlash && path.endsWith("/") && path.length() > 1)
            path = path.substring(0, path.length()-1);

        String oldPath;
        Pattern updir = Pattern.compile("([^/]?[^/.][^/])?/\\.\\./");  // replace dirname/../ with just /
        Pattern updirEnd = Pattern.compile("([^/]?[^/.][^/])?/\\.\\.$");  // replace dirname/.. (at the end of a path) with just /
        do {
            oldPath = path;
            path = path.replace("//", "/");
            path = path.replace("/./", "/");
            path = updir.matcher(path).replaceAll("/");
            path = updirEnd.matcher(path).replaceAll("");
        } while (!oldPath.equals(path));

        if (removeDotSlashAtStart && path.startsWith("./"))
            path = path.substring(2);

        return path;
    }

    /**
     * Returns the true if any files in the given directory match the specified filter, which should be a
     *  a glob expression such as "*.{java,class,jar}"
     * Glob expressions are fairly uniform across platforms. See {@link PathMatcher} for glob details.
     * @param dirPath The path of the directory to search.
     * @param filter The glob expression to match files against.
     * @return True if at least one file in the directory matches the filter, or false otherwise.
     */
    public static boolean anyFilesMatch(String dirPath, @NotNull String filter) {
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(Paths.get(dirPath), filter)) {
            Iterator<Path> iterator = stream.iterator();
            return iterator.hasNext();
        } catch (IOException ex) {
            throw new RuntimeException(String.format("Error reading folder %s: %s", dirPath, ex.getMessage()), ex);
        }
    }

    /**
     * Lists all files in a directory that match the specified filter, which should be a
     *  a glob expression such as "*.{java,class,jar}"
     *  Glob expressions are fairly uniform across platforms. See {@link PathMatcher} for glob details.
     * @param dirPath The path of the directory to search.
     * @param filter The glob expression to match files against.
     * @return A list of Files that match the filter.
     */
    public static List<File> listFiles(String dirPath, @NotNull String filter) {
        Path dir = Paths.get(dirPath);
        List<File> files = new ArrayList<>();
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(dir, filter)) {
            for (Path entry : stream) {
                files.add(entry.toFile());
            }
            return files;
        } catch (IOException ex) {
            throw new RuntimeException(String.format("Error reading folder %s: %s", dir, ex.getMessage()), ex);
        }
    }

    /**
     *    Unlike the other getCanonicalPath overloads, this one simply uses {@link java.io.File#getCanonicalPath()}
     */
    public static String getCanonicalPath(final File file) {
        try {
            return file.getCanonicalPath();
        } catch (IOException ex) {
            return file.getPath();
        }
    }

    private static String[] defaultLocalSearchSubDirs = { ".", "..", "../.." };
    public static File findLocalPath(String name) {return findLocalPath(PathTools.class, name, defaultLocalSearchSubDirs);  }
    public static File findLocalPath(Class refType, String name, String... trySubDirs) {
        String appPath = PathTools.getAppPath(refType); // path to directory containing jar file
        if (appPath == null) return null;
        for (String s : trySubDirs) {
            Path p = Paths.get(appPath, s, name);
            File f = p.toFile();
            if (f.exists()) {
                try {
                    f = f.getCanonicalFile();
                } catch (IOException ex) {
                    // keep file as-is
                }
                return f;
            }
        }
        return null;
    }
}
