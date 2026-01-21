package ur_rna.Utilities;

import java.io.File;
import java.lang.reflect.Field;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Supplier;

/**
 * Loads a Native library or returns helpful error information if unsuccessful.
 */
public class JniLoader {
    private String libraryTitle;
    private String minVersion;
    private List<String> allowedNames = new ArrayList<>();
    private Supplier<String> versionCheck;

    public JniLoader(final String libraryTitle) { this(libraryTitle, null, null, (String[])null); }
    public JniLoader(final String libraryTitle, String... possibleLibNames) { this(libraryTitle, null, null, possibleLibNames); }
    public JniLoader(final String libraryTitle, final Supplier<String> versionCheck, final String requireMinVersion, String... possibleLibNames) {
        this.libraryTitle = libraryTitle;
        this.minVersion = requireMinVersion;
        if (possibleLibNames != null)
            Collections.addAll(allowedNames, possibleLibNames);
        this.versionCheck = versionCheck;
    }

    public List<String> getAllowedLibNames() {
        return allowedNames;
    }
    public void addLibName(String possibleName) {
        allowedNames.add(possibleName);
    }

    /**
     * Detect the architecture of the Java platform (32 or 64 bits
     * -- this is NOT necessarily the architecture of the OS!!)
     * @return returns result32 is the Java VM is 32-bit or result64 if the JVM is 64-bit.
     */
    public static String getJvmBits(String result32, String result64) {
        String jvmBits = System.getProperty("sun.arch.data.model");
        if (jvmBits.equals("64"))
            return result64;
        if (jvmBits.equals("32"))
            return result32;
        return "";
    }
    public static String getJvmBits() { return getJvmBits("_32","_64"); }

    public void load() throws NativeLibLoadError {
        String[] libNames = getLibNames();

        Throwable[] errors = new Throwable[libNames.length];
        String details[] = new String[libNames.length];

        if (tryLoadLibrary(libNames, errors)) return;

        Throwable loadError = null;
        // Inspect and/or modify the exception messages. Make note of any errors that are not simply due to the library not being found.
        for(int i = 0; i < errors.length; i++) {
            details[i] = errors[i] == null ? "Unknown Error" : errors[i].getMessage();
            if (details[i].equals("no " + libNames[i] + " in java.library.path"))
                details[i] = "Not found.";
            else if (loadError == null)
                loadError = errors[i];
        }
        StringBuilder sb = new StringBuilder();
        if (loadError == null) {
            // no real error occurred. Just couldn't find the library.
            sb.append("The required " + libraryTitle + " library is missing.\nIt should be in the same directory as the application, and it should be named:\n    ");
            for (int i = 0; i < libNames.length; i++) {
                if (i != 0)
                    sb.append("  or  ");
                sb.append(System.mapLibraryName(libNames[i]));
            }
            sb.append("\n\nFor troubleshooting purposes, note that the following locations were searched:\n  ")
                    .append(System.getProperty("java.library.path").replace(File.pathSeparator, "\n  "));
            sb.append("\n  (Current directory: ").append(Paths.get(".").toAbsolutePath().normalize().toString()).append(")\n");
        } else {
            sb.append("Unable to load the " + libraryTitle + " library!\n")
                    .append(loadError.getMessage())
                    .append("\n\n-------- Error Details (For all attempted paths): --------\n");

            for(int i = 0; i < errors.length; i++) {
                sb.append(i + 1).append(")  Attempted to load ").append(System.mapLibraryName(libNames[i]))
                        .append("\n      Result: ").append(details[i]).append("\n");
            }
            sb.append("\nAttempted Paths:\n").append(System.getProperty("java.library.path").replace(File.pathSeparator, ", "));
        }

        throw new NativeLibLoadError(sb.toString());
    }

    public String getLoadedLib() {
        String[] loaded = ClassLoaderInfo.getMyLoadedLibraries();
        for (String name : getLibNames()) {
            String find = name.toLowerCase();
            for (String s : loaded) {
                if (s.toLowerCase().contains(find))
                    return s;
            }
        }
        return null;
    }

    private boolean tryLoadLibrary(String[] libNames, Throwable errors[]) {
        // Add the program's directory to the PATH. (i.e. the directory containing the JAR file or the
        setLibPath(System.getProperty("java.library.path") + File.pathSeparator + PathTools.findLocalPath(".")+"/");

        for (int i = 0; i < libNames.length; i++) {
            try {
                System.loadLibrary(libNames[i]);
                // Seems successful so far. Try to call a known method to verify.
                if (versionCheck != null)
                    try {
                        String version = versionCheck.get();
                        if (!Version.isVersionCompatible(minVersion, version))
                            throw new UnsatisfiedLinkError(String.format("The library version is %s, but a version compatible with %s is required.", version, minVersion));
                    } catch (UnsatisfiedLinkError ex) {
                        throw new UnsatisfiedLinkError("The " + libNames[i] + " library is corrupt or has the wrong version or binary format.");
                    }
                return true;
            } catch (Throwable ex) {
                // Error is usually a UnsatisfiedLinkError, but let's catch all of them so we can show them to the user.
                errors[i] = ex;
            }
        }
        return false;
    }
    /**
     * Set the java.library.path System property and ensure that it will be re-evaluated by ClassLoaders
     * @param javaLibraryPath The new setting for java.library.path
     */
    private void setLibPath(String javaLibraryPath) {
        System.setProperty("java.library.path", javaLibraryPath);
        try {
            // Set the private static ClassLoader.sys_paths field to null, so that java.library.path is re-evaluated when LoadLibrary is called.
            Field fieldSysPath = ClassLoader.class.getDeclaredField("sys_paths");
            fieldSysPath.setAccessible(true);
            fieldSysPath.set(null, null);
        } catch (Exception ex) {
            //
        }
    }
    private String[] getLibNames() {
        if (allowedNames.size() == 0)
            return new String[] { libraryTitle, libraryTitle + getJvmBits() };
        else
            return allowedNames.toArray(new String[allowedNames.size()]);
    }

    public static class NativeLibLoadError extends Exception {
        private static final long serialVersionUID = 1L;
        public NativeLibLoadError() { }
        public NativeLibLoadError(String message) {            super(message);        }
        public NativeLibLoadError(Throwable cause) {            super(cause);        }
        public NativeLibLoadError(String message, Throwable cause) {            super(message, cause);        }
    }
}
