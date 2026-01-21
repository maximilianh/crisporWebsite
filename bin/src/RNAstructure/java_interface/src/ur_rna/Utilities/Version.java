package ur_rna.Utilities;

/**
 * Implements semantic versioning for Java.
 * Versions consist of 4 parts, of which only the first two are required, separated by dots (.)
 * I.e.: Major.Minor.Revision.Build
 * Examples:
 *    3.2             (Major and Minor only)
 *    3.1.0.3124      (All components)
 */
public class Version implements Comparable<Version> {
    public final int Major, Minor, Revision, Build;
    public Version(final int major, final int minor, final int revision, final int build) {
        Major = major;
        Minor = minor;
        Revision = revision;
        Build = build;
    }
    public Version(final int major, final int minor, final int revision) {
        this(major, minor, revision, -1);
    }
    public Version(final int major, final int minor) {
        this(major, minor, -1, -1);
    }
    public Version(String version) throws IllegalArgumentException {
        int m=0, n=0, r=-1, b=-1;
        String[] parts = version.split("\\.");

        if (parts.length < 2) throw new IllegalArgumentException("A version number must have at least two components.");

        //if (parts.length > 0)
            m = Integer.parseInt(parts[0]);
        // if (parts.length > 1)
            n = Integer.parseInt(parts[1]);
        if (parts.length > 2)
            r = Integer.parseInt(parts[2]);
        if (parts.length > 3)
            b = Integer.parseInt(parts[3]);
        Major = m; Minor = n; Revision = r; Build = b;
    }
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(Major).append('.').append(Minor);
        if (Revision != -1) {
            sb.append('.').append(Revision);
            if (Build != -1)
                sb.append('.').append(Build);
        }
        return sb.toString();
    }

    @Override
    public int compareTo(final Version o) {
        if (Major != o.Major) return Major < o.Major ? -1 : 1;
        if (Minor != o.Minor) return Minor < o.Minor ? -1 : 1;
        if (Revision != o.Revision) return Revision < o.Revision ? -1 : 1;
        if (Build != o.Build) return Build < o.Build ? -1 : 1;
        return 0;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof Version)) return false;

        Version version = (Version) o;

        if (Major != version.Major) return false;
        if (Minor != version.Minor) return false;
        if (Revision != version.Revision) return false;
        return Build == version.Build;
    }

    @Override
    public int hashCode() {
        int result = Major;
        result = 31 * result + Minor;
        result = 31 * result + Revision;
        result = 31 * result + Build;
        return result;
    }

    public static boolean isVersionCompatible(final Version minVersion, final Version testVersion) {
        // The major version must match identically.
        if (testVersion.Major != minVersion.Major) return false;

        // The minor version of the test must be the same or greater than the reference.
        if (testVersion.Minor < minVersion.Minor) return false;

        // If the minor version of the test is greater, then don't compare the third component.
        // Otherwise, the third component of the test version (if provided) must be at least as high as the minimum version.
        return testVersion.Minor > minVersion.Minor || testVersion.Revision >= minVersion.Revision;
    }

    public static boolean isVersionCompatible(final String minVersion, final String testVersion) {
        String[] ref = minVersion.split("\\.");
        String[] test = testVersion.split("\\.");
        int len = Math.min(ref.length, test.length);

        if (len < 2) throw new IllegalArgumentException("A version number must have at least two components.");

        // The major version must match identically.
        if (Integer.parseInt(test[0]) != Integer.parseInt(ref[0])) return false;

        // The minor version of the test must be the same or greater than the reference
        if (Integer.parseInt(test[1]) < Integer.parseInt(ref[1])) return false;

        // If the minor version of the test is greater, then don't compare the third component.
        if (Integer.parseInt(test[1]) > Integer.parseInt(ref[1])) return true;

        // If there is no 3rd component, we get here if the third component is provided and the minor versions are the same
        return len == 2 || Integer.parseInt(test[2]) >= Integer.parseInt(ref[2]);
    }
    public boolean isEmpty() {
        return Major == 0 && Minor == 0 && Revision == -1 && Build == -1;
    }

    /** Returns a Version object if the passed-in String argument represents a valid version string. Otherwise,
     * returns null.
     * This function prints an error to STDERR (via Throwable.printStackTrace) if the version cannot be parsed.
     *
     * @param version A string representing a version number, e.g. 2.3.1
     *                A version string must have at least two parts (e.g. 6.0) but not more than four parts (e.g. 6.0.2.8888)
     *                The version string can be empty or null, in which case null is returned.
     * @return Null if the version string is invalid, empty, or null. Otherwise returns a version object.
     */
    public static Version tryParse(String version) {
        if (version == null || version.isEmpty()) return null;
        try {
            return new Version(version);
        } catch (IllegalArgumentException ex) {
            ex.printStackTrace();
        }
        return null;
    }

}
