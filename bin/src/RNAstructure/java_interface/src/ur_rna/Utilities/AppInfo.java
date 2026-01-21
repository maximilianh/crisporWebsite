package ur_rna.Utilities;

import ur_rna.Utilities.annotation.ApplicationInfo;
import ur_rna.Utilities.annotation.NotNull;

import java.lang.reflect.Method;
import java.util.Collection;

/**
 * Provides information on the current running application.
 */
public class AppInfo {
    private static Class<?> mainClass;
    private static String appName;
    private static String appTitle;
    private static Version appVersion;
    private static AppInfoProvider provider;

    public static @NotNull String getAppTitle() {
        if (appTitle == null) appTitle = System.getProperty("APP_TITLE");
        if (appTitle == null) appTitle = System.getenv("JAVA_APP_TITLE");
        if (appTitle == null) appTitle = getProvider().getAppTitle();
        return appTitle;
    }
    public static @NotNull String getAppName() {
        if (appName == null) appName = System.getProperty("APP_NAME");
        if (appName == null) appName = System.getenv("JAVA_APP_NAME");
        if (appName == null) appName = getProvider().getAppName();
        return appName;
    }

    public static @NotNull Version getAppVersion() {
        if (appVersion == null) appVersion = Version.tryParse(System.getProperty("APP_VERSION"));
        if (appVersion == null) appVersion = Version.tryParse(System.getenv("JAVA_APP_VERSION"));
        if (appVersion == null) appVersion = getProvider().getAppVersion();
        return appVersion;
    }

    public static void setAppName(String name) { appName = name; }
    public static void setAppTitle(String title) { appTitle = title; }
    public static void setProvider(AppInfoProvider p) { provider = p; }
    public static void setMainClass(Class<?> c) { mainClass = c; }

    /**
     * Attempts to get the AppInfoProvider that has been registered for this application.
     * If one has not been registered, it attempts to get one by calling {@code getAppInfo() } on the
     * application's main class (if it exists).
     * Otherwise a default AppInfoProvider is created and returned. This default provides no valuable information,
     * so it is suggested that each application registers its own AppInfoProvider in main().
     */
    public static @NotNull AppInfoProvider getProvider() {
        if (provider != null) return provider;

        Class<?> main = getMainClass();
        if (main != null) {
            // Get an AppInfoProvider from the getAppInfo() function on the class with the main() startup function (if that function exists).
            try {
                Method m = main.getMethod("getAppInfo");
                return provider = (AppInfoProvider) m.invoke(null);
            } catch (Exception ex) {
                // failed.
            }
            // Get the application info from the @ApplicationInfo annotation attached to the class with the main() startup function.
            try {
                ApplicationInfo note = main.getAnnotation(ApplicationInfo.class);
                if (note != null) return new AnnotationAppInfoProvider(note);
            } catch (Exception ex) {
                // failed.
            }
        }

        // Create and return a default provider.
        provider = new AppInfoProvider() {
            private final Version ver = new Version("0.0");
            @Override public String getAppTitle() { return "JavaApplication"; }
            @Override public String getAppName() { return "JavaApplication"; }
            @Override public Version getAppVersion() { return ver; }
        };

        return provider;
    }

    public static Class<?> getMainClass() {
        if (mainClass != null)
            return mainClass;

        Collection<StackTraceElement[]> stacks = Thread.getAllStackTraces().values();
        for (StackTraceElement[] currStack : stacks) {
            if (currStack.length==0)
                continue;
            StackTraceElement lastElem = currStack[currStack.length - 1];
            if (lastElem.getMethodName().equals("main")) {
                try {
                    String mainClassName = lastElem.getClassName();
                    mainClass = Class.forName(mainClassName);
                    return mainClass;
                } catch (ClassNotFoundException e) {
                    // bad class name in line containing main?!
                    // shouldn't happen
                    e.printStackTrace();
                }
            }
        }
        return null;
    }

    public static class AnnotationAppInfoProvider  implements AppInfoProvider {
        public final ApplicationInfo info;
        public final Version ver;
        public AnnotationAppInfoProvider(ApplicationInfo info) { this.info = info; ver = new Version(info.version()); }
        @Override
        public String getAppTitle() {
            return info.title();
        }
        @Override
        public String getAppName() {
            return info.name();
        }
        @Override
        public Version getAppVersion() {
            return ver;
        }
    }
    public interface AppInfoProvider {
        @NotNull String getAppTitle();
        @NotNull String getAppName();
        @NotNull Version getAppVersion();
    }
}
