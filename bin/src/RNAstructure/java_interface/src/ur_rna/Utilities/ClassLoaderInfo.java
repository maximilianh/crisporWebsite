package ur_rna.Utilities;

import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.Vector;

/**
 * Gets information about native libraries that have been loaded by a {@link ClassLoader}.
 *
 * @author Richard M. Watson
 */
public class ClassLoaderInfo {
    private static final java.lang.reflect.Field reflectedClassLoaderField;
    //static Field[] fields;
    static {
        Field f = null;
        try {
            f = ClassLoader.class.getDeclaredField("loadedLibraryNames");
            f.setAccessible(true);
            //System.out.println(f.getName());
        } catch (NoSuchFieldException ex) {
            ex.printStackTrace();
        }
        reflectedClassLoaderField = f;
        //fields = ClassLoader.class.getDeclaredFields();
//        for (Field f : fields) {
//                  System.out.println("Field: " + f.getName() + " " + f.getModifiers());
//        }
    }
    @SuppressWarnings("unchecked") //calls to reflection functions (such as 'get') result in warnings about unsafe/unchecked code.
    public static String[] getLoadedLibraries(final ClassLoader loader) {
        try {
            final Vector<String> libraries = (Vector<String>) reflectedClassLoaderField.get(loader);
            return libraries.toArray(new String[libraries.size()]);
        } catch (Throwable t) {
            return null;
        }
    }
    public static String[] getLoadedLibraries(final ClassLoader[] loaders) {
        final Vector<String> libraries = new Vector<String>();
        for (ClassLoader l : loaders) {
            String[] list = getLoadedLibraries(l);
            if (list != null)
                libraries.addAll(Arrays.asList(list));
        }
        return libraries.toArray(new String[libraries.size()]);
    }
    public static String[] getAppLoadedLibraries() {
        return getLoadedLibraries(ClassLoader.getSystemClassLoader()); //MyClassName.class.getClassLoader()
    }
    public static String[] getMyLoadedLibraries() {
        return getLoadedLibraries(ClassLoaderInfo.class.getClassLoader());
    }
    public static String[] getAllLoadedLibraries() {
        ClassLoader[] list = new ClassLoader[2];
        list[0] = ClassLoader.getSystemClassLoader();
        list[1] = ClassLoaderInfo.class.getClassLoader();
        return getLoadedLibraries(list);
    }
}
