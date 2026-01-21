package ur_rna.Utilities.annotation;

import ur_rna.Utilities.Version;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Indicates the Application name, title, and version when placed on the startup class (which contains main())
 */
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE) //can use in method only.
public @interface ApplicationInfo {
    String name();
    String title();
    String version();
}
