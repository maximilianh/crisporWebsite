package ur_rna.Utilities.annotation;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Indicates that the decorated target represents an enumerated constant,
 * rather than a pure primitive value (e.g. int, char, String etc).
 * May be used by IDE or compiler to suggest valid values and/or warn of
 * potentially invalid values.
 */

@Retention(RetentionPolicy.SOURCE)
@Target({ElementType.FIELD, ElementType.PARAMETER, ElementType.LOCAL_VARIABLE, ElementType.ANNOTATION_TYPE, ElementType.METHOD})
public @interface MagicConstant {
    /**
     * Indicates that the allowed values for the target are those found in the specified
     * array of constant int or long values.
     */
    long[] intValues() default {};

    /**
     * Indicates that the allowed values for the target are those found in the specified
     * array of constant String values.
     */
    String[] stringValues() default {};

    /**
     * Indicates that the allowed values for the target are those either found in the specified
     * array (of constant int/long values) or a value produced by combining one or more
     * of them using the bitwise OR operator (|).
     */
    long[] flags() default {};

    /**
     * Indicates that the allowed values for the argument are the public static final fields
     * declared in the given class.
     */
    Class valuesFromClass() default void.class;

    /**
     * Indicates that the allowed values for the argument are the public static final fields
     * declared in the given class, or a value produced by combining one or more
     * of them using the bitwise OR operator (|).
     */
    Class flagsFromClass() default void.class;
}
