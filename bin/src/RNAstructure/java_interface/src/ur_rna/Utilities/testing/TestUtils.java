package ur_rna.Utilities.testing;

import java.util.Objects;

/**
 * Utilites to aid in unit tests.
 */
public class TestUtils {
    public interface TestMethod {
        void run() throws Exception;
    }
    public static void assertError(TestMethod testLambda, Class expectedExceptionClass) {
        assertError(testLambda, expectedExceptionClass, true);
    }
    @SuppressWarnings("unchecked")
    public static void assertError(TestMethod testLambda, Class expectedExceptionClass, boolean allowDerivedExceptions) {
        Class actual = null;
        try {
            testLambda.run();
        } catch (Exception ex) {
            actual = ex.getClass();
        }
        boolean equal = Objects.equals(actual, expectedExceptionClass);
        if (!equal && allowDerivedExceptions && actual != null && expectedExceptionClass != null)
            equal = expectedExceptionClass.isAssignableFrom(actual);

        if (!equal) {
            String msg = String.format("The method did not throw the expected Exception.\nExpected: %s\nActual: %s",
                    expectedExceptionClass == null ? "Null" : expectedExceptionClass.getSimpleName(),
                    actual == null ? "Null" : actual.getSimpleName()
                    );
            throw new AssertionError(msg);
        }
    }

    public static void assertInstanceOf(Class expectedClass, Object actualObject) {
        assertInstanceOf(expectedClass, actualObject == null ? null : actualObject.getClass());
    }
    @SuppressWarnings("unchecked")
    public static void assertInstanceOf(Class expected, Class actual) {
        boolean equal = Objects.equals(actual, expected);
        if (!equal && actual != null && expected != null)
            equal = expected.isAssignableFrom(actual);
        if (equal) return;
        String msg = String.format("The actual object is not of the expected type.\nExpected: %s\nActual: %s",
                expected == null ? "Null" : expected.getSimpleName(),
                actual == null ? "Null" : actual.getSimpleName()
        );
        throw new AssertionError(msg);
    }
}
