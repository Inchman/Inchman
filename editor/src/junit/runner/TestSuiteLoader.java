package junit.runner;

/**
 * An interface to define how a test suite should be loaded.
 */
public interface TestSuiteLoader {
	Class load(String suiteClassName) throws ClassNotFoundException;
	Class reload(Class aClass) throws ClassNotFoundException;
}
