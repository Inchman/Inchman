package junit.framework;

/**
 * A <em>Test</em> can be run and collect its results.
 *
 * @see TestResult
 */
public interface Test {
	/**
	 * Counts the number of test cases that will be run by this test.
	 */
	int countTestCases();
	/**
	 * Runs a test and collects its result in a TestResult instance.
	 */
	void run(TestResult result);
}
