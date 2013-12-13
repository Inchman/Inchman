package junit.framework;

/**
 * A Listener for test progress
 */
public interface TestListener {
	/**
 	 * An error occurred.
 	 */
	void addError(Test test, Throwable t);
	/**
 	 * A failure occurred.
 	 */
 	void addFailure(Test test, AssertionFailedError t);  
	/**
	 * A test ended.
	 */
 	void endTest(Test test); 
	/**
	 * A test started.
	 */
	void startTest(Test test);
}
