package junit.runner;

import java.awt.Component; 

import junit.framework.*;

/**
 * A view to show a details about a failure
 */
public interface FailureDetailView {
	/**
	 * Returns the component used to present the TraceView
	 */
	/*Component getComponent();*/
	/**
	 * Shows details of a TestFailure
	 */
	void showFailure(TestFailure failure);
	/**
	 * Clears the view
	 */
	void clear();
}
