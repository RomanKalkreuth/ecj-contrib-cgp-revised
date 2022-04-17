package ec.cgp.functions;

import ec.app.regression.func.VladERCA;
import ec.app.regression.func.VladERCB;
import ec.app.regression.func.VladERCC;

/**
 *
 * Function set A for the Vladislavleva regression benchmarks.
 *
 * 
 * @author Roman Kalkreuth, roman.kalkreuth@tu-dortmund.de,
 *         https://orcid.org/0000-0003-1449-5131,
 *         https://ls11-www.cs.tu-dortmund.de/staff/kalkreuth,
 *         https://twitter.com/RomanKalkreuth
 *         
 *
 */
public class FunctionsRegressionVladA implements Functions {

	static VladERCA vladERCA;
	static VladERCB vladERCB;
	static VladERCC vladERCC;
	
	/** Add */
	static int F_ADD = 0;
	/** Subtract */
	static int F_SUB = 1;
	/** Multiply */
	static int F_MUL = 2;
	/** Safe divide */
	static int F_DIV = 3;
	/** Square */
	static int F_N2 = 4;
	
	static int F_ERCA = 5;
	
    static int F_ERCB = 6;
    
    static int F_ERCC = 7;

	/** Interpret the given function and apply it to the given inputs. */
	public Object callFunction(Object[] inputs, int function, int numFunctions) {
		if (function == F_ADD) {
			return (Double)inputs[0] + (Double)inputs[1];
		} else if (function == F_SUB) {
			return (Double)inputs[0] - (Double)inputs[1];
		} else if (function == F_MUL) {
			return (Double)inputs[0] * (Double)inputs[1];
		} else if (function == F_DIV) {
			if ((Double)inputs[1] == 0) return (Double)inputs[1];
			return (Double)inputs[0] / (Double)inputs[1];
		} else if (function == F_N2) {
			double res = Math.pow((Double)inputs[0],2); 
			if(Double.isInfinite(res) || res >= Double.POSITIVE_INFINITY)
				return (Double)inputs[0];
			else
				return Math.exp((Double)inputs[0]); 
		}
		else if (function == F_ERCA) {
			return Math.pow((Double)inputs[0], vladERCA.value);
		} 
		else if (function == F_ERCB) {
			return (Double)inputs[0] + vladERCB.value;
		} 
		else if (function == F_ERCC) {
			return (Double)inputs[0] * vladERCC.value;
		}
		else throw new IllegalArgumentException("Function #" + function + " is unknown.");
	}

	/**
	 * Return a function name, suitable for display in expressions, for the given
	 * function.
	 */
	public String functionName(int fn) {
		if (fn == F_ADD) return "+";
		if (fn == F_SUB) return "-";
		if (fn == F_MUL) return "*";
		if (fn == F_DIV) return "/";
		if (fn == F_N2) return "SQUARE";
		if (fn == F_ERCA) return "ERCA";
		if (fn == F_ERCB) return "ERCB";
		if (fn == F_ERCC) return "ERCC";
		else return "UNKNOWN FUNCTION";
	}

	/** Return the arity of the given function */
	public int arityOf(int fn) {
		if (fn == F_ADD) return 2;
		if (fn == F_SUB) return 2;
		if (fn == F_MUL) return 2;
		if (fn == F_DIV) return 2;
		if (fn == F_N2) return 1;
		if (fn == F_ERCA) return 1;
		if (fn == F_ERCB) return 1;
		if (fn == F_ERCC) return 1;
		else throw new IllegalArgumentException("Function #" + fn+ " is unknown.");
	}

	/** Return the name, suitable for display, for the given input. */
	public String inputName(int inp, Object val) {
		if (inp == 0) return "x";
		if (inp == 1) return "y";
		if (inp == 2) return "z";
		if (inp == 3) return "w";
		if (inp == 4) return "v";
		else throw new IllegalArgumentException("Input #" + inp+ " is unknown.");
	}

}
