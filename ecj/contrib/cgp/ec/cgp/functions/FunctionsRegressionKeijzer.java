package ec.cgp.functions;

/**
 *
 * Function set for the Keijzer regression benchmarks.
 *
 * 
 * @author Roman Kalkreuth, roman.kalkreuth@tu-dortmund.de,
 *         https://orcid.org/0000-0003-1449-5131,
 *         https://ls11-www.cs.tu-dortmund.de/staff/kalkreuth,
 *         https://twitter.com/RomanKalkreuth
 *         
 *
 */
public class FunctionsRegressionKeijzer implements Functions {

	/** Add */
	static int F_ADD = 0;
	/** Multiply */
	static int F_MUL = 1;
	/** Inversion  */
	static int F_INV = 2;
	/** Negation */
	static int F_NEG = 3;
	/** Square root */
	static int F_SQRT = 4;

	/** Interpret the given function and apply it to the given inputs. */
	public Object callFunction(Object[] inputs, int function, int numFunctions) {
		if (function == F_ADD) {
			return (Double)inputs[0] + (Double)inputs[1];
		} else if (function == F_MUL) {
			return (Double)inputs[0] * (Double)inputs[1];
		} else if (function == F_SQRT) {
			return Math.sqrt((Double)inputs[0]);
		} else if (function == F_NEG) {
			return -1.0*((Double)inputs[0]);
		} else if (function == F_INV) {
			if ((Double)inputs[0] == 0) return (Double)inputs[0];
			else
				return 1.0/((Double)inputs[0]);
		}
		 else throw new IllegalArgumentException("Function #" + function + " is unknown.");
	}

	/**
	 * Return a function name, suitable for display in expressions, for the given
	 * function.
	 */
	public String functionName(int fn) {
		if (fn == F_ADD) return "+";
		if (fn == F_INV) return "INV";
		if (fn == F_MUL) return "*";
		if (fn == F_NEG) return "NEG";
		if (fn == F_SQRT) return "SQRT";
		else return "UNKNOWN FUNCTION";
	}

	/** Return the arity of the given function */
	public int arityOf(int fn) {
		if (fn == F_ADD) return 2;
		if (fn == F_INV) return 1;
		if (fn == F_MUL) return 2;
		if (fn == F_NEG) return 1;
		if (fn == F_SQRT) return 1;
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
