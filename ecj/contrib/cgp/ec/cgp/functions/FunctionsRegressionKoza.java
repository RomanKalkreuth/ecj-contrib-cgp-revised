package ec.cgp.functions;

/**
 *
 * Function set for the Koza regression benchmarks.
 * 
 * @author David Oranchak, doranchak@gmail.com, http://oranchak.com
 * 
 * @author Roman Kalkreuth, roman.kalkreuth@tu-dortmund.de,
 *         https://orcid.org/0000-0003-1449-5131,
 *         https://ls11-www.cs.tu-dortmund.de/staff/kalkreuth,
 *         https://twitter.com/RomanKalkreuth
 *         
 *
 */
public class FunctionsRegressionKoza implements Functions {

	/** Add */
	static int F_ADD = 0;
	/** Subtract */
	static int F_SUB = 1;
	/** Multiply */
	static int F_MUL = 2;
	/** Safe divide */
	static int F_DIV = 3;
	/** Sinus */
	static int F_SIN = 4;
	/** Cosinus */
	static int F_COS = 5;
	/** Power*/
	static int F_POW = 6;
	/** Logarithm */
	static int F_LOG = 7;

	/** Interpret the given function and apply it to the given inputs. */
	public Object callFunction(Object[] inputs, int function, int numFunctions) {
		if (function == F_ADD) {
			return (Double) inputs[0] + (Double) inputs[1];
		} else if (function == F_SUB) {
			return (Double) inputs[0] - (Double) inputs[1];
		} else if (function == F_MUL) {
			return (Double) inputs[0] * (Double) inputs[1];
		} else if (function == F_DIV) {
			if ((Double) inputs[1] == 0)
				return 1.0;
			return (Double) inputs[0] / (Double) inputs[1];
		} else if (function == F_SIN) {
			return Math.sin((Double) inputs[0]);
		} else if (function == F_COS) {
			return Math.cos((Double) inputs[0]);
		} else if (function == F_POW) {

			double res = Math.exp((Double) inputs[0]);
			if (Double.isInfinite(res) || res >= Double.POSITIVE_INFINITY) {
				return (Double) inputs[0];
			} else {
				return Math.exp((Double) inputs[0]);
			}	
		} else if (function == F_LOG) {
			
			double res = Math.log(Math.abs((Double) inputs[0]));
			if (Double.isInfinite(res) || res >= Double.POSITIVE_INFINITY) {
				return (Double) inputs[0];
			}
			else {
				return res;
			}
		}

		else
			throw new IllegalArgumentException("Function #" + function + " is unknown.");
	}

	/**
	 * Return a function name, suitable for display in expressions, for the given
	 * function.
	 */
	public String functionName(int fn) {
		if (fn == F_ADD)
			return "+";
		if (fn == F_SUB)
			return "-";
		if (fn == F_MUL)
			return "*";
		if (fn == F_DIV)
			return "/";
		if (fn == F_SIN)
			return "SIN";
		if (fn == F_COS)
			return "COS";
		if (fn == F_POW)
			return "POW";
		if (fn == F_LOG)
			return "LOG";
		else
			return "UNKNOWN FUNCTION";
	}

	/** Return the arity of the given function */
	public int arityOf(int fn) {
		if (fn == F_ADD) return 2;
		if (fn == F_SUB) return 2;
		if (fn == F_MUL) return 2;
		if (fn == F_DIV) return 2;
		if (fn == F_SIN) return 1;
		if (fn == F_COS) return 1;
		if (fn == F_POW) return 1;
		if (fn == F_LOG) return 1;
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
