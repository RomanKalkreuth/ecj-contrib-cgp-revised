package ec.cgp.problems.regression;

import ec.*;
import ec.app.regression.Benchmarks;
import ec.app.regression.func.KeijzerERC;
import ec.app.regression.func.KornsERC;
import ec.app.regression.func.VladERCA;
import ec.app.regression.func.VladERCB;
import ec.app.regression.func.VladERCC;
import ec.simple.*;
import ec.util.*;
import ec.vector.*;
import ec.cgp.Evaluator;
import ec.cgp.FitnessCGP;
import ec.cgp.problems.ProblemCGP;
import ec.cgp.representation.VectorIndividualCGP;
import ec.cgp.representation.VectorSpeciesCGP;
import ec.multiobjective.*;

import java.util.*;

/**
 * 
 * Regression problem.
 * 
 * @author David Oranchak, doranchak@gmail.com, http://oranchak.com
 *
 */
public class ProblemRegression extends ProblemCGP {

	/** Which function to use. Acceptable values: {1, 2, 3}. */
	public int function;

	static String P_WHICH = "which";

	final double PROBABLY_ZERO = 1.11E-15;
	final double BIG_NUMBER = 1.0e15;

	KeijzerERC keijzerERC;
	KornsERC kornsERC;

	double[][] pagieTraining1;
	double[][] vladTraining4;
	double[][] kornsTraining12;

	/** 50 randomly generated test points used for fitness evaluation */
	static double[] kozaTraining = new double[] { -0.48097666, -0.58232966, -0.076965828, 0.91879868, -0.23519547,
			0.67922012, -0.13699464, 0.01777318, 0.17666018, 0.45179073, 0.38780204, -0.65791783, -0.14948722,
			0.8038485, -0.11198381, 0.66869933, 0.50367412, -0.45215628, 0.42066769, -0.76327382 };

	static double[] nguyenTraining7 = new double[] { 0.23731074, 1.0476574, 0.056487079, 0.70340092, 0.74965303,
			1.738864, 1.4377288, 0.49024159, 0.17674821, 0.71894968, 1.0775117, 1.7473566, 0.73197111, 0.85822845,
			1.5754444, 1.333269, 0.94881026, 0.51224548, 0.89443731, 0.94453242 };

	static double[] nguyenTraining8 = new double[] { 2.9808034, 0.99121923, 0.78747278, 3.285668, 3.3622594, 2.74515,
			0.95739789, 3.3073568, 2.867567, 3.6054944, 1.2326744, 2.3324854, 0.89304816, 2.6587731, 3.4867153,
			3.322763, 1.1421899, 1.755899, 1.2972773, 0.69565992 };

	static double[] nguyenTraining9x = new double[] { 0.93434123, 0.34562961, 0.6609156, 0.94917166, 0.52490733,
			0.45056424, 0.95421317, 0.66350838, 0.75686555, 0.75020616, 0.14699992, 0.059720233, 0.21891122, 0.37390814,
			0.76905837, 0.3997893, 0.76472147, 0.91974642, 0.067912085, 0.64615459 };

	static double[] nguyenTraining9y = new double[] { 0.86087497, 0.5145302, 0.65529209, 0.094423176, 0.68331143,
			0.7875154, 0.51891707, 0.66923185, 0.59261582, 0.39182128, 0.098515389, 0.57886418, 0.85492571, 0.2987171,
			0.7247056, 0.23660522, 0.042563169, 0.11280032, 0.73548888, 0.59618834 };

	static double[] keijzerTraining6 = new double[] { 1.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0,
			14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 31.0, 32.0,
			33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0,
			50.0 };

	/** Evaluate the CGP and compute fitness */
	public void evaluate(EvolutionState state, Individual ind, int subpopulation, int threadnum) {
		if (ind.evaluated)
			return;

		VectorSpeciesCGP s = (VectorSpeciesCGP) ind.species;
		VectorIndividualCGP ind2 = (VectorIndividualCGP) ind;

		double diff = 0.0;
		double error = 0.0;

		Double[] inputs = new Double[s.numInputs];
		double[] x = null;
		double[] y = null;
		double[] z = null;
		double[] v = null;
		double[] w = null;

		double fn = 0.0;

		if (function >= 1 && function <= 12) {

			if (function == 7) {
				x = nguyenTraining7;
			} else if (function == 8) {
				x = nguyenTraining8;
			} else if (function >= 9 && function <= 12) {
				x = nguyenTraining9x;
				y = nguyenTraining9y;
			} else {
				x = kozaTraining;
			}

			for (int i = 0; i < x.length; i++) {

				inputs[0] = x[i]; // one of the randomly-generated independent variables.

				if (function >= 9 && function <= 12) {
					inputs[1] = y[i];
				} else {
					inputs[1] = 1.0;
				}

				/* run the CGP */
				Object[] outputs = Evaluator.evaluate(state, threadnum, inputs, ind2);

				if (function == 1)
					fn = koza1(x[i]);
				else if (function == 2)
					fn = koza2(x[i]);
				else if (function == 3)
					fn = koza3(x[i]);
				else if (function == 4)
					fn = nguyen4(x[i]);
				else if (function == 5)
					fn = nguyen5(x[i]);
				else if (function == 6)
					fn = nguyen6(x[i]);
				else if (function == 7)
					fn = nguyen7(x[i]);
				else if (function == 8)
					fn = nguyen8(x[i]);
				else if (function == 9)
					fn = nguyen9(x[i], y[i]);
				else if (function == 10)
					fn = nguyen10(x[i], y[i]);
				else if (function == 11)
					fn = nguyen11(x[i], y[i]);
				else if (function == 12)
					fn = nguyen12(x[i], y[i]);
			

				error = Math.abs((Double) outputs[0] - fn);

				if (!(error < BIG_NUMBER))
					error = BIG_NUMBER;
				else if (error < PROBABLY_ZERO)
					error = 0.0;

				diff += error;
				
			}
		} else if (function == 13) {

			for (int i = 0; i < pagieTraining1.length; i++) {
				inputs[0] = pagieTraining1[i][0];
				inputs[1] = pagieTraining1[i][1];

				Object[] outputs = Evaluator.evaluate(state, threadnum, inputs, ind2);
				
				fn = pagie1(pagieTraining1[i][0], pagieTraining1[i][1]);
				error = Math.abs((Double) outputs[0] - fn);

				if (!(error < BIG_NUMBER))
					error = BIG_NUMBER;
				else if (error < PROBABLY_ZERO)
					error = 0.0;

				diff += error;
			}

		} else if (function == 14) {

			for (int i = 0; i < keijzerTraining6.length; i++) {

				inputs[0] = keijzerTraining6[i];
				inputs[1] = keijzerERC.value;

				Object[] outputs = Evaluator.evaluate(state, threadnum, inputs, ind2);

				fn = keijzer6(keijzerTraining6[i]);
				error = Math.abs((Double) outputs[0] - fn);

				if (!(error < BIG_NUMBER))
					error = BIG_NUMBER;
				else if (error < PROBABLY_ZERO)
					error = 0.0;

				diff += error;
			}
		}
		else if (function == 15) {

			for (int i = 0; i < kornsTraining12.length; i++) {
				
				inputs[0] = kornsTraining12[i][0];
				inputs[1] = kornsTraining12[i][1];
				inputs[2] = kornsTraining12[i][2];
				inputs[3] = kornsTraining12[i][3];
				inputs[4] = kornsTraining12[i][4];
				inputs[5] = kornsERC.value;
				
				Object[] outputs = Evaluator.evaluate(state, threadnum, inputs, ind2);
				fn = korns12(kornsTraining12[i][0], kornsTraining12[i][4]);
				error = Math.abs((Double) outputs[0] - fn);
				
				
				if (!(error < BIG_NUMBER))
					error = BIG_NUMBER;
				else if (error < PROBABLY_ZERO) 
					error = 0.0;
			
				diff += error;
				
			}	

		}
		
		((FitnessCGP) ind.fitness).setFitness(state, diff, diff <= 0.01); // stop if error is less than 1%.
		ind.evaluated = true;
	}

	/** Koza 1: Fourth-order polynomial. */
	public static double koza1(double x) {
		return x * x * x * x + x * x * x + x * x + x;
	}

	/** Koza 2: Fifth-order polynomial. */
	public static double koza2(double x) {
		return x * x * x * x * x - 2 * x * x * x + x;
	}

	/** Koza 3: Sixth-order polynomial. */
	public static double koza3(double x) {
		return x * x * x * x * x * x - 2 * x * x * x * x + x * x;
	}

	/** Nguyen 4: Sixth-order polynomial. */
	public static double nguyen4(double x) {
		return x * x * x * x * x * x + x * x * x * x * x + x * x * x * x + x * x * x + x * x + x;
	}

	/** Nguyen 5 */
	public static double nguyen5(double x) {
		return Math.sin(Math.pow(x, 2)) * Math.cos(x) - 1;
	}

	/** Nguyen 6 */
	public static double nguyen6(double x) {
		return Math.sin(x) + Math.sin(x + x * x);
	}

	/** Nguyen 7 */
	public static double nguyen7(double x) {
		return Math.log(x + 1) + Math.log(x * x + 1);
	}

	/** Nguyen 8 */
	public static double nguyen8(double x) {
		return Math.sqrt(x);
	}

	/** Nguyen 9 */
	public static double nguyen9(double x, double y) {
		return Math.sin(x) + Math.sin(y * y);
	}

	/** Nguyen 10 */
	public static double nguyen10(double x, double y) {
		return 2 * Math.sin(x) * Math.cos(y);
	}
	
	/** Nguyen 11 */
	public static double nguyen11(double x, double y) {
		return Math.pow(x, y);
	}
	
	/** Nguyen 12 */
	public static double nguyen12(double x, double y) {
		return Math.pow(x, 4.0) - Math.pow(x, 3.0) + ( 0.5 * Math.pow(y, 2.0)) - y;
	}


	/** Keijzer 6 */
	public static double keijzer6(double x) {
		double sum = 0;
		double fx = Math.floor(x);
		for (int i = 1; i < fx + 1; i++) // up to and including floor(x)
			sum += (1.0 / i);
		return sum;
	}

	/** Pagie 1 */
	public static double pagie1(double x, double y) {
		return (1.0 / (1.0 + Math.pow(x, -4.0))) + (1.0 / (1.0 + Math.pow(y, -4.0)));
	}

	public static double vlad4(double[] vars) {
		double sum = 0;
		for (int i = 0; i < 5; i++)
			sum += (vars[i] - 3) * (vars[i] - 3);
		return 10.0 / (5.0 + sum);
	}

	public static double korns12(double x, double w) {
		return 2.0 - (2.1 * (Math.cos(9.8 * x) * Math.sin(1.3 * w)));
	}

	public void setup(EvolutionState state, Parameter base) {
		super.setup(state, base);

		Parameter def = defaultBase();

		function = state.parameters.getInt(base.push(P_WHICH), def.push(P_WHICH), 1);
		if (function == 0)
			state.output.fatal("problem.which must be present and > 0.");
		state.output.exitIfErrors();

		Benchmarks benchmarks = new Benchmarks();

		pagieTraining1 = benchmarks.generateIntervalSpacedSamples(state, new double[] { -5.0, -5.0 },
				new double[] { 5.0, 5.0 }, new double[] { 0.4, 0.4 }, 0);

		vladTraining4 = benchmarks.generateRandomSamples(state, new double[] { 0.05, 0.05, 0.05, 0.05, 0.05 },
				new double[] { 6.05, 6.05, 6.05, 6.05, 6.05 }, 1024, 0);
		
		kornsTraining12 = benchmarks.generateRandomSamples(state, new double[] { -50, -50, -50, -50, -50 },
				new double[] { 50, 50, 50, 50, 50 }, 10000, 0); 
		
		keijzerERC = new KeijzerERC();
		kornsERC = new KornsERC();

	}

}
