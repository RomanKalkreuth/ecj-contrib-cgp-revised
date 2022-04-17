package ec.cgp.representation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import ec.EvolutionState;
import ec.cgp.Stats;
import ec.util.MersenneTwisterFast;
import ec.util.Parameter;
import ec.vector.VectorIndividual;

/**
 * This class is an extension of the CGP IntegerVectorIndividual class. 
 * Contains a set of advanced genetic operators which have been proposed
 * over the last years. 
 * 
 * @author Roman Kalkreuth, roman.kalkreuth@tu-dortmund.de,
 *         https://orcid.org/0000-0003-1449-5131,
 *         https://ls11-www.cs.tu-dortmund.de/staff/kalkreuth,
 *         https://twitter.com/RomanKalkreuth
 *         
 * @author Jakub Husa, ihusa@fit.vut.cz
 * 		   https://www.vut.cz/en/people/jakub-husa-138342?aid_redir=1
 *         
 */
public class AdvancedIntegerVectorIndividual extends IntegerVectorIndividual {

	int maxActiveGenes;

	public ArrayList<Integer> activeFunctionNodes;
	ArrayList<Integer> passiveFunctionNodes;

	/**
	 * TODO Check if the overflow case is really needed
	 */
	public int randomValueFromClosedInterval(int min, int max, int val, MersenneTwisterFast random) {
		int l = 0;
		if (max - min < 0 || max == min) {
			return val;
		} else {

			do {
				l = min + random.nextInt(max - min + 1);
			} while (l == val);

			return l;
		}
	}

	/**
	 * Checks wether a certain gene is active or not. 
	 */
	public boolean geneActive(ArrayList<Integer> activeFunctionNodes, AdvancedIntegerVectorSpecies s, int nodeNum,
			int genePos) {
		return (activeFunctionNodes.contains(nodeNum) || s.phenotype(genePos, genome) == s.GENE_OUTPUT);
	}


	/*
	 * Simple genotypic one point crossover technique. Did not work very well in the past. 
	 * 
	 *  References: 
	 *  Miller (1999) http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.47.5554
	 *  Husa and Kalkreuth (2018) http://dx.doi.org/10.1007/978-3-319-77553-1_13
	 */
	public void onepointCrossover(EvolutionState state, int thread, AdvancedIntegerVectorIndividual ind)
	{
		IntegerVectorSpecies s = (IntegerVectorSpecies) species;// the current species, same for both parents
		AdvancedIntegerVectorIndividual i = (AdvancedIntegerVectorIndividual) ind;// the second parent
		int tmp;
		int point;

		// check that the chromosomes are equally long (true in most cases)
		int len = Math.min(genome.length, i.genome.length);
		if (len != genome.length || len != i.genome.length)
			state.output.warnOnce(
					"Genome lengths are not the same.  Vector crossover will only be done in overlapping region.");

		point = state.random[thread].nextInt((len / s.chunksize));// randomly select the point of crossover
		for (int x = 0; x < point * s.chunksize; x++)// swaps the first halves of the chromosomes
		{
			tmp = i.genome[x];
			i.genome[x] = genome[x];
			genome[x] = tmp;
		}

	}
	
	
	/*
	 * Adaption of discrete/uniform recombination in CGP. 
	 * 
	 * Work currently under review. Paper will be published hopefully soon.
	 * 
	*/
	public void discreteCrossover(EvolutionState state, int thread, AdvancedIntegerVectorIndividual ind)
	{
		AdvancedIntegerVectorSpecies s = (AdvancedIntegerVectorSpecies) species; // the current species, same for both parents
		AdvancedIntegerVectorIndividual i = (AdvancedIntegerVectorIndividual) ind;// the second parent
		
		int tmp = 0;
		
		// Node numbers are stored if two nodes are selected for the swap of the function gene 
		int swapNode1 = 0;
		int swapNode2 = 0;
		
		// Indices for the function genes
		int index1 = 0;
		int index2 = 0;
				
		boolean boundaryExpansion = true;
		
		// Determine the phenotypic length
		int len1 = activeFunctionNodes.size();
		int len2 = i.activeFunctionNodes.size();

		// check that the chromosomes are equally long (true in most cases)
		int min = Math.min(len1, len2);
		int max = Math.max(len1, len2);
	
		// Iterate over the minimum phenotype length
		for (int x = 0; x < min; x++) 
		{
			if(state.random[thread].nextBoolean()) {
				
				if(boundaryExpansion && x == (min - 1) && len1 != len2  ) {
					
					int r = state.random[thread].nextInt(max - x);
					
					if(len1 < len2) {
						swapNode1 = activeFunctionNodes.get(x);
						swapNode2 = i.activeFunctionNodes.get(x + r);
					} else {
						swapNode1 = activeFunctionNodes.get(x + r);
						swapNode2 = i.activeFunctionNodes.get(x);
					}
				
					
				} else {
					swapNode1 = activeFunctionNodes.get(x);
					swapNode2 = i.activeFunctionNodes.get(x);
				}
				
				// calculate the swap indexes
				index1 = (swapNode1 - s.numInputs) * (1 + s.maxArity);
				index2 = (swapNode2 - s.numInputs) * (1 + s.maxArity);
				
				// perform the swaps
				tmp = genome[index1];
				genome[index1] = i.genome[index2];
				i.genome[index2] = tmp;
			}			
		}
	}
	

	/*
	 * Determines a set of active function node by chance with respect to the predefined 
	 * maximum block size.
	 */
	public void determineSwapNodes(int blockSize, ArrayList<Integer> swapNodesList,
			ArrayList<Integer> activeFunctionNodes, EvolutionState state, int thread) {

		int j = 0;
		int randIndex;
		int nodeNumber;

		ArrayList<Integer> possibleNodes = new ArrayList<>(activeFunctionNodes);

		while (j < blockSize) {
			randIndex = state.random[thread].nextInt(possibleNodes.size());
			nodeNumber = possibleNodes.get(randIndex);

			swapNodesList.add(nodeNumber);

			possibleNodes.remove(randIndex);
			j++;
		}
	}

	
	/*
	 * Block crossover swaps blocks of active function genes between two individuals.
	 * 
	 * References: 
	 * Husa and Kalkreuth (2018) http://dx.doi.org/10.1007/978-3-319-77553-1_13
	 */
	public void blockCrossover(EvolutionState state, int thread, AdvancedIntegerVectorIndividual ind, int blockSize) {
		int swapNode1 = 0;
		int swapNode2 = 0;
		int id1_parent = 0;
		int id2_parent = 0;

		boolean debug = false;

		int j = 0;
		int temp = 0; // used to store values during swapping

		// memorize the individual's sepecies and the other individual
		AdvancedIntegerVectorSpecies s = (AdvancedIntegerVectorSpecies) species;
		AdvancedIntegerVectorIndividual i = (AdvancedIntegerVectorIndividual) ind;// the second parent
		
		ArrayList<Integer> swapNodesList1 = new ArrayList<Integer>();
		ArrayList<Integer> swapNodesList2 = new ArrayList<Integer>();

		s.determineActiveFunctionNodes(activeFunctionNodes, s, genome);
		s.determineActiveFunctionNodes(i.activeFunctionNodes, s, genome);


		if ((activeFunctionNodes.size() == 0) || (i.activeFunctionNodes.size() == 0)) {
			return;
		}

		if ((activeFunctionNodes.size() < blockSize) || (i.activeFunctionNodes.size() < blockSize)) {
			blockSize = Math.min(activeFunctionNodes.size(), i.activeFunctionNodes.size());
		}
		

		determineSwapNodes(blockSize, swapNodesList1, activeFunctionNodes, state, thread);
		determineSwapNodes(blockSize, swapNodesList2, i.activeFunctionNodes, state, thread);

		for (j = 0; j < blockSize; j++) {
			swapNode1 = swapNodesList1.get(j);
			swapNode2 = swapNodesList2.get(j);

			// calculate the swap indexes
			id1_parent = (swapNode1 - s.numInputs) * (1 + s.maxArity);
			id2_parent = (swapNode2 - s.numInputs) * (1 + s.maxArity);
		
			// perform the swaps
			temp = genome[id1_parent];
			genome[id1_parent] = i.genome[id2_parent];
			i.genome[id2_parent] = temp;

		}
	}
	
	/*
	 * 
	 */
	public void mergeParts(ArrayList<Integer> legalNodes, ArrayList<Integer> activeNodesPart1,
			ArrayList<Integer> activeNodesPart2, int crossoverPosition, int crossoverPoint, int[] genomePart1,
			int[] genomePart2, int[] newGenome, MersenneTwisterFast rand) {

		AdvancedIntegerVectorSpecies s = (AdvancedIntegerVectorSpecies) species;

		boolean connectedNeighbourNode = false;

		for (int i = 0; i < crossoverPosition; i++) {
			newGenome[i] = genomePart1[i];
		}

		for (int i = crossoverPosition; i < genomePart1.length; i++) {
			newGenome[i] = genomePart2[i];
		}

		for (int i = crossoverPosition; i < newGenome.length; i++) {
			if (s.phenotype(i, newGenome) == s.GENE_FUNCTION) {
				int nodeNum = s.nodeNumber(i, newGenome);

				// is node active ?
				if (activeNodesPart2.contains(nodeNum)) {
					for (int j = 1; j <= s.maxArity; j++) {
						int randArg;

						if (s.connectNeighbourNode && !connectedNeighbourNode) {
							if (crossoverPoint >= 1) {
								newGenome[i + j] = activeNodesPart1.get(crossoverPoint);
							} else {
								newGenome[i + j] = activeNodesPart1.get(0);
							}
							connectedNeighbourNode = true;
						} else if (!legalNodes.contains(newGenome[i + j])) {

							if (rand.nextBoolean(0.5)) {
								randArg = rand.nextInt(s.numInputs);
							} else {

								int max = legalNodes.indexOf(nodeNum);
								int min = s.numInputs;
								randArg = rand.nextInt((max - min)) + min;

							}

							newGenome[i + j] = legalNodes.get(randArg);
						}
					}
				}
			}
		}
		reconnectOutputs(newGenome, legalNodes, rand);
	}
	
	/*
	 * 
	 */
	public void reconnectOutputs(int[] genome, ArrayList<Integer> legalNodes, MersenneTwisterFast rand) {

		IntegerVectorSpecies s = (IntegerVectorSpecies) species;

		for (int i = 0; i < s.numOutputs; i++) {
			int nodeNum = s.nodeNumber(genome.length - i - 1, genome);

			int randArg;

			if (!legalNodes.contains(genome[genome.length - i - 1])) {
				int max = legalNodes.size();
				int min = s.numInputs;
				randArg = rand.nextInt((max - min)) + min;

				genome[genome.length - i - 1] = randArg;

			}

		}
	}

	/*
	 * Determines the list of legal nodes from the given sets of active function nodes.
	 */
	public void createLegalNodes(ArrayList<Integer> legalNodes, ArrayList<Integer> activeNodesPart1,
			ArrayList<Integer> activeNodesPart2, int crossoverNode) {
		IntegerVectorSpecies s = (IntegerVectorSpecies) species;

		// input nodes
		for (int i = 0; i < s.numInputs; i++) {
			legalNodes.add(i);
		}

		// Create legal nodes
		// legal node connections are
		// active nodes in phenotype 1
		for (int i = 0; i < activeNodesPart1.size(); i++) {
			if (activeNodesPart1.get(i) <= crossoverNode) {
				legalNodes.add(activeNodesPart1.get(i));
			}
		}

		// active nodes in phenotype 2 after cut point (c1)
		for (int i = 0; i < activeNodesPart2.size(); i++) {
			if (activeNodesPart2.get(i) > crossoverNode) {
				legalNodes.add(activeNodesPart2.get(i));
			}
		}

	}
	
	public void subgraphCrossover(EvolutionState state, int thread, VectorIndividual ind) {
		
		AdvancedIntegerVectorIndividual ind2 = (AdvancedIntegerVectorIndividual) ind;
		AdvancedIntegerVectorSpecies s = (AdvancedIntegerVectorSpecies) species;

		int[] genome1 = genome;
		int[] genome2 = ind2.genome;

		ArrayList<Integer> activeFunctionNodes1 = new ArrayList<Integer>();
		ArrayList<Integer> activeFunctionNodes2 = new ArrayList<Integer>();
		ArrayList<Integer> legalNodes = new ArrayList<Integer>();

		HashMap<Integer, Integer> activeNodesMap1 = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> activeNodesMap2 = new HashMap<Integer, Integer>();
		
		s.determineActiveFunctionNodes(activeFunctionNodes1, s, genome1);
		s.determineActiveFunctionNodes(activeFunctionNodes2, s, genome2);
		
		MersenneTwisterFast rand = state.random[thread];

		Collections.sort(activeFunctionNodes1);
		Collections.sort(activeFunctionNodes2);

		int c11;
		int c21;

		int c12;
		int c22;

		if (activeFunctionNodes1.size() == 0)
			return;

		if (activeFunctionNodes2.size() == 0)
			return;

		c11 = rand.nextInt(activeFunctionNodes1.size());
		c12 = rand.nextInt(activeFunctionNodes1.size());

		c21 = rand.nextInt(activeFunctionNodes2.size());
		c22 = rand.nextInt(activeFunctionNodes2.size());

		int[] newGenome1 = new int[(s.numNodes * (s.maxArity + 1)) + s.numOutputs];
		int[] newGenome2 = new int[(s.numNodes * (s.maxArity + 1)) + s.numOutputs];

		subgraphCrossover(c11, c21, activeFunctionNodes1, activeFunctionNodes2, legalNodes, genome1, genome2, newGenome1, state,
				thread);
		
		genome = newGenome1;
	}
	
	
	/*
	 * Subgraph recombination determines the active paths of both parents and randomly selects subgraphs of active function nodes
	 * from the phenotype. The selected subgraph are than exchanges between the parents. 
	 * 
	 * References: 
	 * Kalkreuth et al.(2017) http://dx.doi.org/10.1007/978-3-319-55696-3_19
	 * Kalkreuth (2021)  http://dx.doi.org/10.5220/0010110700590070
	 */
	public void subgraphCrossover(int c1, int c2, ArrayList<Integer> activeFunctionNodes1,
			ArrayList<Integer> activFunctioneNodes2, ArrayList<Integer> legalNodes, int[] genome1, int[] genome2,
			int[] newGenome, EvolutionState state, int thread) {

		IntegerVectorSpecies s = (IntegerVectorSpecies) species;
		MersenneTwisterFast rand = state.random[thread];

		int c1Node = activeFunctionNodes1.get(c1);
		int c2Node = activFunctioneNodes2.get(c2);
		int c1Pos = ((activeFunctionNodes1.get(c1) + 1) - s.numInputs) * (s.maxArity + 1);
		int c2Pos = ((activFunctioneNodes2.get(c2) + 1) - s.numInputs) * (s.maxArity + 1);

		if (c1Node <= c2Node) {

			createLegalNodes(legalNodes, activeFunctionNodes1, activFunctioneNodes2, c1Node);
			mergeParts(legalNodes, activeFunctionNodes1, activFunctioneNodes2, c1Pos, c1, genome1, genome2, newGenome,
					rand);
		} else {

			createLegalNodes(legalNodes, activFunctioneNodes2, activeFunctionNodes1, c2Node);
			mergeParts(legalNodes, activFunctioneNodes2, activeFunctionNodes1, c2Pos, c2, genome2, genome1, newGenome,
					rand);
		}
	}

	/*
	 * Single active gene mutation strategy. The genome is randomly mutated by point mutation until
	 * exactly one active has been hit. 
	 * 
	 * References: 
	 *
	 * Goldman and Punch (2014) http://dx.doi.org/10.1109/TEVC.2014.2324539
	 * 
	 */
	public void singleActiveGeneMutation(EvolutionState state, int thread) {
		AdvancedIntegerVectorSpecies s = (AdvancedIntegerVectorSpecies) species;
		ArrayList<Integer> activeFunctionNodes = new ArrayList<Integer>();
		s.determineActiveFunctionNodes(activeFunctionNodes, s, genome);
		int nodeNum;
		int genePos;
		int geneVal;
		boolean hitActiveGene = false;

		do {
			genePos = state.random[thread].nextInt(genome.length);
			geneVal = genome[genePos];
			genome[genePos] = randomValueFromClosedInterval(0, s.computeMaxGene(genePos, genome), geneVal,
					state.random[thread]);

			nodeNum = s.nodeNumber(genePos, genome);

			if (geneActive(activeFunctionNodes, s, nodeNum, genePos)) {
				hitActiveGene = true;
			}

		} while (!hitActiveGene);
	}
	
	/**
	 * 
	 */
	public void multiActiveGeneMutation(EvolutionState state, int thread, int num) {
		AdvancedIntegerVectorSpecies s = (AdvancedIntegerVectorSpecies) species;
		ArrayList<Integer> activeFunctionNodes = new ArrayList<Integer>();
		s.determineActiveFunctionNodes(activeFunctionNodes, s, genome);
		int nodeNum;
		int genePos;
		int geneVal;
		int activeGenesHit = 0;

		do {
			genePos = state.random[thread].nextInt(genome.length);
			geneVal = genome[genePos];
			genome[genePos] = randomValueFromClosedInterval(0, s.computeMaxGene(genePos, genome), geneVal,
					state.random[thread]);

			nodeNum = s.nodeNumber(genePos, genome);

			if (geneActive(activeFunctionNodes, s, nodeNum, genePos)) {
				activeGenesHit++;
			}

		} while (activeGenesHit < num);

	}


	/**
	 * Mutate the genome. Adapted from IntegerVectorIndividual. The acceptable value
	 * range for each position is determined by CGPVectorSpecies.computeMaxGene.
	 */
	public void pointMutation(EvolutionState state, int thread) {
		IntegerVectorSpecies s = (IntegerVectorSpecies) species;
		for (int x = 0; x < genome.length; x++)
			if (state.random[thread].nextBoolean(s.mutationProbability(x))) {
				genome[x] = randomValueFromClosedInterval(0, s.computeMaxGene(x, genome), state.random[thread]);
			}
	}


	
	/*
	 * Insertion mutation activates exactly one inactive function node
	 * 
	 * References: Kalkreuth (2019) http://dx.doi.org/10.5220/0008070100820092
	 * Kalkreuth (2021) http://dx.doi.org/10.1007/978-3-030-70594-7_4
	 */
	public void insertion(EvolutionState state, int thread, int[] genome) {
		AdvancedIntegerVectorSpecies s = (AdvancedIntegerVectorSpecies) species;
		MersenneTwisterFast rand = state.random[thread];

		s.determineActiveFunctionNodes(activeFunctionNodes, s, genome);

		if (activeFunctionNodes.size() >= s.maxNodeInsertion) {
			return;
		}

		passiveFunctionNodes.clear();

		for (int i = (s.numInputs); i <= (s.numNodes + s.numInputs - 1); i++) {
			if (!activeFunctionNodes.contains(i)) {
				passiveFunctionNodes.add(i);
			}
		}

		if (passiveFunctionNodes.size() == 0)
			return;

		int j = rand.nextInt(passiveFunctionNodes.size());

		int mNodeNumber = passiveFunctionNodes.get(j);
		int mNodePosition = s.positionFromNodeNumber(mNodeNumber);

		boolean hasLeftNode = false;
		boolean hasRightNode = false;

		int rightNodeNumber = -1;
		int rightNodePosition = -1;

		Collections.sort(activeFunctionNodes);

		if (activeFunctionNodes.size() > 0) {
			if (mNodeNumber > activeFunctionNodes.get(0))
				hasLeftNode = true;

			if (mNodeNumber < (activeFunctionNodes.get(activeFunctionNodes.size() - 1)))
				hasRightNode = true;
		}

		if (hasRightNode) {
			int i = 0;
			int currentNode = activeFunctionNodes.get(i);
			while (currentNode < mNodeNumber) {
				currentNode = activeFunctionNodes.get(i);
				i++;
			}

			rightNodeNumber = currentNode;
			rightNodePosition = s.positionFromNodeNumber(rightNodeNumber);

		}

		if (hasRightNode) {
			int currentPosition;
			for (int i = 1; i <= s.maxArity; i++) {
				currentPosition = mNodePosition + i;
				genome[currentPosition] = genome[rightNodePosition + i];
			}

			int randomInputNumber = rand.nextInt(s.maxArity);
			genome[rightNodePosition + randomInputNumber + 1] = mNodeNumber;

		} else if (hasLeftNode) {
			int leftNodeNumber = activeFunctionNodes.get(activeFunctionNodes.size() - 1);
			for (int i = 1; i <= s.numOutputs; i++) {
				if (genome[genome.length - i] == leftNodeNumber) {
					genome[genome.length - i] = mNodeNumber;
					break;
				}
			}

			int currentPosition;

			for (int i = 1; i <= s.maxArity; i++) {
				currentPosition = mNodePosition + i;

				if (i == 1) {
					genome[currentPosition] = leftNodeNumber;
				} else {
					if (rand.nextBoolean()) {
						genome[currentPosition] = activeFunctionNodes.get(rand.nextInt(activeFunctionNodes.size()));
					} else {
						genome[currentPosition] = rand.nextInt(s.numInputs);
					}
				}

			}
		} else {
			int output = (rand.nextInt(s.numOutputs)) + 1;
			genome[genome.length - output] = mNodeNumber;

			int currentPosition;
			for (int i = 1; i <= s.maxArity; i++) {
				currentPosition = mNodePosition + i;
				genome[currentPosition] = rand.nextInt(s.numInputs);
			}
		}

		passiveFunctionNodes.remove(j);
		activeFunctionNodes.add(mNodeNumber);

	}

	/*
	 * Deletion mutation deactivates the first active function node 
	 * 
	 * References: 
	 * Kalkreuth (2019) http://dx.doi.org/10.5220/0008070100820092
	 * Kalkreuth (2021) http://dx.doi.org/10.1007/978-3-030-70594-7_4
	 */
	public void deletion(EvolutionState state, int thread, int[] genome) {

		AdvancedIntegerVectorSpecies s = (AdvancedIntegerVectorSpecies) species;
		MersenneTwisterFast rand = state.random[thread];

		s.determineActiveFunctionNodes(activeFunctionNodes, s, genome);

		if (activeFunctionNodes.size() <= s.minNodeDeletion) {
			return;
		}

		Collections.sort(activeFunctionNodes);

		int mNode = activeFunctionNodes.get(0);
		activeFunctionNodes.remove(0);

		int currentNode;
		int currentIndex;
		int randIndex;

		for (int i = 0; i < genome.length; i++) {
			if (s.phenotype(i, genome) == s.GENE_ARGUMENT) {
				currentNode = s.nodeNumber(i, genome);
				if (activeFunctionNodes.contains(currentNode) && genome[i] == mNode) {
					if (rand.nextBoolean()) {
						currentIndex = activeFunctionNodes.indexOf(currentNode);
						if (currentIndex == 0) {
							genome[i] = rand.nextInt(s.numInputs);
						} else {
							randIndex = rand.nextInt(currentIndex);
							genome[i] = activeFunctionNodes.get(randIndex);
						}

					} else {
						genome[i] = rand.nextInt(s.numInputs);
					}
				}
			}

			if (s.phenotype(i, genome) == s.GENE_OUTPUT) {
				if (genome[i] == mNode) {
					if (rand.nextBoolean()) {
						randIndex = rand.nextInt(activeFunctionNodes.size());
						genome[i] = activeFunctionNodes.get(randIndex);
					} else {
						genome[i] = rand.nextInt(s.numInputs);
					}
				}
			}

		}

	}
	
	/**
	 * First, the maximum possible depth for the duplication and inversion mutation is determined. The depth
	 * is than chosen by chance in respect to the maximum.
	 */
	public int stochasticDepth(EvolutionState state, int thread, int maxDepth, int numActiveFunctionNodes) {
		int depth;
		int max;

		if (numActiveFunctionNodes <= maxDepth) {
			max = numActiveFunctionNodes - 1;
		} else {
			max = maxDepth;
		}

		depth = state.random[thread].nextInt(max) + 1;

		return depth;
	}

	/**
	 * Determines a suitable start index for the duplication and inversion mutation by chance 
	 * and in respect to the number of active function nodes. 
	 */
	public int startIndex(EvolutionState state, int thread, int numactiveFunctionNodes, int depth) {
		int startMax = numactiveFunctionNodes - depth;
		int start;

		if (startMax <= 0) {
			start = 0;
		} else {
			start = state.random[thread].nextInt(startMax);
		}

		return start;
	}
	
	/**
	 * Phenotyoic inverison mutation: 
	 * Inverts the order of function genes of a randomly selected set of active nodes. The size of the set 
	 * is determined by chance and in respect to the number of active nodes.
	 *  
	 * Kalkreuth (2022): Phenotypic Duplication and Inversion in Cartesian Genetic Programming applied to Boolean Function Learning
 	 * (accepted for poster presentation at GECCO’22)
	 */
	public void inversion(EvolutionState state, int thread) {
		AdvancedIntegerVectorSpecies s = (AdvancedIntegerVectorSpecies) species;

		if (!state.random[thread].nextBoolean(s.inversionProbability)) {
			return;
		}

		int depth;
		int start;
		int end;
		int leftNode;
		int leftPosition;
		int rightNode;
		int rightPosition;
		int middle;
		int tmp;

		boolean debug = false;
		
		int numactiveFunctionNodes = activeFunctionNodes.size();

		Collections.sort(activeFunctionNodes);

		if (numactiveFunctionNodes <= 1) {
			return;
		}

		depth = stochasticDepth(state, thread, s.maxInversionDepth, numactiveFunctionNodes);
		start = startIndex(state, thread, numactiveFunctionNodes, depth);

		end = start + depth;
		middle = (int) Math.round(depth / 2.0);

		for (int i = 0; i < middle; i++) {
			leftNode = activeFunctionNodes.get(start + i);
			rightNode = activeFunctionNodes.get(end - i);

			leftPosition = s.positionFromNodeNumber(leftNode);
			rightPosition = s.positionFromNodeNumber(rightNode);

			tmp = genome[leftPosition];
			genome[leftPosition] = genome[rightPosition];
			genome[rightPosition] = tmp;
		}
	}

	/**
	 * Phenotyoic duplication mutation: 
	 * Duplicates the function gene of a randomly selected active node to a following sequence of active nodes. 
	 * The size of the sequence is determined by chance and in respect to the number of active nodes.
	 *  
	 * Kalkreuth (2022): Phenotypic Duplication and Inversion in Cartesian Genetic Programming applied to Boolean Function Learning
 	 * (accepted for poster presentation at GECCO’22)
	 */
	public void duplication(EvolutionState state, int thread) {
		AdvancedIntegerVectorSpecies s = (AdvancedIntegerVectorSpecies) species;

		if (!state.random[thread].nextBoolean(s.duplicationProbability)) {
			return;
		}
		
		int depth;
		int start;
		int end;
		int position;
		int node;
		int tmp;
		int function;

		boolean debug = false;
		int numactiveFunctionNodes = activeFunctionNodes.size();
		
		Collections.sort(activeFunctionNodes);

		if (numactiveFunctionNodes <= 1) {
			return;
		}

		depth = stochasticDepth(state, thread, s.maxDuplicationDepth, numactiveFunctionNodes);
		start = startIndex(state, thread, numactiveFunctionNodes, depth);
		end = start + depth;

		node = activeFunctionNodes.get(start);
		position = s.positionFromNodeNumber(node);
		function = genome[position];

		for (int i = start + 1; i <= end; i++) {
			node = activeFunctionNodes.get(i);
			position = s.positionFromNodeNumber(node);
			genome[position] = function;
		}
	}


	public void defaultMutate(EvolutionState state, int thread) {
		
		AdvancedIntegerVectorSpecies s = (AdvancedIntegerVectorSpecies) species;
		
		if(s.mutationType == s.C_POINT) {
			pointMutation(state, thread);
		} else if (s.mutationType == s.C_SINGLE) {
			singleActiveGeneMutation(state, thread);
		} else if (s.mutationType == s.C_MULTI) {
			multiActiveGeneMutation(state, thread, s.mutateActiveGenes);
		}
		
		if(s.insertionProbability > 0.0f) {
			insertion(state, thread, genome);
		}
		
		if(s.deletionProbability > 0.0f) {
			deletion(state, thread, genome);	
		}
		
		if(s.inversionProbability > 0.0f || s.duplicationProbability > 0.0f) {
			s.determineActiveFunctionNodes(activeFunctionNodes, s, genome);	
		}
		
		if(s.inversionProbability > 0.0f) {
			inversion(state, thread);	
		}
		
		if(s.duplicationProbability > 0.0f) {
			duplication(state, thread);
		}

	}
	
	/** Make a full copy of this individual. */
	public Object clone() {
		AdvancedIntegerVectorIndividual myobj = (AdvancedIntegerVectorIndividual) (super.clone());

		if (activeFunctionNodes != null) {
			myobj.activeFunctionNodes = new ArrayList<Integer>();
		}
		
		if (passiveFunctionNodes != null) {
			myobj.passiveFunctionNodes = new ArrayList<Integer>();
		}

		return myobj;
	}
	
	
	public void setup(final EvolutionState state, final Parameter base) {
		super.setup(state, base); // actually unnecessary unless
		activeFunctionNodes = new ArrayList<Integer>();
		passiveFunctionNodes = new ArrayList<Integer>();
		
	}


}
