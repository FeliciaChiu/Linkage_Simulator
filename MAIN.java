package DifferentialEvolution;

/*
* x is a n dimension real vector which represents a chromosome/individual of a population: double [] chromosome;
* control param: 
	n = numGenes
	NP = numPopulation (>10n)
	CR = crossover probability [0,1] (~0.9)
	F = differential weight [0,2] (~0.9)
*main flow
	Step 1: initialize all Individuals/chromosomes with random position in search space
	Step 2: loop over generations()
		2-1: foreach x in Population { 
			2-1-1: pick up 3 random Individuals/chromosomes a, b, c (a is usually best)
			2-1-2: pick random index R = random(0, n-1)
			2-1-3: compute new chromosome y as: 
				foreach (i in (0, n-1)) {
					generate random number r = Math.rand();
					if (r<CR || i==R) {
						set y[i] = a[i] + F*(b[i]-c[i]);
					} else { y[i] = x[i];}
				}
			2-1-4: replacement: 
				if (x.fitness()<=y.fitness()) {
					replace y as x: replace population(x) = y;
				}
		}
	Step 3: pick most fit chromosome x in Population as best solution

 */
public class MAIN {
	public static void DifferentialEvolution(Simulator fSim) {
        int numPopulation = 10000;        // target: 10000
        int numGenerations = 20;         // target: 1000
        double 	crossover_probability = 0.05;
        double 	mutation_probability = 0.05;
        int 	numGen = fSim.getNumParameter(); 
        double 	mutateF = 0.8;
        Population oldPop = new Population(crossover_probability, mutateF,fSim);
        oldPop.generateRandomPopulation(numPopulation);
        oldPop.evaluateFitnessAll();
        oldPop.statsFitnessAll();
        
        for (int i = 0; i <= numGenerations; i++ ) {
			oldPop.evolveOneGeneration();
			// oldPop.evaluateFitnessAll();
			oldPop.statsFitnessAll();
        }
	}
	
	public static void main(String[] args) {
		//SimulatorFlapWingII fSim = new SimulatorFlapWingII();
		// SimulatorFourBar
		// SimulatorTheoJansen
		// SimulatorFlapWingI
		// SimulatorFlapWingII
		
		DifferentialEvolution(new SimulatorFlapWingIII());
		//DifferentialEvolution(new SimulatorFlapWingII());
		//DifferentialEvolution(new SimulatorFlapWingII());
		//DifferentialEvolution(new SimulatorFlapWingII());
		//DifferentialEvolution(new SimulatorTheoJansen());
		//TestSimulator.testSimulatorFourBar();
		//TestSimulator.testSimulatorStrandBeast ();
		//TestSimulator.testSimulatorFlapWingI ();
		//TestSimulator.testSimulatorFlapWingII ();
	}

}
