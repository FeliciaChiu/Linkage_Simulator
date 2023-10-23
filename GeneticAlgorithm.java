/*
main procedure
{
    // S1: [Start] Generate random population of n indivisuals
    Population pop = new Population(); // empty population
    pop.GenerateRandomPopulation();

    for(int g=0;g<numGeneration;g++) {
        // S2: [Fitness] Evaluate the fitness f(x) of each indivisual x in the population
        pop.EvaluateFitnessAll();
        
        //S3: [New population] Create a new population by repeating following steps until the new population is complete
        Population new_pop = new Population ();
        for(amount of new child to be added) {
            //  S3-1: [Selection] Select two parent indivisuals from a population according to their fitness (the better fitness, the bigger chance to be selected)
            parent1 = pop.Select();
            parent2 = pop.Select();

            //  S3-2: [Crossover] With a crossover probability cross over the parents to form a new offspring (children). If no crossover was performed, offspring is an exact copy of parents.
            child = pop.Crossover(parent1,parent2,crossover_probability);
            
            //  S3-3: [Mutation] With a mutation probability mutate new offspring at each locus (position in chromosome).
            child.Mutate(mutation_probability);
            
            //  S3-4: [Accepting] Place new offspring in a new population
            new_pop.Accept(child); // adding to the popluation pool
        }
        // S4: [Replace] Use new generated population for a further run of algorithm
        pop = new_pop;
        // S5: [Test] If the end condition is satisfied, stop, and return the best solution in current population
    }
}
*/

import java.util.ArrayList; 

abstract class Gene
{
	// member: _toppology and _diemnsion of mechanism
    
	//*gene() : empty
	public Gene()	{ }
	
	//*GenerateRandom(): Generate random _toppology and _diemnsion
	public void generateRandom()	{	}
    
	//*Mutate(mutation_probability): mutate _toppology  and _diemnsion	
    public void mutate(double mutation_probability)	{	}
}

class BaseFrameGene extends Gene
{
	private double _lenB;
	private double _thetaB;
	
	public BaseFrameGene()
	{
		_lenB = 0.0;
		_thetaB = 0.0;
	}
	
	//*GenerateRandom(): Generate random _toppology and _diemnsion
	public void generateRandom()
	{
		
	}
    
	//*Mutate(mutation_probability): mutate _toppology  and _diemnsion	
    public void mutate(double mutation_probability)
	{
		
	}	
}

class CrankGene extends Gene
{
	private double _lenC;
	
	public CrankGene()	{
		_lenC = 0.0;
	}
	
	//*GenerateRandom(): Generate random _toppology and _diemnsion
	public void generateRandom()	{
		
	}
    
	//*Mutate(mutation_probability): mutate _toppology  and _diemnsion	
    public void mutate(double mutation_probability)	{
		
	}	
}

/* Topology gene
   construct 3rd point from two previous points.
   so, we need 	_fromPoint1 (index of _chromosome), and length 1
				_fromPoint2 (index of _chromosome), and length 2
*/

class TopoGene extends Gene
{
	private double _len1;
	private double _len2;
	private int _fromPoint1;
	private int _fromPoint2;
	
	public TopoGene() 	{
		_len1 = _len2 = 0.0;
		_fromPoint1 = _fromPoint2 = 0;
	}
	
	//*GenerateRandom(): Generate random _toppology and _diemnsion
	public void generateRandom() {
		
	}
    
	//*Mutate(mutation_probability): mutate _toppology  and _diemnsion	
    public void mutate(double mutation_probability) {
		
	}	
}


class Simulator
{
    private ArrayList <Gene> _chromosome;
	public Simulator(ArrayList <Gene> chromosome) 	{
		// MechSimulator(chromosome): setting chromosome
        _chromosome = chromosome;
	}
	
    public double evaluateFitness()	{
		return 0.0;
	}
}


class Individual
{
	// member: _chromosome (list of genes)
	private ArrayList <Gene> _chromosome;
	private double _fitness;		// store result of evaluateFitness()
     
	// *indivisual(): empty genes;
	public Individual()	{
		_chromosome = new ArrayList <Gene> ();
		_fitness = 0.0;
	}

	public double getFitness() {
		return _fitness;
	}
	// *GenerateRandom(): Generate random genes for th chromosome
	public void generateRandom()	{
		// BaseFrameGene for point B
		BaseFrameGene bGene = new BaseFrameGene();
		bGene.generateRandom();
		_chromosome.add(bGene);
		// CrankGene for point C
		CrankGene cGene = new CrankGene();
		cGene.generateRandom();
		_chromosome.add(cGene);
		// TopoGene for point D
		TopoGene dGene = new TopoGene();
		dGene.generateRandom();
		_chromosome.add(dGene);
		// TopoGene for point E
		TopoGene eGene = new TopoGene();
		eGene.generateRandom();
		_chromosome.add(eGene);
	}
    
	// *evaluateFitness(): evalute fitness : using Mechinism(_chromosome).evaluateFitness()
	public double evaluateFitness()	{
		Simulator sim = new Simulator(_chromosome);
		_fitness = sim.evaluateFitness();
		return _fitness;
	}
	
	// *crossover(parent2,crossover_probability): crossover with parent and return a new child.
    public Individual crossover(Individual parent2, double crossover_probability)	{
		return null;
	}
	
	// *Mutate(mutation_probability): foreach gene x in chromosome { x.Mutate(mutation_probability); }
    public void mutate(double mutation_probability)	{
		for(Gene x:_chromosome) {
			x.mutate(mutation_probability);
		}
	}
}

class Population
{
	// member: _pool: colletion of  indivisual
	private ArrayList <Individual> _pool;
    
	// *population():  create empty collections
	public Population()	{
		_pool = new ArrayList <Individual> ();
	}
	
    // *GenerateRandomPopulation(): Generate random population of n indivisual
	public void generateRandomPopulation(int n)	{
		for (int i = 0; i<n; i++) {
			Individual member = new Individual();
			member.generateRandom();	// generate random genes for chromosome
			_pool.add(member);
		}
	}
    // *EvaluateFitnessAll() : foreach indivisual x in _pool { x.EvaluateFitness(); }
	public void evaluateFitnessAll()	{
		for (Individual x:_pool) {
			x.evaluateFitness();
		}
	}

	// rank individuals by fitness: sort 
	public void rankIndividuals() {
		// sort fitness:Individual.getFitness()
	}
	
    // *Select(): randomly select a parent from _pool, and return it.
	public Individual select() 	{
		// select from ranking pool
		return null;
	}
	
    // *Crossover(parent1,parent2,crossover_probability): return parent1.crossover(parent2,crossover_probability).
	public Individual crossover(Individual parent1, Individual parent2, double crossover_probability)	{
		Individual child = parent1.crossover(parent2,crossover_probability);
		return child;
    }
	
	// *Accept(child): _pool.add(child);			
	public void accept(Individual child) {
		_pool.add(child);
	}
}

public class GeneticAlgorithm {
    public static void main(String [] args) {
        /*S1: [Start] Generate random population of n indivisuals
        (suitable solutions for the problem)
        loop over one generation {
        S2: [Fitness] Evaluate the fitness f(x) of each indivisual x in the population
        S3: [New population] Create a new population by repeating following steps until the new population is complete
            S3-1: [Selection] Select two parent indivisuals from a population according to their fitness (the better fitness, the bigger chance to be selected)
            S3-2: [Crossover] With a crossover probability cross over the parents to form a new offspring (children). If no crossover was performed, offspring is an exact copy of parents.
            S3-3: [Mutation] With a mutation probability mutate new offspring at each locus (position in chromosome).
            S3-4: [Accepting] Place new offspring in a new population
        S4: [Replace] Use new generated population for a further run of algorithm
        S5: [Test] If the end condition is satisfied, stop, and return the best solution in current population
        */
        int numPopulation = 1000;                       // target: 10000
        int numGenerations = 100;                       // target: 1000
        double crossover_probability = 0.05;
        double mutation_probability = 0.05;
        Population oldPop = new Population();
        oldPop.generateRandomPopulation(numPopulation);
        for (int i = 0; i <= numGenerations; i++ ) {
           
			oldPop.evaluateFitnessAll();
			oldPop.rankIndividuals();
            Population newPop = new Population();

            for (int j = 0; j< numPopulation; j++) {
                Individual parent1 = oldPop.select();
                Individual parent2 = oldPop.select();
                Individual child = oldPop.crossover(parent1, parent2, crossover_probability);
                child.mutate(mutation_probability);
                newPop.accept(child);
            }
            oldPop = newPop;
        }
        // now have most evolved population in oldPop
    }
}




