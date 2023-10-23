package DifferentialEvolution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class Population {
	private ArrayList<double []> _pool;
	private double[]   _fitness;
	private int    _numGene;
	private double _crossOverProb;
	private double _mutateF;
	private Simulator _simulator;
	
	public Population(double crossOverProb,double mutateF,Simulator simulator)	{
		_crossOverProb 	= crossOverProb;
		_mutateF 		= mutateF;
		_pool    = new ArrayList<double []> ();
		_fitness = null;
		_simulator = simulator;
		_numGene  = _simulator.getNumParameter();
	}
	private double [] generateRandomChromesom() {
		double [] chrom = new double [_numGene];
		//Individual member = new Individual();
		//member.generateRandomIndividual();	// generate random genes for chromosome
		
		for(int i=0;i<_numGene;i++) {
			//double num = Math.random()*2.0*Math.PI;
			double num = Math.random()*360.0; // 0-360.0 degree
			chrom[i] = num;
		}
		chrom[2] = chrom[1] * (Math.random()*0.8 +0.1);
		_simulator.generateRandomFeasibleChrom(chrom);
		return chrom;
	} 
    // *GenerateRandomPopulation(): Generate random population of n indivisual
	public void generateRandomPopulation(int numPop)	{
		_fitness = new double[numPop];
		for (int j = 0; j<numPop; j++) {
			double [] chrom = generateRandomChromesom();
			_pool.add(chrom);
			_fitness[j] = 0.0;
		}
	}
    // *EvaluateFitnessAll() : foreach indivisual x in _pool { x.EvaluateFitness(); }
	private double computeFitness(double [] chromX) {
		//SimulatorFourBar sim= new SimulatorFourBar(); 
		_simulator.setChromosome(chromX);
		double fit = _simulator.evaluateFitness();
		return fit;
	}
	
	public void evaluateFitnessAll()	{
		for(int j=0;j<_pool.size();j++) {
			double [] chromX = _pool.get(j);
			double fit = computeFitness(chromX);
			_fitness[j] = fit;
		}
	}
	
	public void statsFitnessAll()	{
		int numSurive   = 0;
		int numNearOK	= 0;
		int numPop = _pool.size();
		double  maxFit =-10000;
		double  [] maxChrom = null;
		
		double [] toSort = new double [_pool.size()];
		for(int j=0;j<numPop;j++) {
			double [] chromX = _pool.get(j);
			double fit = _fitness[j] ;
			toSort[j] = -fit;
			if(fit >=  0.0)	{	numSurive ++;			}
			if(fit > -1.0) {	numNearOK++;			}
			if(maxFit < fit) {
				maxFit = fit;
				maxChrom = chromX;
			}
		}
		Arrays.sort(toSort);
		int numOK=0;
		for(int j=0;j<numPop;j++) {
			double fit = -toSort[j];
			toSort[j] = fit;
			if(fit<-1.0) {  break; }
			numOK++;
		}
		//double maxFit = toSort[0];
		double q1Fit = toSort[numOK/4];
		double q2Fit = toSort[numOK/2];
		double q3Fit = toSort[numOK*3/4];
		double minFit = toSort[numOK-1];
				
		System.out.format("// #Surival=%.3f (%%); #OK=%d; M-Q1-Q2-Q3-m = [%.3f, %.3f, %.3f, %.3f, %.3f]\n",
							(100.0*numSurive)/numPop,
							numOK,
							maxFit,  
							q1Fit,
							q2Fit,
							q3Fit,
							minFit);
		String strMaxChrom = "let chrom1 = [";
		for(int i=0;i<maxChrom.length;i++) {	strMaxChrom += String.format("%.3E",maxChrom[i]) + ", ";	}
		strMaxChrom += "];";
		System.out.println(strMaxChrom);
	}
	
	public void evolveOneGeneration() {
		for(int j =0; j<_pool.size();j++) {
			double[] chromX = _pool.get(j);
			// 2-1-1: pick up 3 random Individuals/chromosomes a, b, c (a is usually best)
			int n1= (int)(Math.random()*_pool.size());
			int n2= (int)(Math.random()*_pool.size());
			int n3= (int)(Math.random()*_pool.size());
			int nA = n1, nB= n2, nC= n3;
			
			if(_fitness[n1]>=_fitness[n2]) {
				if(_fitness[n1]>=_fitness[n3]) { // n1> n3, n1 >n2
					nA=n1;  nB=n2; nC=n3;
				} else { // n3 biggest
					nA=n3;  nB=n1; nC=n2;
				}
			} else { // f[n1]< f[n2]
				if(_fitness[n2]>=_fitness[n3]) { // n2 biggest
					nA=n2;  nB=n1; nC=n3;
				} else { // n3 biggest
					nA=n3;  nB=n2; nC=n1;
				}
			}
			
			double [] chromA = _pool.get(nA);
			double [] chromB = _pool.get(nB);
			double [] chromC = _pool.get(nC);
			
			//2-1-2: pick random index R = random(0, n-1)
			int idxR= (int)(Math.random()*_numGene);

			// 2-1-3: compute new chromosome y as: 
			double [] chromY = new double[_numGene];
			for (int i=0;i<_numGene;i++) {
				double r = Math.random();
				if (r<this._crossOverProb || i==idxR) {
					double testGene = chromA[i] + _mutateF*(chromB[i]-chromC[i]);
					if(testGene<=0.0) { testGene = chromA[i]; }
					chromY[i] = testGene;
				} else { 
					chromY[i] = chromX[i];
				}
			}
			// 2-1-4: replacement:
			// compute chromY 's fitness
			double fitY = computeFitness(chromY);
			if (_fitness[j]<=fitY) {
				//replace y as x: replace population(x) = y;
				_pool.set(j,chromY);
				_fitness[j] = fitY;
			}
		}		
	}
}
