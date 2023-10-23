package DifferentialEvolution;

abstract class Simulator {
	public class Point 	{
		public double x;
		public double y;
		public Point(double x, double y){
			this.x = x;
			this.y = y;
		}
		public String toString(){
			return "(" + this.x + ", " + this.y + ")";
		}
	}
	/*
	public class Bar 	{
		public Point p1;
		public Point p2;
		public Bar(Point p1, Point p2){
			this.p1 = p1;
			this.p2 = p2;
		}
		public double length(){
			double dx = this.p1.x - this.p2.x;
			double dy = this.p1.y - this.p2.y;
			return Math.sqrt(dx*dx + dy*dy);
		}
	}
	*/
	protected Point   fromPoint(Point point, double distance, double angle){
		angle = angle * Math.PI / 180; // change to radian angle
		double x = point.x + distance * Math.cos(angle);
		double y = point.y + distance * Math.sin(angle);
		return new Point(x, y);
	}
	protected double  distance(Point p1, Point p2){
		// return distance of p1, p2
		double dx = p1.x - p2.x;
		double dy = p1.y - p2.y;
		return Math.sqrt(dx*dx + dy*dy);
	}
	protected boolean ccw(Point p1, Point p2, Point px){
		// is 3 point counter clock wise?
		double ax = p1.x; 
		double ay = p1.y;
		double bx = p2.x; 
		double by = p2.y;
		double cx = px.x; 
		double cy = px.y;
		return ((bx - ax) * (cy - ay) - (by - ay) * (cx - ax)) < 0;
	}
	
	protected Point link3rdPoint(Point p1, double len1, Point p2, double len2){
		// link 3rd point from 2 points and 2 lengths 
		// SSS uniqueness, 2 possibilities
		// find counter clock wise point 
		// step 1: determine if lengths are valid
		double len0 = distance(p1, p2);
		if(len0>=(len1+len2) || len0<=Math.abs(len1-len2)){ // triangle inequality
			return null;
		}
		double dx = p2.x - p1.x;
		double dy = p2.y - p1.y;
		double d = Math.sqrt(dx*dx + dy*dy); // // # distance between the centers
		double a = (len1*len1 - len2*len2 + d*d) / (2 * d); // # distance from first center to the radical line
		// # intersection of centerline with radical line
		double x = p1.x + dx * a / d; 
		double y = p1.y + dy * a / d;
		Point M = new Point (x, y);
		double h = Math.sqrt(len1*len1 - a*a); // # distance from the midline to any of the points
		// # There are two results, but only one (the correct side of the line) must be chosen
		double rx = -dy * h / d;
		double ry = dx * h / d;
		x = M.x + rx;
		y = M.y + ry;
		Point R1 = new Point (x, y);
		boolean test1 = ccw(p1, p2, R1);     // let test2 = ccw(p1, p2, R2);
		if(!test1) {
			R1.x = M.x - rx;
			R1.y = M.y - ry;
		}
		return R1;
	}	
	
	//private ArrayList <Gene> _chromosome;
	protected Point []  _locus;
	protected double [] _chromosome;
	
	protected Simulator() 	{
		_chromosome = null;
		_locus 		= null;
	}
	
	public void setChromosome(double [] chromosome) {
		 _chromosome = chromosome;
	}
	
	abstract protected int   getNumParameter();
	abstract protected Point computePosition(double angle);

	protected int  computeLocus() {
		normalizeChromosome();
		int count =0;
		_locus = new Point[360];
		for(int i=0; i<360; i++){
			Point p = computePosition((double) i);
			_locus[i] = p;
			if(p==null) { count++; 
			} else {
			  //System.out.println("angle [" + i + "] : (" +  p.x + " , " + p.y + ")" );
			}
		}
		return count;
	}
	
	protected double computeLiftForce(){
		// use _locus and drag force equation to estimate lift force
		// lift force ~ 0.5*rho*C*vy^2*Area
		// vy = _locus[i].y - _locus[i-1].y
		// Area = abs(_locus[i].x)
		// if vy<0, lift force>0
		// else lift force <0
		// summation lift force [i] as impulse to compute average lift force
		double totalLiftForce = 0.0;
		if (_locus == null){
			return -1.0E100;
		}
		for(int i = 0; i <360; i++){
			double dy = (i==359)?(_locus[0].y - _locus[359].y):(dy = _locus[i+1].y - _locus[i].y);
			double area = Math.abs(_locus[i].x);
			double force = dy*dy*area;
			if (dy>0){
				force = -force;
			}
			totalLiftForce += force;
		}
		return totalLiftForce;
	}
	protected double computeWeight(){
		// weight ~ total length of Bars (compute from _mu array)
		double total = 1.0;
		return total;
	}
	
	public void normalizeChromosome() {
		// find the total length
		// _chromosome [0] is thetaB; others are length
		double totalLength=0;
		for(int i=1;i<_chromosome.length;i++) {
			totalLength += _chromosome[i];
		}
		// normalize
		for(int i=1;i<_chromosome.length;i++) {
			_chromosome[i] /= totalLength;
		}		
	}

	public double evaluateFitness(){
		int failcount = computeLocus();
		if(failcount>0) { return -failcount; }
		double liftForce = computeLiftForce();
		double weight = computeWeight(); 
		double fitness = liftForce / weight;
		return fitness;
	}
}

class SimulatorFourBar extends Simulator {
	
	public SimulatorFourBar() {
		super();
	}
	protected int   getNumParameter() {
		return 17;
	}
	protected Point computePosition(double angle) {
		Point p0  = new Point(0,0);
		Point p1  = fromPoint(p0, _chromosome[1], _chromosome[0]);
		Point p2  = fromPoint(p0, _chromosome[2], angle);
		Point p3  = link3rdPoint(p1, _chromosome[3], p2, _chromosome[4]); if(p3==null)	{ return null; }
		Point p4  = link3rdPoint(p3, _chromosome[5], p2, _chromosome[6]); if(p4==null)	{ return null; }
		return p4;
	}
}

class SimulatorTheoJansen extends Simulator {
	
	public SimulatorTheoJansen() {
		super();
	}
	protected int   getNumParameter() {
		return 13;
	}
	protected Point computePosition(double angle) {
		Point p0 = new Point(0,0);
		Point p1 = fromPoint(p0, _chromosome[1], _chromosome[0]);
		Point p2 = fromPoint(p0, _chromosome[2], angle);
		Point p3 = link3rdPoint(p1, _chromosome[3], p2, _chromosome[4]);   if(p3==null)	{ return null; }
		Point p4 = link3rdPoint(p1, _chromosome[5], p3, _chromosome[6]);   if(p4==null)	{ return null; }
		Point p5 = link3rdPoint(p2, _chromosome[7], p1, _chromosome[8]);   if(p5==null)	{ return null; }
		Point p6 = link3rdPoint(p5, _chromosome[9], p4, _chromosome[10]);  if(p6==null)	{ return null; }
		Point p7 = link3rdPoint(p5, _chromosome[11], p6, _chromosome[12]); if(p7==null)	{ return null; }		
		return p7;
	}
}

class SimulatorFlapWingI extends Simulator {
	
	public SimulatorFlapWingI() {
		super();
	}
	protected int   getNumParameter() {
		return 15;
	}		
	protected Point computePosition(double angle) {
		Point p0  = new Point(0,0);
		Point p1  = fromPoint(p0, _chromosome[1], _chromosome[0]);
		Point p2  = fromPoint(p0, _chromosome[2], angle);
		Point p3 = link3rdPoint(p1, _chromosome[3], p2, _chromosome[4]);   if(p3==null)	{ return null; }
		Point p4 = link3rdPoint(p1, _chromosome[5], p3, _chromosome[6]);   if(p4==null)	{ return null; }
		Point p5 = link3rdPoint(p2, _chromosome[7], p1, _chromosome[8]);   if(p5==null)	{ return null; }
		Point p6 = link3rdPoint(p5, _chromosome[9], p1, _chromosome[10]);   if(p6==null)	{ return null; }
		Point p7 = link3rdPoint(p6, _chromosome[11], p4, _chromosome[12]);   if(p7==null)	{ return null; }
		Point p8 = link3rdPoint(p6, _chromosome[13], p7, _chromosome[14]);   if(p8==null)	{ return null; }
		return p8;
	}
}

class SimulatorFlapWingII extends Simulator {
	
	public SimulatorFlapWingII() {
		super();
	}
	protected int   getNumParameter() {
		return 11;
	}
	protected Point computePosition(double angle) {
		Point p0  = new Point(0,0);
		Point p1  = fromPoint(p0, _chromosome[1], _chromosome[0]);
		Point p2  = fromPoint(p0, _chromosome[2], angle);
		Point p3 = link3rdPoint(p2, _chromosome[3], p1, _chromosome[4]);   if(p3==null)	{ return null; }
		Point p4 = link3rdPoint(p2, _chromosome[5], p3, _chromosome[6]);   if(p4==null)	{ return null; }
		Point p5 = link3rdPoint(p4, _chromosome[7], p1, _chromosome[8]);   if(p5==null)	{ return null; }
		Point p6 = link3rdPoint(p4, _chromosome[9], p5, _chromosome[10]);   if(p6==null)	{ return null; }
		return p6;
	}
}

class TestSimulator {
	public static void testSimulatorFourBar () {
		SimulatorFourBar fbMec = new SimulatorFourBar();
		//double [] gene = {-3.14, 50.0, 10.0, 30.0, 65.0, 40.0, 90.0};
		double [] gene = {0, 16.8, 13.2, 26.7, 26.7, 50,50};
		fbMec.setChromosome(gene);
		double fit = fbMec.evaluateFitness();
		System.out.println("Fitness = " +100000*fit);
	}
	public static void testSimulatorStrandBeast () {
		SimulatorTheoJansen fbMec = new SimulatorTheoJansen();
		double [] gene = {1.57, 33.6, 11.23, 24.5, 36.2, 26.2, 33.5, 39.7, 24.4, 26.3, 24.8, 37.1, 45.2 };
		fbMec.setChromosome(gene);
		double fit = fbMec.evaluateFitness();
		System.out.println("Fitness = " +100000*fit);		
	}
	public static void testSimulatorFlapWingI () {
		SimulatorFlapWingI fbMec = new SimulatorFlapWingI();
		double [] gene = {0.75, 11.0, 4.7, 7.34, 9.9, 7.89, 7.88, 13.2, 9.4, 25.4, 33.0, 11.8, 30.2, 95.21, 95.9};
		fbMec.setChromosome(gene);
		double fit = fbMec.evaluateFitness();
		System.out.println("Fitness = " +100000*fit);			
	}	
	public static void testSimulatorFlapWingII () {
		SimulatorFlapWingII fbMec = new SimulatorFlapWingII();
		double [] gene = {1.57, 36.2,12.3,27.6,26.7,84.4,58.4,26.5,58.9,76.3,99.3};
		fbMec.setChromosome(gene);
		double fit = fbMec.evaluateFitness();
		System.out.println("Fitness = " +100000*fit);			
	}
}

package DifferentialEvolution;

import java.util.ArrayList;

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
			double num = Math.random()*2.0*Math.PI;
			chrom[i] = num;
		}
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
		int count=0;
		int countOK=0;
		double maxFit =-10000;
		for(int j=0;j<_pool.size();j++) {
			double [] chromX = _pool.get(j);
			double fit = computeFitness(chromX);
			_fitness[j] = fit;
			if(fit >=0) {
				//System.out.println("Fitness[" + j+ "]= " +100000*fit);
				count ++;
			}
			if(fit > -36.0) {
				countOK++;
			}
			if(maxFit<fit) {
				maxFit = fit;
			}
		}
		System.out.println("Num of Surival="+count + "; num of OK="+ countOK +   "; MaxFit = " + +100000*maxFit);
	}
	
	public void evolveOnGeenrate() {
		for(int j =0; j<_pool.size();j++) {
			double[] chromX = _pool.get(j);
			// 2-1-1: pick up 3 random Individuals/chromosomes a, b, c (a is usually best)
			int nA= (int)(Math.random()*_pool.size());
			double [] chromA = _pool.get(nA);
			int nB= (int)(Math.random()*_pool.size());
			double [] chromB = _pool.get(nB);
			int nC= (int)(Math.random()*_pool.size());
			double [] chromC = _pool.get(nC);
			
			//2-1-2: pick random index R = random(0, n-1)
			int idxR= (int)(Math.random()*_numGene);
			/*
			2-1-3: compute new chromosome y as: 
				foreach (i in (0, n-1)) {
					generate random number r = Math.rand();
					if (r<CR || i==R) {
						set y[i] = a[i] + F*(b[i]-c[i]);
					} else { y[i] = x[i];}
				}
			*/
			double [] chromY = new double[_numGene];
			for (int i=0;i<_numGene;i++) {
				double r = Math.random();
				if (r<this._crossOverProb || i==idxR) {
					chromY[i] = chromA[i] + _mutateF*(chromB[i]-chromC[i]);
				} else { chromY[i] = chromX[i];}
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
        int numPopulation = 100000;        // target: 10000
        int numGenerations = 100;         // target: 1000
        double 	crossover_probability = 0.05;
        double 	mutation_probability = 0.05;
        int 	numGen = fSim.getNumParameter(); 
        double 	mutateF = 0.8;
        Population oldPop = new Population(crossover_probability, mutateF,fSim);
        oldPop.generateRandomPopulation(numPopulation);
        oldPop.evaluateFitnessAll();
        /*
        for (int i = 0; i <= numGenerations; i++ ) {
			oldPop.evaluateFitnessAll();
			oldPop.evolveOnGeenrate();
        }
        */
        
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//SimulatorFlapWingII fSim = new SimulatorFlapWingII();
		// SimulatorFourBar
		// SimulatorTheoJansen
		// SimulatorFlapWingI
		// SimulatorFlapWingII
		
		DifferentialEvolution(new SimulatorFlapWingII());
		DifferentialEvolution(new SimulatorFlapWingII());
		DifferentialEvolution(new SimulatorFlapWingII());
		DifferentialEvolution(new SimulatorFlapWingII());
		DifferentialEvolution(new SimulatorFlapWingII());
		//TestSimulator.testSimulatorFourBar();
		//TestSimulator.testSimulatorStrandBeast ();
		//TestSimulator.testSimulatorFlapWingI ();
		//TestSimulator.testSimulatorFlapWingII ();
	}

}
