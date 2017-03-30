package FP;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.IntStream;

import MCDP.benchmark.Statistics;
import MCDP.model.MCDPData;
import MCDP.model.MCDPModel;
import MCDP.model.MCDPRandomSolution;
import MCDP.model.Solution;

public class MonarchOptimization {
	// Parametros de la metaheuristica
	private int numberPoblation;
	private int NP1;
	private int NP2;
	private int numberIteration;
	private double SMax;
	private double BAR;
	private double peri;
	private double p;
	
	private int vector_fitness[];
	Grafico grafico = null;
	private String directoryName = null;
	private ArrayList<Solution> poblation; // Un arreglo de soluciones
	private Solution bestSolution, tempSolution,worstSolution,tempSolution2,auxSolution; // Mejor Solucion
	int tempFitness = 0;
	private int cont =0;
	private double promedio = 0f;
	int backFitness=0;
	int iterationAutonomousSearch=5;
	int poblationIncrease=5;

	private int matrizSimilitud[][];
	private int modoDelta = 1;
	private double varMin = 9999f;
	private double varMax = 0f;
	// Dataset (benchmark)
	private MCDPData data;

	// Estadisticas
	// veces que un random genera una solucion que cumple las restricciones
	private long numAcceptedMoves;
	private long numRejectedMoves;

	Random rn;
	double d;
	Scanner sc = new Scanner(System.in);
	
	public MonarchOptimization(int numberPoblation, int numberIteration, MCDPData data, double SMax,
			double BAR,double peri,double p, String directory) {
		this.numberPoblation = numberPoblation;
		this.numberIteration = numberIteration;
		this.data = data;
		this.poblation = new ArrayList<Solution>();
		this.bestSolution = new Solution();
		this.worstSolution = new Solution();
		this.numAcceptedMoves = 0;
		this.numRejectedMoves = 0;
		this.SMax = SMax;
		this.BAR = BAR;
		this.peri = peri;
		this.p = p;
		this.vector_fitness = new int[numberIteration];
		this.directoryName = directory;
	}

	public double[] run( int ejecucion, int ejecuciones) {
		
		matrizSimilitud = new int[data.M][data.M];
		calcularSimilitudMaquinas();
		generateInitialPoblation();
		chooseBestSolutionInPoblation();
		int iteration = 0;
		int iterationOpt = 0;
		int optimo = 9999999;
		double[] var = new double[3];
		int iterationEstancadap = 0;
		rn = new Random();
		//PrintSolutions(0);
		while (iteration < this.numberIteration) {

		/*	try {
				System.in.read();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	PrintSolutions(0);*/
			NP1 = (int) Math.round(p*(poblation.size())); //calcular poblacion land1
			NP2 = poblation.size() - NP1; 		// calcular poblacion land2
			numberPoblation=poblation.size();
				for(int i=0;i<NP1;i++){		//Para la primera poblacion		
				tempSolution = new Solution(poblation.get(i).getMachine_cell(), poblation.get(i).getPart_cell(),poblation.get(i).getFitness());
				GenerateMovementMigrationOperator(i);	//primer movimiento
			
					if (tempFitness < poblation.get(i).getFitness()) {  //actualizar segun el algoritmo ""greedy""
						poblation.get(i).setMachine_cell(tempSolution.getMachine_cell());
						poblation.get(i).setPart_cell(tempSolution.getPart_cell());
						poblation.get(i).setFitness(tempFitness);
					}
				}
				chooseBestSolutionInPoblation(); //actualizar la mejor solucion
				chooseWorstSolutionInPoblation(); //actualizar la peor solucion
			
				for(int i=NP1;i<poblation.size();i++){
					tempSolution = new Solution(poblation.get(i).getMachine_cell(), poblation.get(i).getPart_cell(),poblation.get(i).getFitness()); 
					tempSolution2 = new Solution(poblation.get(i).getMachine_cell(), poblation.get(i).getPart_cell(),poblation.get(i).getFitness());
					GenerateMovementButterflyAdjustingOperator(i,iteration);
					generateCrossoverOperator(i);

					if (tempSolution2.getFitness()<poblation.get(i).getFitness()&&tempSolution2.getFitness() < tempFitness) { //dependiendo de cual sea mejor es la que se queda
						poblation.get(i).setMachine_cell(tempSolution2.getMachine_cell());
						poblation.get(i).setPart_cell(tempSolution2.getPart_cell());
						poblation.get(i).setFitness(tempSolution2.getFitness());					
					}
					if(tempFitness<poblation.get(i).getFitness()){
						poblation.get(i).setMachine_cell(tempSolution.getMachine_cell());
						poblation.get(i).setPart_cell(tempSolution.getPart_cell());
						poblation.get(i).setFitness(tempFitness);					
					}
				}
			vector_fitness[iteration] = bestSolution.getFitness();
			iteration++;
			if(optimo> bestSolution.getFitness()){
				optimo=bestSolution.getFitness();
				iterationOpt=iteration;	
				iterationEstancadap=0;
				
			}else{
				iterationEstancadap++;
			}
			//PrintSolutions(iteration);
		//	iterationEstancadap=AutonomousSearchP(iterationEstancadap);	
		}	
	//	Statistics.createConvergenciGraph(data.getIdentificator(), vector_fitness, directoryName,ejecucion,ejecuciones);
		var[0]=iterationOpt;
		var[1]=varMin;
		var[2]=varMax;
		return var;		 
	}
	
	void GenerateMovementMigrationOperator(int butterfly){
		boolean constraintOK = false;
		int[] vectorMaquinaActual = new int[data.M];
		int[] vectorAux = new int[data.M];
		vectorMaquinaActual = matrixToVector(tempSolution.getMachine_cell());
	//	System.out.println("\nAnterior: "+Arrays.toString(matrixToVector(tempSolution.getMachine_cell()))+" Fit "+tempSolution.getFitness());
		while (constraintOK == false) {	
			vectorAux = vectorMaquinaActual;
			for (int i = 0; i < data.M; i++) { //para todos los elementos de la mariposa
				int randomNum = ThreadLocalRandom.current().nextInt(-100, 100);
				double r=randomNum*peri;
				 Random randomGenerator = new Random();
				 
				if(r<=p){ //moviemiento 1
					 int randomInt = randomGenerator.nextInt(NP1); //random de la primera poblacion
					 vectorAux[i] = GetVectorValue(i,poblation.get(randomInt).getMachine_cell());					 
				}else{//moviemiento 2
					int randomInt = ThreadLocalRandom.current().nextInt(NP1, poblation.size());//random de la segunda poblacion
					 vectorAux[i] = GetVectorValue(i,poblation.get(randomInt).getMachine_cell());				
				}
			}
			
			constraintOK=CheclConstraints(butterfly);
		}
		tempSolution = vectorToMatrix(tempSolution,vectorAux); 
		constraintOK=CheclConstraints(butterfly);
		//System.out.println("despues:  "+Arrays.toString(matrixToVector(tempSolution.getMachine_cell()))+" Fit "+tempFitness);
		
			
	}	

		
	@SuppressWarnings("static-access")
	void GenerateMovementButterflyAdjustingOperator(int butterfly,int iteration){
		boolean constraintOK = false;
		tempFitness = 0;
		Vuelo_levy L = new Vuelo_levy();
		double step_levy;
		int[] vectorMaquinaBest = new int[data.M];
		int[] vectorMaquinaActual = new int[data.M];
		int[] vectorAux = new int[data.M];
		vectorMaquinaActual = matrixToVector(tempSolution.getMachine_cell());
		vectorMaquinaBest = matrixToVector(bestSolution.getMachine_cell());	
		while (constraintOK == false) {	
			vectorAux = vectorMaquinaActual;
			do {
				step_levy = L.levy_step(0.1f, 1);				
			} while (Double.isNaN(step_levy));
			
			for (int i = 0; i < data.M; i++) { //para todos los elementos de la mariposa
				double r=Math.random();
				if(r<=p){ //moviemiento 1
			
				vectorAux[i] = vectorMaquinaBest[i]; //sacar de la mejor mariposa
				}
				else{//moviemiento 2				
					int randomInt = ThreadLocalRandom.current().nextInt(NP1, poblation.size());//random de la segunda poblacion
					 vectorAux[i] =GetVectorValue(i,poblation.get(randomInt).getMachine_cell());
					 if(r>BAR){
						 double alpha = SMax/iteration;
					 vectorAux[i] = aproximar(vectorAux[i]+alpha*(step_levy-0.5));
					 }
				}
			}
			
			constraintOK = CheclConstraints(butterfly);
		}
		tempSolution = vectorToMatrix(tempSolution,vectorAux);
		constraintOK=CheclConstraints(butterfly);
	}
	
	void  generateCrossoverOperator(int butterfly){
		int[] vectorMaquinaActual = new int[data.M];
		int[] vectorNewButterfly = new int[data.M];
		int[] vectorMaquinaPreviousMov = new int[data.M]; //sacado del movimiento previo
		auxSolution = tempSolution;
		boolean constraintOK = false;
		float CR=0f;
		float movimiento=0;
		CR = calculateCR(butterfly);
		vectorMaquinaActual = matrixToVector(poblation.get(butterfly).getMachine_cell());
		vectorMaquinaPreviousMov = matrixToVector(tempSolution.getMachine_cell());
		while (constraintOK == false) {	
			CR = calculateCR(butterfly);
			for (int i = 0; i < data.M; i++) { 
				movimiento=vectorMaquinaPreviousMov[i]*(1-CR)+vectorMaquinaActual[i]*CR;
				vectorNewButterfly[i] =  aproximar((float)(movimiento-1));
			}
			tempSolution = vectorToMatrix(tempSolution,vectorNewButterfly);
			constraintOK = CheclConstraints(butterfly);			
		}
		tempSolution2 = tempSolution;
		tempSolution=auxSolution;		
		calculateFitness(tempSolution2);
		constraintOK = CheclConstraints(butterfly);
	}
	
	float calculateCR(int butterfly){
		float CR=0f;
		CR = 0.8f+0.2f*((float)(poblation.get(bestButterflyOnSecondPoblation()).getFitness()-bestSolution.getFitness())/(float)(worstSolution.getFitness()-bestSolution.getFitness()));
		return CR;		
	}
	
	
	boolean CheclConstraints(int butterfly){
		boolean constraintOK = false;
		// Check constraint
					MCDPModel boctorModel = new MCDPModel(data.A, data.M, data.P, data.C, data.mmax,
							tempSolution.getMachine_cell(), tempSolution.getPart_cell());
					constraintOK = boctorModel.checkConstraint();

					if (constraintOK == true) {
						tempFitness = boctorModel.calculateFitness();
						this.numAcceptedMoves++;
						return true;
					} else {
						constraintOK = repararSolucion();
						if (constraintOK == false) {
							tempSolution = poblation.get(butterfly);
						} else {
							tempFitness = boctorModel.calculateFitness();
							return true;
						}
						this.numRejectedMoves++;
					}
		return false;
	}
	
	int[] matrixToVector(int[][] Matrix){
		int[] vector = new int[data.M];
		
		for (int i = 0; i < data.M; i++) {
			for (int j = 0; j < data.C; j++) {
				if (Matrix[i][j] == 1) {
					vector[i] = j + 1;
				}
			}
		}
		return vector;
	}
	
	int GetVectorValue(int machine,int[][] Matrix){
		for (int j = 0; j < data.C; j++) {
			if (Matrix[machine][j] == 1) {
				return  j + 1;
			}
		}
		return (0);
	}
	
	private Solution vectorToMatrix(Solution tempSolution,int [] vectorMaquinaActual){
		int[][] tempMachine_cell = new int[data.M][data.C];
		
		for (int i = 0; i < data.M; i++) {
			tempMachine_cell[i][vectorMaquinaActual[i] - 1] = 1;
		}
		tempSolution.setMachine_cell(tempMachine_cell);
		// Posteriormente generamos manualmente la matriz PxC
		for (int j = 0; j < data.P; j++) {
			for (int k = 0; k < data.C; k++) {
				tempSolution.getPart_cell()[j][k] = 0;
			}
		}
		for (int j = 0; j < data.P; j++) {
			int[] tempPart = new int[data.M];
			int[] cellCount = new int[data.C];

			for (int k = 0; k < data.C; k++) {
				for (int i = 0; i < data.M; i++) {
					tempPart[    i] = tempSolution.getMachine_cell()[i][k] * data.A[i][j];
				}
				cellCount[k] = IntStream.of(tempPart).sum();
			}
			int maxIndex = 0;
			for (int i = 1; i < cellCount.length; i++) {
				int newNumber = cellCount[i];
				if ((newNumber > cellCount[maxIndex])) {
					maxIndex = i;
				}
			}
			tempSolution.getPart_cell()[j][maxIndex] = 1;
		}		
		return tempSolution;
	}	
	
	
	
	
	void calculateFitness(Solution solution){		
		MCDPModel boctorModel = new MCDPModel(data.A, data.M, data.P, data.C, data.mmax,
				tempSolution2.getMachine_cell(), tempSolution2.getPart_cell());
		
		tempSolution2.setFitness(boctorModel.calculateFitness()) ;
	}
	
	int bestButterflyOnSecondPoblation(){
		int fitness=99999;
		int posicion=0;
		for(int i=NP1;i<poblation.size();i++){
			if(poblation.get(i).getFitness()<fitness){
				fitness=poblation.get(i).getFitness();
				posicion=i;
			}
		}
		return posicion;
	}
	
	
public int AutonomousSearchP(int iteraciones){
		String parametros[];
		/*parametros= ParamAutonomous.split("-");
		int CantIntEstan = Integer.parseInt(parametros[0]);
		int CantModoEstan = Integer.parseInt(parametros[1]);
		double step=Double.parseDouble(parametros[2]);*/
		int CantIntEstan = 3;
		int CantModoEstan = 4;
		double step=0.3f;
		double var;
		var=SMax;
		switch (modoDelta){
		case 0: 
			if(iteraciones%CantIntEstan==0&&iteraciones>=CantIntEstan&&var+step<=3.01f){ //aumentar delta de probabilidad
				var = var+step;	
			if(var>varMax)
				varMax=var;
				}
				if(iteraciones>CantModoEstan){ //Si van dos aumentos y aun no mejora
					modoDelta=1;						 //Se cambia al modo de disminuir
					iteraciones=0;
				}
								break;
		case 1: if(iteraciones>=CantIntEstan&&var-step>=0.099f){ //disminuir delta de probabilidad
			var = var-step;		
			if(var<varMin)
				varMin=var;
			}
			if(iteraciones>CantModoEstan){ //Si van dos aumentos y aun no mejora
				modoDelta=0;						 //Se cambia al modo de aumentar
				iteraciones=0;
			}
			
			break;
		}
	//	System.out.println("Iteracion estancada: "+iteraciones+" Delta: "+delta+" modo: ");
		
		SMax=var;
		return iteraciones;
		
	}
	

	public void calcularSimilitudMaquinas() {
		for (int i = 0; i < data.M; i++) {
			for (int j = 0; j < data.M; j++) {
				if (i != j) {
					matrizSimilitud[i][j] = contarSimilitud(i, j);
					matrizSimilitud[j][i] = matrizSimilitud[i][j];
				}
			}
		}
	}

	private int contarSimilitud(int i, int j) {
		int similitud = 0;

		for (int k = 0; k < data.P; k++) {
			if ((data.A[i][k] == 1) && (data.A[i][k] == data.A[j][k])) {
				similitud++;
			}
		}
		return similitud;
	}

	/**
	 * Esta función genera randomicamente las soluciones iniciales. Generar
	 * poblacion aleatorea es EXPLORACION
	 */
	private void generateInitialPoblation() {
		for (int i = 0; i < numberPoblation; i++) {
			// Inicialite procedure
			boolean constraintOK = false;
			// crear una solucion segun los datos leidos
			MCDPRandomSolution randomSolution = new MCDPRandomSolution(data.A, data.M, data.P, data.C, data.mmax);
			int randomSolutionFitness = 0;
			// Estoy en el ciclo hasta generar una solucion randomica que
			// satisfaga las restricciones
			while (constraintOK == false) {
				// Create random solution
				randomSolution.createRandomSolution();
				// Check constraint
				MCDPModel boctorModel = new MCDPModel(data.A, data.M, data.P, data.C, data.mmax,
						randomSolution.getMachine_cell(), randomSolution.getPart_cell());
				constraintOK = boctorModel.checkConstraint();

				if (constraintOK == true) {
					//System.out.println("Paso");
					randomSolutionFitness = boctorModel.calculateFitness();
					this.numAcceptedMoves++;
					break;
				} else {
					//System.out.println("Error");
					this.numRejectedMoves++;
				}
			}

			// Create Solution
			Solution s = new Solution(randomSolution.getMachine_cell(), randomSolution.getPart_cell(),
					randomSolutionFitness);

			// Add Solution in poblation
			poblation.add(s);
		}
	}
	
	private void addRandomSolutionToPoblation(){
		boolean constraintOK = false;
		// crear una solucion segun los datos leidos
		MCDPRandomSolution randomSolution = new MCDPRandomSolution(data.A, data.M, data.P, data.C, data.mmax);
		int randomSolutionFitness = 0;

		// Estoy en el ciclo hasta generar una solucion randomica que
		// satisfaga las restricciones
		while (constraintOK == false) {
			// Create random solution
			randomSolution.createRandomSolution();
			// Check constraint
			MCDPModel boctorModel = new MCDPModel(data.A, data.M, data.P, data.C, data.mmax,
					randomSolution.getMachine_cell(), randomSolution.getPart_cell());
			constraintOK = boctorModel.checkConstraint();

			if (constraintOK == true) {
				//System.out.println("Paso");
				randomSolutionFitness = boctorModel.calculateFitness();
				this.numAcceptedMoves++;
				break;
			} else {
				//System.out.println("Error");
				this.numRejectedMoves++;
			}
		}

		// Create Solution
		Solution s = new Solution(randomSolution.getMachine_cell(), randomSolution.getPart_cell(),
				randomSolutionFitness);

		// Add Solution in poblation
		poblation.add(s);
	}
	
	private void deleteRandomSolutionToPoblation(){
		  Random randomGenerator = new Random();
		 int randomInt = randomGenerator.nextInt(poblation.size());
		
		 do{
			// System.out.println("deleting");
			 if(poblation.get(randomInt)!=bestSolution){
				// System.out.println("Best solution: "+bestSolution+ " solucion eliminada: "+poblation.get(randomInt).getFitness()+" Random: "+randomInt);
				 poblation.remove(randomInt);
				 
			 }
				 
			 randomInt = randomGenerator.nextInt(poblation.size());
		 }while(poblation.get(randomInt)==bestSolution);
	}
	

	@SuppressWarnings("unused")
	private void toConsolePoblation() {
		for (int i = 0; i < numberPoblation; i++) {
			System.out.println(">> Poblation > Solution [" + (i + 1) + "]");
			poblation.get(i).toConsoleMachineCell();
			poblation.get(i).toConsolePartCell();
			poblation.get(i).toConsoleFitness();
			System.out.println("");
		}
	}

	public boolean repararSolucion() {
		boolean constraintOK = false;
		int contAsignaciones = 0, colConMenosAsig = 0, maqMenosAfin = 0, flag = 0, mejorUbicacionCelda;

		MCDPModel boctorModel = new MCDPModel(data.A, data.M, data.P, data.C, data.mmax, tempSolution.getMachine_cell(),
				tempSolution.getPart_cell());

		constraintOK = boctorModel.consistencyConstraint_1();

		if (!constraintOK) {
			for (int i = 0; i < data.M; i++) {
				flag = 0;
				for (int j = 0; j < data.C; j++) {
					if (tempSolution.getMachine_cell()[i][j] == 1) {
						if (flag == 1) {
							mejorUbicacionCelda = buscarCeldaCorrecta(i);
							for (int k = 0; k < data.C; k++) {
								tempSolution.getMachine_cell()[i][k] = 0;
							}
							tempSolution.getMachine_cell()[i][mejorUbicacionCelda] = 1;
							continue;
						} else {
							flag = 1;
						}
					}
				}
			}
		}

		boctorModel = new MCDPModel(data.A, data.M, data.P, data.C, data.mmax, tempSolution.getMachine_cell(),
				tempSolution.getPart_cell());
		constraintOK = boctorModel.consistencyConstraint_3();

		if (!constraintOK) {
			for (int j = 0; j < data.C; j++) {
				contAsignaciones = 0;
				for (int i = 0; i < data.M; i++) {
					if (tempSolution.getMachine_cell()[i][j] == 1) {
						contAsignaciones++;
					}
				}
				if (contAsignaciones > data.mmax) {
					int dif = contAsignaciones - data.mmax;
					while (dif > 0) {
						colConMenosAsig = buscarColumnaMenorAsig(j);
						maqMenosAfin = buscarMaqMenorSim(j);
						for (int k = 0; k < data.C; k++) {
							tempSolution.getMachine_cell()[maqMenosAfin][k] = 0;
						}
						tempSolution.getMachine_cell()[maqMenosAfin][colConMenosAsig] = 1;
						dif--;
					}
				} else {

				}
			}
		}
		for (int j = 0; j < data.P; j++) {
			for (int k = 0; k < data.C; k++) {
				tempSolution.getPart_cell()[j][k] = 0;
			}
		}

		for (int j = 0; j < data.P; j++) {
			int[] tempPart = new int[data.M];
			int[] cellCount = new int[data.C];

			for (int k = 0; k < data.C; k++) {
				for (int i = 0; i < data.M; i++) {
					tempPart[i] = tempSolution.getMachine_cell()[i][k] * data.A[i][j];
				}
				cellCount[k] = IntStream.of(tempPart).sum();
			}
			int maxIndex = 0;
			for (int i = 1; i < cellCount.length; i++) {
				int newNumber = cellCount[i];
				if ((newNumber > cellCount[maxIndex])) {
					maxIndex = i;
				}
			}
			tempSolution.getPart_cell()[j][maxIndex] = 1;
		}

		boctorModel = new MCDPModel(data.A, data.M, data.P, data.C, data.mmax, tempSolution.getMachine_cell(),
				tempSolution.getPart_cell());
		constraintOK = boctorModel.checkConstraint();
		return constraintOK;
	}

	public int buscarMaqMenorSim(int columActual) {
		int contSim = 0, menSimilitud = 1000000000, maquina = 0;

		for (int j = 0; j < data.M; j++) {
			contSim = 0;
			if (tempSolution.getMachine_cell()[j][columActual] == 1) {
				for (int k = 0; k < data.M; k++) {
					contSim = contSim + matrizSimilitud[j][k];
				}
				if (menSimilitud > contSim) {
					menSimilitud = contSim;
					maquina = j;
				}
			}
		}

		return maquina;
	}

	public int buscarColumnaMenorAsig(int columActual) {
		int contAsig = 0, menAsig = 1000000000, colum = 0;
		for (int i = columActual; i < data.C; i++) {
			contAsig = 0;
			for (int j = 0; j < data.M; j++) {
				if (tempSolution.getMachine_cell()[j][i] == 1) {
					contAsig++;
				}
			}
			if (menAsig > contAsig) {
				menAsig = contAsig;
				colum = i;
			}
		}
		return colum;
	}

	public int buscarCeldaCorrecta(int maquina) {
		int similitud = 0, simAnterior = 0, celda = 0;
		for (int j = 0; j < data.C; j++) {
			for (int i = 0; i < maquina; i++) {
				if (tempSolution.getMachine_cell()[i][j] == 1) {
					similitud = similitud + matrizSimilitud[maquina][i];
				}
			}
			if (similitud > simAnterior) {
				simAnterior = similitud;
				celda = j;
			}
		}
		return celda;
	}


	private int aproximar(double movimiento) {
		int aproximado = 0;
		promedio=promedio+movimiento;
		cont++;
		
		//aproximado = IntervalDiscretization.IntervalDoubleValue(SShaped.S5(movimiento), data.C, 0, 1)+1;
		//System.out.println("VALOR MOVIMIENTO: ("+movimiento+") APROXIMADO("+aproximado+")"+"SShaped"+VShaped.V1(movimiento));
		//System.out.println("Aproximado con vshape: "+aproximado);
		
		aproximado = Math.round((float) movimiento);
		//aproximado = (int) VShaped.V1(movimiento);
		if (aproximado < 1) {
			aproximado = 1;
		} else if (aproximado > data.C) {
			aproximado = data.C;
		}
	//System.out.println("movimiento: "+movimiento+"Aproximado con aproximar: "+aproximado);*/
		return aproximado;
	}

	public int binarizacion(double numDiscreto) {
		/*
		 * double random = rn.nextDouble(); if (random <= numDiscreto) { return
		 * 1; } return 0;
		 */
		if (numDiscreto > 0.2) {
			return 1;
		}
		// return Math.round((float) numDiscreto);
		return 0;
	}

	private void chooseBestSolutionInPoblation() {
		// Escoger temporalmente el primer elemento como la mejor solucion.
		double[][] maquina_celda = new double[data.M][data.C];
		double[][] parte_celda = new double[data.P][data.C];

		bestSolution.setMachine_cell(poblation.get(0).getMachine_cell());
		bestSolution.setPart_cell(poblation.get(0).getPart_cell());
		bestSolution.setFitness(poblation.get(0).getFitness());

		for (int i = 0; i < data.M; i++) {
			for (int j = 0; j < data.C; j++) {
				maquina_celda[i][j] = poblation.get(0).getMachine_cell()[i][j];
			}
		}

		for (int i = 0; i < data.P; i++) {
			for (int j = 0; j < data.C; j++) {
				parte_celda[i][j] = poblation.get(0).getPart_cell()[i][j];
			}
		}

		bestSolution.setDoubleMachine_cell(maquina_celda);
		bestSolution.setDoublePart_cell(parte_celda);

		for (int i = 1; i < numberPoblation; i++) {
			if (poblation.get(i).getFitness() < bestSolution.getFitness()) {
				// Escoger una nueva mejor solucion
				bestSolution.setMachine_cell(poblation.get(i).getMachine_cell());
				bestSolution.setPart_cell(poblation.get(i).getPart_cell());
				bestSolution.setFitness(poblation.get(i).getFitness());
				for (int k = 0; k < data.M; k++) {
					for (int j = 0; j < data.C; j++) {
						maquina_celda[k][j] = poblation.get(i).getMachine_cell()[k][j];
					}
				}
				for (int k = 0; k < data.P; k++) {
					for (int j = 0; j < data.C; j++) {
						parte_celda[k][j] = poblation.get(i).getPart_cell()[k][j];
					}
				}
				bestSolution.setDoubleMachine_cell(maquina_celda);
				bestSolution.setDoublePart_cell(parte_celda);
			}
		}
	}
	private void PrintSolutions(int k){
		System.out.println("Set "+k+" de soluciones");
		for(int i=0;i<poblation.size();i++){
			System.out.println(Arrays.toString(matrixToVector(poblation.get(i).getMachine_cell()))+" Fitness " + poblation.get(i).getFitness());
		}
		System.out.println("Mejor Solucion: "+bestSolution.getFitness());
	}


	private void chooseWorstSolutionInPoblation() {
		// Escoger temporalmente el primer elemento como la mejor solucion.
		double[][] maquina_celda = new double[data.M][data.C];
		double[][] parte_celda = new double[data.P][data.C];

		worstSolution.setMachine_cell(poblation.get(0).getMachine_cell());
		worstSolution.setPart_cell(poblation.get(0).getPart_cell());
		worstSolution.setFitness(poblation.get(0).getFitness());

		for (int i = 0; i < data.M; i++) {
			for (int j = 0; j < data.C; j++) {
				maquina_celda[i][j] = poblation.get(0).getMachine_cell()[i][j];
			}
		}

		for (int i = 0; i < data.P; i++) {
			for (int j = 0; j < data.C; j++) {
				parte_celda[i][j] = poblation.get(0).getPart_cell()[i][j];
			}
		}

		worstSolution.setDoubleMachine_cell(maquina_celda);
		worstSolution.setDoublePart_cell(parte_celda);

		for (int i = 1; i < numberPoblation; i++) {
			if (poblation.get(i).getFitness() > worstSolution.getFitness()) {
				// Escoger una nueva mejor solucion
				worstSolution.setMachine_cell(poblation.get(i).getMachine_cell());
				worstSolution.setPart_cell(poblation.get(i).getPart_cell());
				worstSolution.setFitness(poblation.get(i).getFitness());
				for (int k = 0; k < data.M; k++) {
					for (int j = 0; j < data.C; j++) {
						maquina_celda[k][j] = poblation.get(i).getMachine_cell()[k][j];
					}
				}
				for (int k = 0; k < data.P; k++) {
					for (int j = 0; j < data.C; j++) {
						parte_celda[k][j] = poblation.get(i).getPart_cell()[k][j];
					}
				}
				worstSolution.setDoubleMachine_cell(maquina_celda);
				worstSolution.setDoublePart_cell(parte_celda);
			}
		}
	}

	private void toConsoleBestSolution() {
		System.out.println(">> Mejor solucion");
		bestSolution.toConsoleMachineCell();
		bestSolution.toConsolePartCell();
		bestSolution.toConsoleFitness();
	}

	public void toConsoleFinalReport() {
		System.out.println("===============================");
		System.out.println(">> Reporte final");

		toConsoleBestSolution();
		System.out.println(">> Numero de movimientos aceptados: " + this.numAcceptedMoves);
		System.out.println(">> Numero de movimientos fallados : " + this.numRejectedMoves);
	}

	public Solution getBestSolution() {
		return bestSolution;
	}

	@SuppressWarnings("unused")
	private void toConsoleSingleSolutio(int i) {

		System.out.println(i + "Poblacion del vector soluciones");
		poblation.get(i).toConsoleMachineCell();
		poblation.get(i).toConsolePartCell();
		poblation.get(i).toConsoleFitness();
		System.out.println("");

	}
}