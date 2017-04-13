import java.awt.		Toolkit;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.Locale;

import org.apache.log4j.Logger;

import FP.MonarchOptimization;
import FP.SendMail;
import MCDP.benchmark.Benchmark;
import MCDP.benchmark.Statistics;
import MCDP.model.MCDPData;
import MCDP.model.Solution;

public class Main {

	public static void main(String[] args) throws Exception {
		// Crear parametros iniciales de la metaheuristica
		int numberPoblation = 50;
		int numberIteration = 100;
		float SMax = 1.0f; //Max step
		float BAR = 0.416f;// butterfly adjusting rate
		float peri = 1.2f; //migration period
		float p = 0.416f; // migration ratio 
		//parametros ejecucion
		int executions = 31;
		int best_fitness = 999999999;
		float mean_fitness = 0f;
		int optimal_global = 0;
		int numIteracion = 0;
		double[] variables = new double[3];
		float iterationOptAvg = 0f;
		long ejecutionTimeAvg = 0;
		long tiempoInicio,tiempoFin = 0;
		String ParamAutonomous = null;
		Logger log = Logger.getLogger(Main.class);
		String finLectura="Stop";
		String parametros[];
		File f = new File("Parametros.txt");
		BufferedReader  br= null;
		 br = new BufferedReader(new InputStreamReader( new FileInputStream(f)));
		while(finLectura!=null){//hasta que no haya nuevos parametros
			numIteracion=0;
		if(f.exists() && !f.isDirectory()) { 
			parametros = (br.readLine().split("="));
			numberIteration=Integer.parseInt(parametros[1]);
			parametros = (br.readLine().split("="));
			numberPoblation=Integer.parseInt(parametros[1]);
			parametros = (br.readLine().split("="));
			SMax=Float.parseFloat(parametros[1]);
			parametros = (br.readLine().split("="));
			BAR=Float.parseFloat(parametros[1]);
			parametros = (br.readLine().split("="));
			peri=Float.parseFloat(parametros[1]);
			parametros = (br.readLine().split("="));
			p=Float.parseFloat(parametros[1]);
			parametros = (br.readLine().split("="));
			executions=Integer.parseInt(parametros[1]);
			br.readLine();
			ParamAutonomous=br.readLine();
			finLectura=br.readLine();
		}

		log.info("Read all filenames");
		Benchmark benchmark = new Benchmark();
		ArrayList<String> dataFiles = benchmark.readSetFileBenchmark("/resources/MCDP_BENCHMARK_FILES.txt");
		
		ArrayList<MCDPData> modelSet = benchmark.getSetModelsByNames(dataFiles);
		System.out.println("Read all filenames");
		Iterator<MCDPData> iterator = modelSet.iterator();

		String timeStamp = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss").format(new Date());
		String benchmarkFileConfig = "n_" + numberPoblation + "_SM_" + SMax + "_Bar_" + BAR + "_peri_"+ peri + "_p_"+ p + "_it_"
				+ numberIteration;

		String currentDirectory = timeStamp + "_" + benchmarkFileConfig+" "+ParamAutonomous;
		String folder = "Results";
		String directory = folder + "/" + currentDirectory;

		File directoryFile = new File(folder);
		boolean b = directoryFile.mkdir();
		if (b) {
			System.out.printf("Successfully created new folder: %s%n", folder);
		} else {
			System.out.printf("Failed to create new folder: ", folder);
			if (directoryFile.exists()) {
				System.out.printf("folder " + folder + " already exists\n");
			}
		}

		directoryFile = new File(directory);
		b = directoryFile.mkdir();
		if (b) {
			System.out.printf("Successfully created new directory: %s%n\n", directory);
		} else {
			System.out.printf("Failed to create new directory: %s%n\n", directory);
		}	
		long inicioejecucion = System.currentTimeMillis();
		while (iterator.hasNext()) {

			MCDPData model = iterator.next();

			optimal_global = obtenerOptimo("src/resources/" + model.getIdentificator());
			model.setBestSGlobal(optimal_global);

			long startBenchmark = System.currentTimeMillis();
			SimpleDateFormat sdf = new SimpleDateFormat("MMM dd,yyyy HH:mm:ss.SSS", Locale.ENGLISH);
			Date resultdate = new Date(startBenchmark);
			System.out.println("Process Start, with " + executions + " executions: " + sdf.format(resultdate));
			ejecutionTimeAvg=0;
			for (int i = 0; i < executions; i++) {
			
				MonarchOptimization metaheuristic = new MonarchOptimization(numberPoblation, numberIteration, model, SMax,
						BAR,peri,p, directory);
				tiempoInicio=System.currentTimeMillis();
				variables = metaheuristic.run(ParamAutonomous,i,executions);
				tiempoFin=System.currentTimeMillis();
				
				iterationOptAvg = (float) (variables[0] + iterationOptAvg);
				Solution bestSolution = metaheuristic.getBestSolution();
				mean_fitness = mean_fitness + bestSolution.getFitness();
				if (best_fitness > bestSolution.getFitness()) {
					best_fitness = bestSolution.getFitness();
				}
				ejecutionTimeAvg=ejecutionTimeAvg+(tiempoFin-tiempoInicio);
			
			}
			mean_fitness = mean_fitness / executions;
			iterationOptAvg = iterationOptAvg / executions;
			ejecutionTimeAvg = ejecutionTimeAvg / executions;
			System.out.println("Mean Best Fitness:[" + mean_fitness + "] " + "Best Solution: [" + best_fitness + "]");
			Statistics.createTable(numIteracion, model.mmax, model.getBestSGlobal(), best_fitness, mean_fitness,model.getIdentificator(), currentDirectory,model.getC(),iterationOptAvg,ejecutionTimeAvg,variables,model.M,model.P);
			mean_fitness = 0;
			
			numIteracion++;

			System.out.println("Problem [" + model.getIdentificator() + "]");
			double RPD = ((double)(optimal_global - best_fitness) / optimal_global)*100;
			System.out.println("Best solution global " + optimal_global+" RPD: "+RPD);
			best_fitness = 999999999;
			long endBenchmark = System.currentTimeMillis();
			Date resultdate2 = new Date(endBenchmark);
			System.out.println("Process End: " + sdf.format(resultdate2));
			System.out.println("==========================================================");
		}
		long finejecucion = System.currentTimeMillis();
		System.out.println("Run time= "+((float)(finejecucion-inicioejecucion)/60000)+" minutos");
		Statistics.writeRunTime(numIteracion, currentDirectory, "Run time= "+((float)(finejecucion-inicioejecucion)/60000)+" minutos");
		}
	
		br.close();
		Toolkit.getDefaultToolkit().beep();
	//	SendMail.main(args);
	}

	public static int obtenerOptimo(String directorio) throws IOException {
		String dato[];
		String line32 = Files.readAllLines(Paths.get(directorio)).get(21);
		dato = line32.split("=");
		return Integer.parseInt(dato[1]);
	}

}
