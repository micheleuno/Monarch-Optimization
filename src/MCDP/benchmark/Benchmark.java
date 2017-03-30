package MCDP.benchmark;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

import org.apache.log4j.Logger;

import MCDP.model.MCDPData;


public class Benchmark
{
	/**
	 * 
	 */
	private ArrayList<String> fileNameSet;
	
	/**
	 * Other parameters
	 */
	Logger log;
	
	public Benchmark()
	{
		fileNameSet = new ArrayList<String>();
		
		// Inicializate log.
		log = Logger.getLogger(Benchmark.class);
	}
	
	// leo a traves de un listado de nombres de archivo, los modelos con datos.
	public ArrayList<MCDPData> getSetModelsByNames(ArrayList<String> filenames) throws IOException
	{
		ArrayList<MCDPData> models = new ArrayList<MCDPData>();
		Iterator<String> nombreIterator = filenames.iterator();

		while (nombreIterator.hasNext())
		{
			String filename = nombreIterator.next();
			//System.out.println(filename);
			Hashtable<String,String> hashTable = new Hashtable<String,String>();
			
			PropertiesRead p = new PropertiesRead();
			hashTable = p.getPropEssentialValues(filename);

			Convert convert = new Convert();			

			// Get Values for MCDP Model.
			int M = Integer.parseInt(hashTable.get("Machines"), 10);
			int P = Integer.parseInt(hashTable.get("Parts"), 10);
			int C = Integer.parseInt(hashTable.get("Cells"), 10);
			int Mmax = Integer.parseInt(hashTable.get("Mmax"), 10);
			@SuppressWarnings("unused")
			int Sum = Integer.parseInt(hashTable.get("Sum"), 10);
			String string = hashTable.get("Matrix").toString();
			int A[][] = convert.StringToMatrixInt(string, M, P);		
			// Set values into model.
			MCDPData model = new MCDPData(M, P, C, Mmax, A,filename);

			// Set identificator Model (KEY).
			int pos = filename.lastIndexOf(".");
			if (pos > 0)
			{
				filename = filename.substring(0, pos);
			}
			//model.setIdentificator(filename);

			// Add model into model set.
			models.add(model);	
			//filename = new String();
		}

		return models;
	}
	
	public ArrayList<String> readSetFileBenchmark(String file) throws Exception
	{
		try
		{
			InputStream in = getClass().getResourceAsStream(file); 
			BufferedReader br = new BufferedReader(new InputStreamReader(in));	
			
			String str;
			
			while ((str = br.readLine()) != null)
			{
				log.info(str);
				fileNameSet.add(str);
			}
			br.close();	
		}
		catch (IOException io)
		{
			throw new Exception(io);
		} 
		catch (Exception e)
		{
			log.error("", e);
		}
		
		setFileNameSet(fileNameSet);
		
		return fileNameSet;
	}
	
	public ArrayList<String> getFileNameSet()
	{
		return fileNameSet;
	}

	public void setFileNameSet(ArrayList<String> fileNameSet)
	{
		this.fileNameSet = fileNameSet;
	}
}
