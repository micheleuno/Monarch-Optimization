package MCDP.benchmark;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Hashtable;
import java.util.Properties;

public class PropertiesRead
{
	
	public PropertiesRead()
	{
		// TODO Auto-generated constructor stub
	}
	
	// Get Essential Values for MCDP
	public Hashtable<String,String> getPropEssentialValues(String fileName) throws IOException
	{		 
		Properties prop = new Properties();
		String propFileName = "/resources/"+fileName;
	
		Hashtable<String,String> hashTable = new Hashtable<String,String>();
 
		InputStream inputStream = getClass().getResourceAsStream(propFileName);
 
		if (inputStream != null)
		{
			prop.load(inputStream);
		}
		else
		{
			throw new FileNotFoundException("property file '" + propFileName + "' not found in the classpath");
		}
 
		// get the property value and print it out
		String machines = prop.getProperty("Machines");
		String parts = prop.getProperty("Parts");
		String cells = prop.getProperty("Cells");
		String mmax = prop.getProperty("Mmax");
		String sum = prop.getProperty("Sum");
		String matrix = prop.getProperty("Matrix");
		
		//System.out.println(result + "\nMachines "+machines+" Parts "+parts+" Matrix "+matrix);
		
		hashTable.put("Machines", machines);
		hashTable.put("Parts", parts);
		hashTable.put("Cells", cells);
		hashTable.put("Mmax", mmax);
		hashTable.put("Sum", sum);
		hashTable.put("Matrix", matrix);
			
		return hashTable;
	}

	// Get File details values
	public Hashtable<String,String> getPropFileDetailsValues(String fileName) throws IOException
	{	
		Properties prop = new Properties();
		String propFileName = fileName;
		Hashtable<String,String> hashTable = new Hashtable<String,String>();
 
		InputStream inputStream = getClass().getClassLoader().getResourceAsStream(propFileName);
 
		if (inputStream != null)
		{
			prop.load(inputStream);
		}
		else
		{
			throw new FileNotFoundException("property file '" + propFileName + "' not found in the classpath");
		}

		// get the property value and print it out
		String fileDetails = prop.getProperty("File details");
		String version = prop.getProperty("Version");
		String createdBy = prop.getProperty("Created by");
		String dateCreated = prop.getProperty("Date created");
		String modifiedBy = prop.getProperty("Modified by");
		String dateModified = prop.getProperty("Date modified");
		
		hashTable.put("FileDetails", fileDetails);
		hashTable.put("Version", version);
		hashTable.put("CreatedBy", createdBy);
		hashTable.put("DateCreated", dateCreated);
		hashTable.put("ModifiedBy", modifiedBy);
		hashTable.put("DateModified", dateModified);
			
		return hashTable;
	}
	
	// Get Benchmark Info
	public Hashtable<String,String> getPropBenchmarkInfoValues(String fileName) throws IOException
	{	
		Properties prop = new Properties();
		String propFileName = fileName;
		Hashtable<String,String> hashTable = new Hashtable<String,String>();
 
		InputStream inputStream = getClass().getClassLoader().getResourceAsStream(propFileName);
 
		if (inputStream != null)
		{
			prop.load(inputStream);
		}
		else
		{
			throw new FileNotFoundException("property file '" + propFileName + "' not found in the classpath");
		}

		// get the property value and print it out
		String benchmarkInfo = prop.getProperty("Benchmark info");
		String references = prop.getProperty("References");
		String bestSolution = prop.getProperty("Best Solution");
		
		hashTable.put("BenchmarkInfo", benchmarkInfo);
		hashTable.put("References", references);
		hashTable.put("BestSolution", bestSolution);
			
		return hashTable;
	}
	
	// Get all values
	public Hashtable<String,String> getPropAllValues(String fileName) throws IOException
	{	
		Properties prop = new Properties();
		String propFileName = fileName;
		Hashtable<String,String> hashTable = new Hashtable<String,String>();
 
		InputStream inputStream = getClass().getClassLoader().getResourceAsStream(propFileName);
 
		if (inputStream != null)
		{
			prop.load(inputStream);
		}
		else
		{
			throw new FileNotFoundException("property file '" + propFileName + "' not found in the classpath");
		}
 
		// get the property value and print it out
		String fileDetails = prop.getProperty("File details");
		String version = prop.getProperty("Version");
		String createdBy = prop.getProperty("Created by");
		String dateCreated = prop.getProperty("Date created");
		String modifiedBy = prop.getProperty("Modified by");
		String dateModified = prop.getProperty("Date modified");
		
		// get the property value and print it out
		String machines = prop.getProperty("Machines");
		String parts = prop.getProperty("Parts");
		String cells = prop.getProperty("Cells");
		String mmax = prop.getProperty("Mmax");
		String sum = prop.getProperty("Sum");
		String matrix = prop.getProperty("Matrix");
		
		// get the property value and print it out
		String benchmarkInfo = prop.getProperty("Benchmark info");
		String references = prop.getProperty("References");
		String bestSolution = prop.getProperty("Best Solution");
		
		hashTable.put("FileDetails", fileDetails);
		hashTable.put("Version", version);
		hashTable.put("CreatedBy", createdBy);
		hashTable.put("DateCreated", dateCreated);
		hashTable.put("ModifiedBy", modifiedBy);
		hashTable.put("DateModified", dateModified);
		
		hashTable.put("BenchmarkInfo", benchmarkInfo);
		hashTable.put("References", references);
		hashTable.put("BestSolution", bestSolution);
		
		hashTable.put("Machines", machines);
		hashTable.put("Parts", parts);
		hashTable.put("Cells", cells);
		hashTable.put("Mmax", mmax);
		hashTable.put("Sum", sum);
		hashTable.put("Matrix", matrix);
			
		return hashTable;
	}

}
