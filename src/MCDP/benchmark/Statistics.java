package MCDP.benchmark;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.hssf.util.CellReference;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;

public class Statistics {
	public static void createConvergenciGraph(String identificator, int[] datos, String directoryName,int ejecucion, int ejecuciones) {
		try {

			HSSFWorkbook workbook = new HSSFWorkbook();
			HSSFSheet sheet = workbook.createSheet("Hoja 1");
			
			String directory = directoryName + "/" + identificator + ".xls";
			File file = new File(directory);	
			if(file.exists()){
				 workbook = new HSSFWorkbook(new FileInputStream(directoryName + "/" + identificator + ".xls"));
				 sheet = workbook.getSheet("Hoja 1");
				int fila = 0;
				Row row;
				Cell cell;
				row = sheet.getRow(0);

				cell = row.createCell(ejecucion+1);
				cell.setCellValue("Fitness");

				while (fila < datos.length) {
					row = sheet.getRow(fila+1);					
					cell = row.createCell(ejecucion+1);
					cell.setCellValue(datos[fila]);
					fila++;
				}
			
				if(ejecucion==ejecuciones-1){
					fila = 0;
					row = sheet.getRow(fila);		
					cell = row.createCell(ejecucion+2);
					cell.setCellValue("Promedios");
					while (fila < datos.length) {
						row = sheet.getRow(fila+1);					
						cell = row.createCell(ejecucion+2);
					
						String strFormula= "AVERAGE(B"+(fila+2)+":"+CellReference.convertNumToColString(cell.getColumnIndex()-1)+(fila+2)+")";
						cell.setCellType(HSSFCell.CELL_TYPE_FORMULA);
						cell.setCellFormula(strFormula);
						fila++;
					}
					
					
				}
				
				
				
				
				
				
				FileOutputStream out = new FileOutputStream(new File(directoryName + "/" + identificator + ".xls"));
				workbook.write(out);
				out.close();
				workbook.close();
			}else{
				int fila = 0;
				Row row;
				Cell cell;
				row = sheet.createRow(0);

				cell = row.createCell(0);
				cell.setCellValue("Iteracion");
				cell = row.createCell(1);
				cell.setCellValue("Fitness");

				while (fila < datos.length) {
					row = sheet.createRow(fila + 1);
					cell = row.createCell(0);
					cell.setCellValue(fila + 1);
					cell = row.createCell(1);
					cell.setCellValue(datos[fila]);
					fila++;
				}
				FileOutputStream out = new FileOutputStream(new File(directoryName + "/" + identificator + ".xls"));
				workbook.write(out);
				out.close();
				workbook.close();
				
			}
			
		
			
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void createTable(int fila, int mmax, int bestGlobal, int bestFitness, float meanFitness,String nProblema, String directoryName, int numCell, float iterationOptAvg,long averageTime,double[]variables,int maquinas, int partes) {
		System.out.println(directoryName);
		try {
			String dato[] = nProblema.split("_");
			String num_problema = dato[2].substring(dato[2].length() - 2);
			String nomb_autor = dato[1];
			double RPD = 0;
			if(bestGlobal!=0){
				RPD = ((double)(bestFitness - bestGlobal) / bestGlobal)*100;
			}
			String directory = "Sumary [" + directoryName + "].xls";
			File file = new File(directory);

			if (file.exists()) {
				HSSFWorkbook workbook = new HSSFWorkbook(new FileInputStream("Sumary [" + directoryName + "].xls"));
				HSSFSheet sheet = workbook.getSheet("Hoja 1");

				Row row;
				Cell cell;
				row = sheet.createRow(fila + 1);

				cell = row.createCell(0);
				cell.setCellValue(fila + 1);
				cell = row.createCell(1);
				cell.setCellValue(Integer.parseInt(num_problema)+"-"+nomb_autor);
				cell = row.createCell(2);
				cell.setCellValue(numCell);
				cell = row.createCell(3);
				cell.setCellValue(mmax);
				cell = row.createCell(4);
				cell.setCellValue(bestGlobal);
				cell = row.createCell(5);
				cell.setCellValue(bestFitness);
				cell = row.createCell(6);
				cell.setCellValue(meanFitness);
				cell = row.createCell(7);
				cell.setCellValue(RPD);
				cell = row.createCell(8);
				cell.setCellValue(iterationOptAvg);
				cell = row.createCell(9);
				cell.setCellValue(averageTime);
				cell = row.createCell(10);
				cell.setCellValue(variables[1]);
				cell = row.createCell(11);
				cell.setCellValue(variables[2]);
				cell = row.createCell(12);
				cell.setCellValue(maquinas);
				cell = row.createCell(13);
				cell.setCellValue(partes);

				FileOutputStream out = new FileOutputStream(new File("Sumary [" + directoryName + "].xls"));
				workbook.write(out);
				out.close();

				workbook.close();
			} else {
				HSSFWorkbook workbook = new HSSFWorkbook();
				HSSFSheet sheet = workbook.createSheet("Hoja 1");

				Row row;
				Cell cell;
				row = sheet.createRow(fila);

				cell = row.createCell(0);
				cell.setCellValue("Instance");
				cell = row.createCell(1);
				cell.setCellValue("Autor Problem");
				cell = row.createCell(2);
				cell.setCellValue("num Cell");
				cell = row.createCell(3);
				cell.setCellValue("Mmax");
				cell = row.createCell(4);
				cell.setCellValue("Optimun Value");
				cell = row.createCell(5);
				cell.setCellValue("Best Optimum");
				cell = row.createCell(6);
				cell.setCellValue("Mean Fitness");
				cell = row.createCell(7);
				cell.setCellValue("RPD");
				cell = row.createCell(8);
				cell.setCellValue("Iteration");
				cell = row.createCell(9);
				cell.setCellValue("Time (ms)");
				cell = row.createCell(10);
				cell.setCellValue("varMin");
				cell = row.createCell(11);
				cell.setCellValue("VarMax");
				cell = row.createCell(12);
				cell.setCellValue("maquinas");
				cell = row.createCell(13);
				cell.setCellValue("partes");

				row = sheet.createRow(fila + 1);

				cell = row.createCell(0);
				cell.setCellValue(fila + 1);
				cell = row.createCell(1);
				cell.setCellValue(Integer.parseInt(num_problema)+"-"+nomb_autor);
				cell = row.createCell(2);
				cell.setCellValue(numCell);
				cell = row.createCell(3);
				cell.setCellValue(mmax);
				cell = row.createCell(4);
				cell.setCellValue(bestGlobal);
				cell = row.createCell(5);
				cell.setCellValue(bestFitness);
				cell = row.createCell(6);
				cell.setCellValue(meanFitness);
				cell = row.createCell(7);
				cell.setCellValue(RPD);
				cell = row.createCell(8);
				cell.setCellValue(iterationOptAvg);
				cell = row.createCell(9);
				cell.setCellValue(averageTime);
				cell = row.createCell(10);
				cell.setCellValue(variables[1]);
				cell = row.createCell(11);
				cell.setCellValue(variables[2]);
				cell = row.createCell(12);
				cell.setCellValue(maquinas);
				cell = row.createCell(13);
				cell.setCellValue(partes);
				
				FileOutputStream out = new FileOutputStream(new File("Sumary [" + directoryName + "].xls"));
				workbook.write(out);
				out.close();

				workbook.close();

			}
		} catch (FileNotFoundException e) {

			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public static void writeRunTime(int fila, String directoryName,String var){
		String directory = "Sumary [" + directoryName + "].xls";
		File file = new File(directory);
		try {
		if (file.exists()) {
			HSSFWorkbook workbook = new HSSFWorkbook(new FileInputStream("Sumary [" + directoryName + "].xls"));
			HSSFSheet sheet = workbook.getSheet("Hoja 1");

			Row row;
			Cell cell;
			row = sheet.createRow(fila + 1);
			cell = row.createCell(0);
			cell.setCellValue(var);
			
			
			
			
			
			FileOutputStream out = new FileOutputStream(new File("Sumary [" + directoryName + "].xls"));
			workbook.write(out);
			out.close();

			workbook.close();
		}
		} catch (FileNotFoundException e) {

			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
