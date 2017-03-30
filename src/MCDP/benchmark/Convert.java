package MCDP.benchmark;


public class Convert
{

	public Convert()
	{
		
	}
	
	public int[][] StringToMatrixInt(String string, int rows, int columns)
	{
		String delims = " ";
		int[][] matrix = new int[rows][columns];

		// Cut String an put into matrix.
		String[] tokens = string.split(delims);
		int tokenCount = tokens.length;
		
		if (tokenCount != (rows * columns))
		{
			System.out.println("Error in matrix size");
			return null;
		}
		
		int i = 0;
		for (int j = 0; j < rows; j++)
		{
			for (int k = 0; k < columns; k++){
				matrix[j][k] = Integer.parseInt(tokens[i]);
				i++;
			}
		}			

		return matrix;
	}
	
	public boolean[][] StringToMatrixBoolean(String string, int rows, int columns)
	{
		String delims = " ";
		boolean[][] matrix = new boolean[rows][columns];

		// Cut String an put into matrix.
		String[] tokens = string.split(delims);
		int tokenCount = tokens.length;
		
		if (tokenCount != (rows * columns))
		{
			System.out.println("Error in matrix size");
			return null;
		}
		
		int i = 0;
		for (int j = 0; j < rows; j++)
		{
			for (int k = 0; k < columns; k++){
				matrix[j][k] = Boolean.parseBoolean(tokens[i]);
				i++;
			}
		}			

		return matrix;
	}

}
