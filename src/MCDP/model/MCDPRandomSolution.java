package MCDP.model;
import java.util.Random;
import java.util.stream.IntStream;

public class MCDPRandomSolution
{
	// Solution matrix
	private int[][] machine_cell;
	private int[][] part_cell;
	
	// Data
	private int[][] A;
	private int M;
	private int P;
	private int C;
	@SuppressWarnings("unused")
	private int Mmax;
	
	// Constructor
	public MCDPRandomSolution(int A[][], int M, int P, int C, int Mmax)
	{
		this.A = A;
		this.M = M;
		this.P = P;
		this.C = C;
		this.Mmax = Mmax;
	}
	
	public void createRandomSolution()
	{
		// create random permutation of X e machine_cell
	    this.machine_cell = new int[M][C];
	    this.part_cell = new int[P][C];
	
	    // Generar matripart_cell M*C randomicamente
	    for (int i = 0; i < M; i++)
	    {
	        Random random = new Random();
	        int cell;
	
	        // The arramachine_cell is initialipart_celled with part_celleros.
	        for (int k = 0; k < C; k++)
	        {
	            machine_cell[i][k] = 0;
	        }
	
	        // Get random value between 1 to Number Cell.
	        cell = (int) (random.nextDouble() * C);
	        machine_cell[i][cell] = 1;
	    }

	    //Se crea la matripart_cell piepart_cella×celda para el momento en que se necesita
	    for (int j = 0; j < P; j++) //Rellenar la matripart_cell piepart_cella×celda de "0"
	    {
	        for (int k = 0; k < C; k++)
	        {
	            part_cell[j][k] = 0;
	        }
	    }

	    for (int j = 0; j < P; j++)
	    {	
	        int[] tempPart = new int[M];
	        int[] cellCount = new int[C];
	
	        for (int k = 0; k < C; k++)
	        {
	            for (int i = 0; i < M; i++)
	            {
	                // Esto hace una multiplicación para revisar si: P(j) es subconjunto de C(k)
	                tempPart[i] = machine_cell[i][k] * A[i][j];
	            }
	            cellCount[k] = IntStream.of(tempPart).sum();
	        }
	        
	        // Extraer el índice de la posición con el número más grande.
	        int maxIndex = 0;
	        for (int i = 1; i < cellCount.length; i++)
	        {
	            int newNumber = cellCount[i];
	            if ((newNumber > cellCount[maxIndex]))
	            {
	                maxIndex = i;
	            }
	        }
	        // Puts un 1 in the new cell
	        part_cell[j][maxIndex] = 1;
	    }
	}
	
	// Get and Set
	public int[][] getMachine_cell() {
		return machine_cell;
	}

	public void setMachine_cell(int[][] machine_cell) {
		this.machine_cell = machine_cell;
	}

	public int[][] getPart_cell() {
		return part_cell;
	}

	public void setPart_cell(int[][] part_cell) {
		this.part_cell = part_cell;
	}

}
