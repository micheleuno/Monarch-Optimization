package MCDP.model;
import java.util.Arrays;

public class MCDPModel
{
	private int[][] A;
	private int M;
	private int P;
	private int C;
	private int Mmax;
	
	// Solution matrix
	private int[][] machine_cell;
	private int[][] part_cell;
	
	public MCDPModel(int A[][], int M, int P, int C, int Mmax, int[][] machine_cell, int[][] part_cell)
	{
		this.A = A;
		this.M = M;
		this.P = P;
		this.C = C;
		this.Mmax = Mmax;
		this.machine_cell = machine_cell;
		this.part_cell = part_cell;
	}
	


	/**
	 * 
	 * @return true si esta correta la restriccion, en caso contrario retorna falso.
	 */
	public boolean consistencyConstraint_1()
	{
		for (int i = 0; i < this.M; i++)
		{
			int flag = 0;
			int zeros = 0;
			for (int k = 0; k < this.C; k++)
			{
				// Check that there is only one true value in the column.
				if (machine_cell[i][k] == 1)
				{
					if (flag == 1)
					{
						return false;
					}
					else
					{
						flag = 1;
					}
				}
				else
				{
					zeros = zeros + 1;
				}
			}
			if (zeros == C)
			{
				return false;
			}
			else
			{
				zeros = 0;
			}
		}
		
		return true;
	}
	
	/**
	 * 
	 * @return true si esta correcta la restriccion, en caso contrario retorna falso.
	 */	
	public boolean consistencyConstraint_2()
	{
		for (int j = 0; j < this.P; j++)
		{
			int flag = 0;
			int zeros = 0;
			for (int k = 0; k < this.C; k++)
			{
				// Check that there is only one true value in the column.
				if (part_cell[j][k] == 1)
				{
					if (flag == 1)
					{
						return false;
					}
					else
					{
						flag = 1;
					}
				}
				else
				{
					zeros = zeros + 1;
				}
			}
			
			if (zeros == this.C)
			{
				return false;
			}
			else
			{
				zeros = 0;
			}
		}
		
		return true;
	}
	/**
	 * 
	 * @return true si esta correta la restriccion, en caso contrario retorna falso.
	 */	
	public boolean consistencyConstraint_3()
	{
		int count = 0;
		
		for (int k = 0; k < this.C; k++)
		{
			for (int i = 0; i < this.M; i++)
			{
				// Count 1 in the Column
				if (machine_cell[i][k] == 1)
				{
					count++;
				}
			}
			
			// Check Summation <= Mmax
			if (count > this.Mmax)
			{
				return false;
			}
			
			count = 0;
		}
		
		return true;
	}
	
	public boolean checkConstraint()
	{
		boolean c1 = consistencyConstraint_1();
		boolean c2 = consistencyConstraint_2();
		boolean c3 = consistencyConstraint_3();
		
		if ((c1 == false) )
		{
			return false;
		}
		
		if ( (c2 == false) )
		{
			return false;
		}
		
		if ( (c3 == false))
		{
			return false;
		}
		
		return true;
	}
	
	// Objective function
	public int calculateFitness()
	{
		int sum = 0;
		
		for (int k = 0; k < this.C; k++)
		{
			for (int i = 0; i < this.M; i++)
			{
				for (int j = 0; j < this.P; j++)
				{
					sum = sum + (A[i][j] * part_cell[j][k]) * (1 - machine_cell[i][k]);			
				}
			}
		}
		return sum;
	}
	
	/**
	 *  Esta parte son funciones utilitarias
	 */
	public void convertToFinalMatrix()
	{
		int[][] A_Solution_Temp = new int[M][P];
		int[][] finalSolution_A = new int[M][P];
		// Sort Machine Cell
		int ii = 0;
		int jj = 0;
		
		// Ordenar Maquinas
		for (int k = 0; k < C; k++)
		{
			for (int i = 0; i < M; i++)
			{
				// Si Y[maquina][celda] == 1
				if (machine_cell[i][k] == 1)
				{
					// copiar toda la fila de la maquina Y[i][k] en la matriz A2 (final)
					for (int j = 0; j < P; j++)
					{
						A_Solution_Temp[ii][j] = A[i][j];
					}
					ii++;
				}
			}
		}
		
		// Ordenar Partes.
		for (int k = 0; k < C; k++)
		{
			for (int j = 0; j < P; j++)
			{
				// Si Z[parte][celda] == 1
				if (part_cell[j][k] == 1)
				{
					// copiar toda la fila de la parte Z[j][k] en la matriz A_Solution (final),
					// acá se realiza la copia utilizando la matriz A_Solution.
					for (int i = 0; i < M; i++)
					{
						finalSolution_A[i][jj] = A_Solution_Temp[i][j];
					}
					jj++;
				}
			}
		}
		
		System.out.println(">> Matriz A de incidencias re-ordenada en celdas:");
		for (int i = 0; i < M; i++)
		{
			System.out.println(Arrays.toString(finalSolution_A[i]));
		}
	}

}
