package MCDP.model;
import java.util.Arrays;

public class Solution
{
	// Solution matrix
	private int[][] machine_cell;
	private int[][] part_cell;	

	private double[][] doubleMachine_cell;
	private double[][] doublePart_cell;
	
	// Fitness solution
	private int fitness;
	
	// Constructor
	public Solution(int[][] machine_cell, int[][] part_cell, int fitness)
	{
		this.machine_cell = new int[machine_cell.length][machine_cell[0].length];
		this.part_cell    = new int[part_cell.length][part_cell[0].length];
		this.doubleMachine_cell =  new double[machine_cell.length][machine_cell[0].length];
		this.doublePart_cell = new double[part_cell.length][part_cell[0].length];
		
        // Copy original values
        for (int i = 0; i < machine_cell.length; i++)
        {
            System.arraycopy(machine_cell[i], 0, this.machine_cell[i], 0, machine_cell[0].length);
        }
        for (int i = 0; i < machine_cell.length; i++){
        	for(int j = 0; j < machine_cell[i].length; j++){
        		this.doubleMachine_cell[i][j]=machine_cell[i][j];        		
        	}
        }
        for (int i = 0; i < part_cell.length; i++)
        {
            System.arraycopy(part_cell[i], 0, this.part_cell[i], 0, part_cell[0].length);
        }
        for (int i = 0; i < part_cell.length; i++)
        {
        	for(int j = 0; j < part_cell[i].length; j++){
        		this.doublePart_cell[i][j]=part_cell[i][j];
        	}
        }
        
		this.fitness = fitness;
	}
	
	public Solution(int[][] machine_cell, int[][] part_cell)
	{	
		this.machine_cell = new int[machine_cell.length][machine_cell[0].length];
		this.part_cell    = new int[part_cell.length][part_cell[0].length];
		
        // Copy original values
        for (int i = 0; i < machine_cell.length; i++)
        {
            System.arraycopy(machine_cell[i], 0, this.machine_cell[i], 0, machine_cell[0].length);
        }

        for (int i = 0; i < part_cell.length; i++)
        {
            System.arraycopy(part_cell[i], 0, this.part_cell[i], 0, part_cell[0].length);
        }
	}
	
	public Solution()
	{
	}
	
	// Get and Set
	public int[][] getMachine_cell() {
		return machine_cell;
	}

	public void setMachine_cell(int[][] machine_cell)
	{
		this.machine_cell = new int[machine_cell.length][machine_cell[0].length];
		
        // Copy original values
        for (int i = 0; i < machine_cell.length; i++)
        {
            System.arraycopy(machine_cell[i], 0, this.machine_cell[i], 0, machine_cell[0].length);
        }
	}

	public int[][] getPart_cell()
	{
		return part_cell;
	}

	public void setPart_cell(int[][] part_cell)
	{
		this.part_cell    = new int[part_cell.length][part_cell[0].length];
		
        for (int i = 0; i < part_cell.length; i++)
        {
            System.arraycopy(part_cell[i], 0, this.part_cell[i], 0, part_cell[0].length);
        }
	}
	
	public int getFitness() {
		return fitness;
	}

	public void setFitness(int fitness) {
		this.fitness = fitness;
	}
	
	public double[][] getDoubleMachine_cell() {
		return doubleMachine_cell;
	}

	public void setDoubleMachine_cell(double[][] doubleMachine_cell) {
		this.doubleMachine_cell = doubleMachine_cell;
	}

	public double[][] getDoublePart_cell() {
		return doublePart_cell;
	}

	public void setDoublePart_cell(double[][] doublePart_cell) {
		this.doublePart_cell = doublePart_cell;
	}

	// To Console
	public void toConsoleMachineCell()
	{
		System.out.println("   > Machine Cell: "+Arrays.deepToString(machine_cell));
	}
	
	public void toConsolePartCell()
	{
		System.out.println("   > Part Cell   : "+Arrays.deepToString(part_cell));
	}
	
	public void toConsoleFitness()
	{
		System.out.println("   > Fitness     : "+fitness);
	}
	
	
}
