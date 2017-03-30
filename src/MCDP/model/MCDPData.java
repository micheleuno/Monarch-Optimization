package MCDP.model;

public class MCDPData
{
	private String identificator;
	private int bestSGlobal;

	public int M ;
	public int P ;
	public int C ;
	public int mmax ;
	public int [][] A ;
	
	
	public MCDPData( int M, int P, int C, int Mmax, int[][] A, String identificator) {
		this.M=M;
		this.P=P;
		this.C=C;
		this.mmax=Mmax;
		this.A=A;
		this.identificator=identificator;
	}
	public int getM() {
		return M;
	}
	public void setM(int m) {
		M = m;
	}
	public int getP() {
		return P;
	}
	public void setP(int p) {
		P = p;
	}
	public int getC() {
		return C;
	}
	public void setC(int c) {
		C = c;
	}
	public int getMmax() {
		return mmax;
	}
	public void setMmax(int mmax) {
		this.mmax = mmax;
	}
	public int[][] getA() {
		return A;
	}
	public void setA(int[][] a) {
		A = a;
	}
	
	public String getIdentificator() {
		return identificator;
	}
	public void setIdentificator(String identificator) {
		this.identificator = identificator;
	}
	public int getBestSGlobal() {
		return bestSGlobal;
	}
	public void setBestSGlobal(int bestSGlobal) {
		this.bestSGlobal = bestSGlobal;
	}

}
