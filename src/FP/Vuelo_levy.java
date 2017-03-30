package FP;


public class Vuelo_levy {
	static double sigma2;
	static double L;
	static double U;
	static double V;
	static double s;

	public static double levy_step(double delta,double varsBusqueda){
		sigma2=Math.pow(((gamma(1.0+delta)/(delta*gamma((1.0+delta)/2.0)))*((Math.sin(Math.PI*delta/2.0))/(Math.pow(2.0,((delta-1)/2.0))))),1.0/delta) ;   
		V=gaussrandom(0.0,1.0);
		U=gaussrandom(0.0,sigma2);
		s=U/Math.pow(Math.abs(V),1.0/delta);
		L=((delta*gamma(delta)*Math.sin(Math.PI*delta/2.0))/Math.PI)*(1.0/(Math.pow(s,1.0+delta))); 
		return L;
	}

	public static double  gaussrandom(double xort,double r){
		int iset=0;
		double gset=0;
		double fac,rsq,v1,v2;
		double aret=0;
		if(iset==0){
			do{
				v1=2.0*Math.random()-1.0;
				v2=2.0*Math.random()-1.0;
				rsq=v1*v1+v2*v2;  
			}while(rsq>=1.0 || rsq==0.0);	   
			fac=Math.sqrt(-2.0*Math.log(rsq)/rsq);
			gset=v1*fac;
			iset=1;
			aret=v2*fac;
		}else{
			iset=0;
			aret=gset;
		}
		return aret/8.0*r+xort;	
	}

	static double logGamma(double x){
		double tmp =(x-0.5)*Math.log(x+4.5)-(x+4.5);
		double ser = 1.0 + 76.18009173    / (x + 0)   - 86.50532033    / (x + 1)+ 24.01409822    / (x + 2)   -  1.231739516   / (x + 3)+  0.00120858003 / (x + 4)   -  0.00000536382 / (x + 5);
		return tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
	}

	static double gamma(double x) {
		return Math.exp(logGamma(x));
	}

}
