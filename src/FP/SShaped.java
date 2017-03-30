package FP;



/**
 * Created by Salazar on 5/18/16.
 */
public final class SShaped
{
    // S1
    public static double S1(float x)
    {
        return 1/(1 + Math.pow(Math.E, (-2 * x)));
    }

    public static double S1(double x)
    {
        return 1/(1 + Math.pow(Math.E, (-2 * x)));
    }

    // S2
    public static double S2(float x)
    {
        return 1/(1 + Math.pow(Math.E, (-1 * x)));
    }

    public static double S2(double x)
    {
        return 1/(1 + Math.pow(Math.E, (-1 * x)));
    }

    // S3
    public static double S3(float x)
    {
        return 1/(1 + Math.pow(Math.E, ((-1 * x)/2) ));
    }

    public static double S3(double x)
    {
        return 1/(1 + Math.pow(Math.E, ((-1 * x)/2) ));
    }

    // S4
    public static double S4(float x)
    {
        return 1/(1 + Math.pow(Math.E, ((-1 * x)/3) ));
    }

    public static double S4(double x)
    {
        return 1/(1 + Math.pow(Math.E, ((-1 * x)/3) ));
    }

    public static void main(String[] args)
    {
        TCATest();
    }
    //S5
    public static double S5(float x)
    {
        return 1/(1 + Math.pow(Math.E, (-2 * (x-2.015))));
    }
    
    public static double S5(double x)
    {
        return 1/(1 + Math.pow(Math.E, (-2 * (x-2.015))));
    }
    // S6
    public static double S6(float x)
    {
        return 1/(1 + Math.pow(Math.E, (-1 * x-2.015)));
    }

    public static double S6(double x)
    {
        return 1/(1 + Math.pow(Math.E, (-1 * x-2.015)));
    }

    // S7
    public static double S7(float x)
    {
        return 1/(1 + Math.pow(Math.E, ((-1 * x-2.015)/2) ));
    }

    public static double S7(double x)
    {
        return 1/(1 + Math.pow(Math.E, ((-1 * x-2.015)/2) ));
    }

    // S8
    public static double S8(float x)
    {
        return 1/(1 + Math.pow(Math.E, ((-1 * x-2.015)/3) ));
    }

    public static double S8(double x)
    {
        return 1/(1 + Math.pow(Math.E, ((-1 * x-2.015)/3) ));
    }

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    public static void TCATest()
    {
        final double[] a = new double[]
        {
            0.001,
            0.01
        };


    }

    public static void basicTest()
    {
        System.out.println("S1(-1): " + S1(-1));
        System.out.println("S1(0): " + S1(0));
        System.out.println("S1(1): " + S1(1));
        System.out.println("S1(-0.999999): " + S1(-0.999999));
        System.out.println("S1(0.999999): " + S1(0.999999));
        System.out.println("S1(-0.111111): " + S1(-0.111111));
        System.out.println("S1(0.111111): " + S1(0.111111));
        System.out.println("S1(-1.6939411510425213): " + S1(-1.6939411510425213));
        System.out.println("S1(1.6939411510425213): " + S1(1.6939411510425213));

        System.out.println("---------------");

        System.out.println("S4(-1): " + S4(-1));
        System.out.println("S4(0): " + S4(0));
        System.out.println("S4(1): " + S4(1));
        System.out.println("S4(-0.999999): " + S4(-0.999999));
        System.out.println("S4(0.999999): " + S4(0.999999));
        System.out.println("S4(-0.111111): " + S4(-0.111111));
        System.out.println("S4(0.111111): " + S4(0.111111));
        System.out.println("S4(-1.6939411510425213): " + S4(-1.6939411510425213));
        System.out.println("S4(1.6939411510425213): " + S4(1.6939411510425213));

        System.out.println("---------------");
        System.out.println("S4(0.8): " + S4(0.8));
        System.out.println("S4(-0.1): " + S4(-0.1));
        System.out.println("S4(0): " + S4(0));
        System.out.println("S4(1.1): " + S4(1.1));
        System.out.println("S4(0.5): " + S4(0.5));
    }
}
