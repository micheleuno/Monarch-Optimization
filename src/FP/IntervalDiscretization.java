package FP;



/**
 * Created by Salazar on 5/18/16.
 *
 * Esta clase fue creada para MCDP
 */
public final class IntervalDiscretization
{
    /**
     *
     * @param value value
     * @param numberIntervals number of intervals
     * @param lowerBoundRange lower bound range
     * @param upperBoundRange upper bound range
     * @return un intervalo entre 0 a N-1  (N: Número de intervalos)
     */
    public static int IntervalDoubleValue(double value, int numberIntervals,
                                          double lowerBoundRange, double upperBoundRange)
    {
        // Input validations.

        if (lowerBoundRange > upperBoundRange)
        {
            throw new NumberFormatException("Upper bound range it's should be greater than lower bound range.");
        }

        if (lowerBoundRange == upperBoundRange)
        {
            throw new NumberFormatException("Upper bound range and lower bound range should be different.");
        }

        if (numberIntervals < 0)
        {
            throw new NumberFormatException("Number of intervals must be greater than zero.");
        }

        /*
        if (value > upperBoundRange)
        {
            throw new NumberFormatException("Value it's should be less than upper bound range.");
        }

        if (value < lowerBoundRange)
        {
            throw new NumberFormatException("Value it's should be greater than lower bound range.");
        }*/

        // Calcular el intervalo
        double  rangeNumber = (upperBoundRange - lowerBoundRange) / numberIntervals;
        //System.out.println("range number: "+ rangeNumber);

        int interval;
        double tempRange = lowerBoundRange;
        for (interval = 0; interval < numberIntervals; interval++)
        {
            boolean inInterval = check(value, lowerBoundRange, tempRange + rangeNumber);
            //System.out.println("["+lowerBoundRange+";"+(tempRange + rangeNumber)+"] Numero a discretizar: " +value);
            if (inInterval == true)
            {
                return interval;
            }
            tempRange = tempRange + rangeNumber;
        }
        return interval;
    }

    private static boolean check(double value, double lowerBound, double upperBound)
    {
        if (value >= lowerBound && value <= upperBound)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    public static void main(String[] args)
    {
        System.out.println("resp: " + IntervalDoubleValue(1, 3, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(11, 3, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(110, 3, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(1100, 3, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(11000, 3, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(110000, 3, 0, 1));
        System.out.println("para 3 celdas");
        System.out.println("resp: " + IntervalDoubleValue(0.8898976, 3, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(0.1898976, 3, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(0.0898976, 3, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(0.5, 3, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(0.16, 3, 0, 1));
        System.out.println("para 16 maquinas");
        System.out.println("resp: " + IntervalDoubleValue(0.8898976, 16, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(0.1898976, 16, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(0.0898976, 16, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(0.5, 16, 0, 1));
        System.out.println("resp: " + IntervalDoubleValue(0.1, 16, 0, 1));
    }
}
