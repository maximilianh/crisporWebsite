package ur_rna.StructureEditor.services;

import ur_rna.StructureEditor.models.Bond;
import ur_rna.StructureEditor.models.Nuc;

import java.awt.*;
import java.util.Collection;
import java.util.Map;

/**
 * Applies color to nucleobases in a scene.
 */
public abstract class SceneColorizer {
    protected int colorMode; // one of the ColorMode constants

    protected SceneColorizer(int colorMode) {
        this.colorMode = colorMode;
    }
    public static SceneColorizer ColorRemover(final int colorMode) { return SingleColor(null, colorMode); }
    public static SceneColorizer SingleColor(final Color color, final int colorMode) {
        return new SceneColorizer(colorMode) {
            @Override
            public boolean color(Nuc n) { return setColor(n, color);}
        };
    }
    public static SceneColorizer ColorBySymbol(final Map<String,Color> colors, final int colorMode) {
        return new SceneColorizer(colorMode) {
            @Override
            public boolean color(Nuc n) {
                Color c = colors.get(n.symbol);
                return c != null && setColor(n, c);
            }
        };
    }
    public static SceneColorizer ColorByDataValue(final double[] values, final NumberRangeList ranges, final int colorMode) {
        return new SceneColorizer(colorMode) {
            @Override
            public boolean color(Nuc n) {
                int index = n.indexInScene();
                if (index >= values.length || Double.isNaN(values[index])) return false;
                double d = values[index];
                Color c = ranges.lookup(d);
                return c != null && setColor(n, c);
            }
        };
    }

    public static SceneColorizer ColorByProbability(final PairProbabilityData values, final NumberRangeList ranges, final int colorMode) {
        return new SceneColorizer(colorMode) {
            @Override
            public boolean color(Nuc n) {
                double prob;
                Bond b = n.getPairBond();
                if (b == null)
                    prob = values.getSingle(n.indexInScene());
                else
                    prob = values.getPair(b.n5.indexInScene(), b.n3.indexInScene());

                Color c = ranges.lookup(prob);
                return c != null && setColor(n, c);
            }
        };
    }

    /** Return true if any nucs were affected. */
    public boolean color(Collection<Nuc> nucs) {
        boolean changed = false;
        for(Nuc n : nucs)
            changed |= color(n);
        return changed;
    }
    public abstract boolean color(Nuc n);

    protected boolean setColor(Nuc n, Color c) { return setColor(n, c, colorMode); }
    protected boolean setColor(final Nuc n, Color c, final int colorMode) {
        // Note: If c == null and n.style == null colorsEqual returns false (and so this function returns false).
        // Therefore either c != null or n.style != null
        // If the former is true, n.style already exists so calling style() won't change anything.
        // If the latter is true, c is not null, so calling n.style() is necessary to create the style and set the color to non-null c.
        switch (colorMode) {
            case ColorMode.None:
                return false;
            case ColorMode.Fill:
                if (!colorsEqual(n.style==null?null:n.style.fillColor, c)) { n.style().fillColor = c; return true; }
                break;
            case ColorMode.Text:
                if (!colorsEqual(n.style==null?null:n.style.textColor, c)) { n.style().textColor = c; return true; }
                break;
            case ColorMode.Outline:
                if (!colorsEqual(n.style==null?null:n.style.outlineColor, c)) { n.style().outlineColor = c; return true; }
                break;
            case ColorMode.Bond:
                if (!colorsEqual(n.style==null?null:n.style.bondColor, c)) { n.style().bondColor = c; return true; }
                break;
            default:
                // could be a combination of colors. Set all that apply.
                // note that "|" is used instead of "||" this is important because we want to run all cases, which would not be done if short-circuit evaluation were used.
                return (0!=(colorMode&ColorMode.Fill) && setColor(n, c, ColorMode.Fill))
                        | (0!=(colorMode&ColorMode.Text) && setColor(n, c, ColorMode.Text))
                        | (0!=(colorMode&ColorMode.Outline) && setColor(n, c, ColorMode.Outline))
                        | (0!=(colorMode&ColorMode.Bond) && setColor(n, c, ColorMode.Bond));
        }
        return false;
    }

    private static boolean colorsEqual(Color c1, Color c2) {
        return c1==c2 || c1 != null && c1.equals(c2);
    }

    public static class NumberRangeList {
        public double[] limits;
        public Color[] colors;
        public Comparison comparison;
        public DoubleComparer cmp;
        public Color defaultColor;

        public NumberRangeList(){}
        public NumberRangeList(final double[] limits, final Color[] colors, final String comparisonSymbol, final Color defaultColor) { this(limits, colors, Comparison.fromSymbol(comparisonSymbol), defaultColor);}
        public NumberRangeList(final double[] limits, final Color[] colors, final Comparison cmp, final Color defaultColor) {
            this.limits = limits;
            this.colors = colors;
            this.comparison = cmp;
            this.cmp = comparison.getComparer();
            this.defaultColor = defaultColor;
        }
        public Color lookup(final double d) {
            for (int i = 0; i < limits.length; i++)
                if (cmp.compare(d, limits[i]))
                    return colors[i];

//            switch (cmp) {
//                case Greater:
//                    for (int i = 0; i < limits.length; i++)
//                        if (d > limits[i])
//                            return colors[i];
//                    break;
//                case Less:
//                    for (int i = 0; i < limits.length; i++)
//                        if (d < limits[i])
//                            return colors[i];
//                    break;
//                case GreaterOrEqual:
//                    for (int i = 0; i < limits.length; i++)
//                        if (d >= limits[i])
//                            return colors[i];
//                    break;
//                case LessOrEqual:
//                    for (int i = 0; i < limits.length; i++)
//                        if (d <= limits[i])
//                            return colors[i];
//                    break;
//            }
            return defaultColor;
        }

        // Sort the rules in case the user hasn't listed them correctly.
        // For < and <= comparisons, sort the limits in ascending order (from least to greatest) so that
        //   when a data value is analyzed, it will match the smallest limit that it is less than (or equal to).
        // Conversely, for > and >= comparisons, sort the values in descending order so that
        //   a data value will only match the largest value that is is greater than (or equal to).
        public void sortForComparison() {
            boolean sortAscending = comparison == Comparison.Less | comparison == Comparison.LessOrEqual;
            // Uses insertion sort. Simple, and probably runs in linear time in most cases.
            for (int i = 1; i < limits.length; i++) {
                double d = limits[i];
                Color c = colors[i];
                int j = i - 1;
                while (j >= 0 && (sortAscending == limits[j] > d)) {
                    limits[j + 1] = limits[j];
                    colors[j + 1] = colors[j];
                    j = j - 1;
                }
                limits[j + 1] = d;
                colors[j + 1] = c;
            }
        }
    }

    // Data is listed in the double[] in upper-triangular order. For N=5 (number of nucleobases) the array looks like this:
    // P(1:2) P(1:3) P(1:4) P(1:5)
    // P(2:3) P(2:4) P(2:5)
    // P(3:4) P(3:5)
    // P(4:5)
    // (i.e. 5*4/2 = 10 elements)
    // So to calculate the position of P(I, J)  (where I < J -- if not, swap I, J)
    // Position of first element in row I (0-based) is:
    // ROW(I) = TRIANGLE(N) - TRIANGLE(N-I)   (where TRIANGLE(X) is the number of elements in triangular array -- i.e. X*(X-1)/2
    //        = (N-1)N/2    - (N-I-1)(N-I)/2
    //        = (N-1)N/2    - (N-I-1)(N-I)/2
    //        = (  N^2 - N  -  (N^2-2NI+I^2 +I -N)  ) / 2
    //        = (  2NI - I^2 - I )/2
    //        = (2N-I-1)I/2
    // The first element in row I is P(I, I+1) so the position of P(I, J) is ROW(I) + J - (I+1)
    // so P(I, J) (where I < J) is at (2*N-I-1)*I/2+J-I-1    (for 0-based I and J)
    //                                (2*N-I)*(I-1)/2+J-I-1  (for 1-based I and J -- as is the case with user input)
    public static class PairProbabilityData {
        private int N;
        private double[] values;
        public PairProbabilityData(int nucCount) {
            N = nucCount;
            values = new double[nucCount * (nucCount-1) / 2];
        }
        /**
         * Get the probability of pairing between Nucs I and J.
         * The indices are assumed to be 0-based.
         */
        public double getPair(int i, int j) {
            if (i > j) { int t = i; i = j; j = t; } //swap i and j if i > j
            return values[(2*N-i-1)*i/2+j-i-1];
        }
        /**
         * Set the probability of pairing between Nucs I and J.
         * The indices are assumed to be 0-based.
         */
        public void setPair(int i, int j, double value) {
            if (i > j) { int t = i; i = j; j = t; } //swap i and j if i > j
            values[(2*N-i-1)*i/2+j-i-1] = value;
        }
        /**
         * Get the probability that Nuc I is single-stranded.
         * This is the same as 1.0 - SUM(getPair(I, J)) for all J != I
         */
        public double getSingle(int i) {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += getPair(i, j);
            for (int j = N-1; j > i; j--)
                sum += getPair(j, i);
            return 1 - sum;
        }
    }

    private interface DoubleComparer {
        boolean compare(double lhs, double rhs);
    }
    public enum Comparison {
        Less("<", (a, b)->a < b),
        Greater(">", (a, b)->a > b),
        LessOrEqual("<=", (a, b)->a <= b),
        GreaterOrEqual(">=", (a, b)->a >= b);
        public final String symbol;
        public final DoubleComparer comparer;
        Comparison(String textSymbol, DoubleComparer cmp) { symbol=textSymbol; comparer = cmp; }
        private static Comparison[] values = Comparison.values();
        public static Comparison fromSymbol(String symbol) {
            for(Comparison c : values)
                if (c.symbol.equals(symbol)) return c;
            return null;
        }
        public boolean compare(double lhs, double rhs) { return comparer.compare(lhs, rhs); }
        public DoubleComparer getComparer() { return comparer; }
    }

    public interface ColorMode {
        int None = 0,
            Fill = 1,
            Text = 2,
            Outline = 4,
            Bond = 8;
    }
}
