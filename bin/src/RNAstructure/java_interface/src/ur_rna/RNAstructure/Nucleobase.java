package ur_rna.RNAstructure;

/**
 * Nucleobase symbols and corresponding numeric values (as per the RNA backend)
 */
public class Nucleobase {
    public final static int Unknown = 0;
    public final static int A = 1, C = 2, G = 3, T = 4;
    /** Represents and intermolecular break between separate strands ("I") */
    public final static int Break = 5;

    public final static char[] rnaSymbols = new char[] { 'X' ,'A', 'C', 'G', 'U', 'I' };
    public final static char[] rnaSymbolsLc = new char[] { 'x' ,'a', 'c', 'g', 'u', 'i' };
    public final static char[] dnaSymbols = new char[] { 'X' ,'A', 'C', 'G', 'T', 'I' };
    public final static char[] dnaSymbolsLc = new char[] { 'x' ,'a', 'c', 'g', 't', 'i' };
    public final static char[] allowedSymbols = new char[] { 'X' ,'A', 'C', 'G', 'U', 'T', 'I', 'x' ,'a', 'c', 'g', 'u', 't', 'i' };

    public final static String getSymbol(int backendValue, boolean isRNA, boolean upper) {
        if (backendValue < 0)
            backendValue = 0;
        char[] symbols = isRNA ? (upper ? rnaSymbols : rnaSymbolsLc) : (upper ? dnaSymbols : dnaSymbolsLc);
        if (backendValue >= symbols.length)
            backendValue = 0;
        return Character.toString(symbols[backendValue]);
    }
    public final static String getSymbol(int backendValue, boolean isRNA) { return getSymbol(backendValue, isRNA, true); }
    public final static String getSymbol(int backendValue) { return getSymbol(backendValue, true, true);   }

    public static char getVerifiedSymbol(String symbol) { return getVerifiedSymbol(symbol, 'X'); }
        public static char getVerifiedSymbol(String symbol, char defaultIfInvalid) {
        if (symbol != null && symbol.length() == 1) {
            char c = symbol.charAt(0);
            for (char allowedSymbol : allowedSymbols)
                if (allowedSymbol == c)
                    return c;
        }
        return defaultIfInvalid;
    }


//    if (base[0]=='A'||base[0]=='a') numseq[count]=1;
//    else if (base[0]=='C'||base[0]=='c') numseq[count]=2;
//    else if (base[0]=='G'||base[0]=='g') numseq[count]=3;
//    else if (base[0]=='U'||base[0]=='u'||base[0]=='T'||base[0]=='t') numseq[count]=4;
//    else if (base[0]=='I') numseq[count]=5;
//    else numseq[count]=0;
        }
