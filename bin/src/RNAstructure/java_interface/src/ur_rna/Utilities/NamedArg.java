package ur_rna.Utilities;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Richard M. Watson
 */
public class NamedArg {
    public final String name;
    public final Object value;
    public NamedArg(String name, Object value) {
        this.name = name;
        this.value = value;
    }
    public static Map<String,Object> toMap(NamedArg[] args, boolean includeUnnamed) {
        return toMap(args, true, "#");
    }
    public static Map<String,Object> toMap(NamedArg[] args, boolean includeUnnamed, String prefixUnnamed) {
        Map<String,Object> map = new HashMap<>(args.length);
        for (int i = 0; i < args.length; i++)
            if (args[i].hasName())
                map.put(args[i].name, args[i].value);
            else if (includeUnnamed)
                map.put(prefixUnnamed + i, args[i].value);
        return map;
    }
    public boolean hasName() { return !Strings.isEmpty(name); }
    public int intValue() { return Convert.toInt(value); }
    public float floatValue() { return Convert.toFloat(value); }
    public double doubleValue() { return Convert.toDouble(value); }
    public char charValue() { return Convert.toChar(value); }
    public boolean boolValue() { return Convert.toBool(value); }
    public String stringValue() { return Convert.toString(value); }
}
