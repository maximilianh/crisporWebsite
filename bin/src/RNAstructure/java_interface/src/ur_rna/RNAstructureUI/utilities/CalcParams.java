package ur_rna.RNAstructureUI.utilities;

import ur_rna.Utilities.ObjTools;

import java.util.HashMap;

/**
 * @author Richard M. Watson
 */
public class CalcParams extends HashMap<String, Object> {
    protected boolean verifyType = true;
    public void setVerify(boolean enableVerification) { verifyType = enableVerification; }
    public boolean getVerify() { return verifyType; }

    public String getStr(String paramName) { return getStr(paramName, ""); }
    public String getStr(String paramName, String defaultValue) {
        return ObjTools.toStr(super.getOrDefault(paramName, defaultValue));
    }
    public int getInt(String paramName) { return getInt(paramName, 0); }
    public int getInt(String paramName, int defaultValue) {
        Object value = super.get(paramName);
        doVerify(paramName, value, Number.class);
        return value instanceof Number ? ((Number) value).intValue() : defaultValue;
    }
    public double getDbl(String paramName) { return getDbl(paramName, 0d); }
    public double getDbl(String paramName, double defaultValue) {
        Object value = super.get(paramName);
        doVerify(paramName, value, Number.class);
        return value instanceof Number ? ((Number) value).doubleValue() : defaultValue;
    }
    public boolean getBool(String paramName) { return getBool(paramName, false); }
    public boolean getBool(String paramName, boolean defaultValue) {
        return ObjTools.asBool(super.getOrDefault(paramName, defaultValue));
    }

    private void doVerify(String name, final Object value, final Class<?> type) {
        if (!verifyType || value == null || type.isInstance(value)) return;
        throw new IllegalArgumentException(String.format(
                "The calculation parameter '%s' requires a %s value, but a %s was given instead.",
                name, type.getSimpleName(), value.getClass().getSimpleName()));
    }
}
