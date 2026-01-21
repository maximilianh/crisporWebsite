package ur_rna.RNAstructureUI.utilities;

/**
 * @author Richard M. Watson
 */
public abstract class RunFunc implements Runnable {
    public RunFunc(Object param) { this.param = param; }
    public RunFunc() { this.param = null; }
    protected final Object param;
    protected Object result;
    @Override
    public void run() {
        this.result = run(this.param);
    }
    public abstract Object run(Object param);
    public Object getResult() { return result; }
    public boolean getBoolResult() { return result instanceof Boolean ? (Boolean)result : false; }
    public int getIntResult() { return (Integer)result; }
    public String getStrResult() { return result instanceof String ? (String)result : result == null ? null : result.toString(); }
}
