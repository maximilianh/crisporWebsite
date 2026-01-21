package ur_rna.Utilities;

/**
 * @author Richard M. Watson
 */
@FunctionalInterface
public interface ActionHandler<T extends EventArgs> {
    void handle(Object source, T e);
}
