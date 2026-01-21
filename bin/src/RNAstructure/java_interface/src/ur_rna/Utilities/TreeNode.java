package ur_rna.Utilities;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Represents an abstract Tree data structure.
 * Extend TreeNode to introduce custom properties.
 * e.g. {@code class MyTree extends TreeNode<MyTree> { ... } }
 */
public abstract class TreeNode<T extends TreeNode> {
    List<T> children;
    public List<T> getChidren() {
        if (children == null)
            return Collections.emptyList();
        return children;
    }
    public boolean hasChildren() {
        return children != null && children.size() != 0;
    }
    public void add(T child) {
        if (children == null)
            children = new ArrayList<>();
    }
}