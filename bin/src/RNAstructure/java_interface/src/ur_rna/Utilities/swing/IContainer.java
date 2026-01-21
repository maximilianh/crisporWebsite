package ur_rna.Utilities.swing;

import java.awt.*;

/**
 *
 */
public interface IContainer {
    Component[] getComponents();
    int getComponentCount();
    Component getComponent(int index);


    static Component[] getComponents(Object o) {
        if (o instanceof IContainer) return ((IContainer)o).getComponents();
        if (o instanceof Container) return ((Container)o).getComponents();
        throw new UnsupportedOperationException();
    }
    static int getComponentCount(Object o) {
        if (o instanceof IContainer) return ((IContainer)o).getComponentCount();
        if (o instanceof Container) return ((Container)o).getComponentCount();
        throw new UnsupportedOperationException();
    }
    static Component getComponent(Object o, int index) {
        if (o instanceof IContainer) return ((IContainer)o).getComponent(index);
        if (o instanceof Container) return ((Container)o).getComponent(index);
        throw new UnsupportedOperationException();
    }
}
