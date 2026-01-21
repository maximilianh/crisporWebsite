package ur_rna.StructureEditor.services.drawing;

import java.awt.*;

/**
 * Represents a swing component (e.g. JPanel etc) on which an RNA scene can be drawn.
 */
public interface ICanvas {
    View2D getView();
    Component getComponent();
    void repaintRequired();
}
