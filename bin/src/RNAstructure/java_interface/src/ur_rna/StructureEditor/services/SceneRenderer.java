package ur_rna.StructureEditor.services;

import ur_rna.StructureEditor.models.DrawSettings;
import ur_rna.StructureEditor.models.IScreenObject;
import ur_rna.StructureEditor.services.drawing.View2D;

import java.awt.*;
import java.util.Collection;

public interface SceneRenderer {
    Collection<IScreenObject> render(Graphics2D graphics, View2D view, boolean includeInteractive);
    DrawSettings getSettings();
    Rectangle calcBounds(Graphics2D g, View2D view, boolean includeInteractive);
}
