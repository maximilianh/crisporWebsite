package ur_rna.StructureEditor.services;

import ur_rna.StructureEditor.models.RnaScene;
import ur_rna.StructureEditor.models.SceneUpdateInfo;
import ur_rna.Utilities.annotation.Nullable;

/**
 * Event that indicates an update has occurred in the specified RnaScene.
  */
public class SceneUpdateEvent {
    public final SceneUpdateCategory type;
    public final RnaScene scene;
    public final Object source;
    @Nullable
    public final String description;
    public SceneUpdateEvent(final Object source, final SceneUpdateInfo info, final RnaScene scene) {
        this.source = source;
        this.type = info.category;
        this.scene = scene;
        this.description = info.description;
    }
}
