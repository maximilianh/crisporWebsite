package ur_rna.StructureEditor.models;

import java.util.HashMap;
import java.util.Map;

/**
 *
 */
public class HistoryState {
    public RnaSceneState state;
    public SceneUpdateInfo update;
    public Map<String, Object> stateInfo = new HashMap<>();
    public HistoryState(final RnaSceneState state, final SceneUpdateInfo info) {
        this.state = state;
        this.update = info;
    }
    public void put(String key, Object value) {
        stateInfo.put(key, value);
    }
    @SuppressWarnings("unchecked")
    public <T> T get(String key, T defaultValue) {
        Object value = stateInfo.get(key);
        if (value == null)
            return defaultValue;
        return (T)value;
    }
    public Object get(String key) {
        return stateInfo.get(key);
    }
}
