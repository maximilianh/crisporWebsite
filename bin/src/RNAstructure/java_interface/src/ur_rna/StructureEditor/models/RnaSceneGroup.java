package ur_rna.StructureEditor.models;

import ur_rna.StructureEditor.FileType;
import ur_rna.StructureEditor.services.fileIO.RnaFileIO;
import ur_rna.Utilities.PathTools;
import ur_rna.Utilities.Strings;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * A group of related RnaScene objects.
 */
public class RnaSceneGroup extends ArrayList<RnaScene> {
    private String filePath;
    private FileType originalType = FileType.AnyOpenable;
    private String title;
    public DrawSettings settings;
    public String getFilePath() {
        return filePath;
    }
    public FileType getFileType() { return originalType; }
    public String getTitle() {
        return title;
    }
    public void setTitle(String title) {
        this.title = title;
    } //structure.SetSequenceLabel(title)
    public int structureCount() { return this.size(); }
    public void setSource(final String path, final FileType type) {
        filePath = path;
        originalType = type;
    }
    public static RnaSceneGroup from(final RnaScene... scenes) {
        RnaSceneGroup group = new RnaSceneGroup();
        group.addAll(Arrays.asList(scenes));
        group.guessTitle();
        return group;
    }
    public void guessTitle() {
        if (!Strings.isEmpty(title)) return;
        for(RnaScene r : this)
            if (!Strings.isEmpty(r.title)) {
                title =  RnaFileIO.stripEnergyTitle(r.title);
                return;
            }
        if (!Strings.isEmpty(filePath))
            title = PathTools.getBaseName(filePath);
    }
    public RnaSceneGroup subset(final RnaScene... scenes) {
        RnaSceneGroup copy = this.copyProps();
        copy.addAll(Arrays.asList(scenes));
        copy.guessTitle();
        return copy;
    }
    public RnaSceneGroup subset(final int... scenes) {
        RnaSceneGroup copy = this.copyProps();
        for (int i = 0; i < scenes.length; i++) {
            copy.add(this.get(i));
        }
        copy.guessTitle();
        return copy;
    }
    private RnaSceneGroup copyProps() {
        RnaSceneGroup copy = new RnaSceneGroup();
        copy.filePath = this.filePath;
        copy.originalType = this.originalType;
        copy.title = this.title;
        copy.settings = this.settings;
        return copy;
    }
}
