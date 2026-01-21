package ur_rna.StructureEditor.models;

import ur_rna.StructureEditor.services.SceneUpdateCategory;

/**
 * Represents the type of update that was applied to the RnaScene (e.g. to cause a change to the undo/redo history)
 */
public class SceneUpdateInfo {
    // The following updates cause changes to the scene.
    public static SceneUpdateInfo Initial = new SceneUpdateInfo(SceneUpdateCategory.Structure, "Initial State"),
    SequenceSize = new SceneUpdateInfo(SceneUpdateCategory.Structure, "Added/Removed Bases"), // Add/remove nucleotides. Changes the size of the structure.
    StrandDivisions = new SceneUpdateInfo(SceneUpdateCategory.Structure, "Joined/Divided Strands"), // Divide a strand into two strands  or join two strands into one.
    SequenceText = new SceneUpdateInfo(SceneUpdateCategory.Structure, "Modified Sequence Text"),  // Change nucleotide symbols etc (but not the number of nucs)
    BondType = new SceneUpdateInfo(SceneUpdateCategory.Structure, "Changed Type of Bonds"), // Change Bond Type (Normal, Pseudo etc)
    Bonds = new SceneUpdateInfo(SceneUpdateCategory.Structure, "Added/Removed Basepairs"), // Add/Remove Bonds
            AddBonds=Bonds.subType("Added Basepairs"),
            RemBonds=Bonds.subType("Removed Basepairs"),
    Layout = new SceneUpdateInfo(SceneUpdateCategory.Layout, "Moved Bases"),
            Loop = Layout.subType("Resized Loop"),
            Rotate = Layout.subType("Rotated Bases"),
            Flip = Layout.subType("Flipped Bases"),
    FormatBases = new SceneUpdateInfo(SceneUpdateCategory.Style, "Changed Color/Style of Bases"),
    GotoHistory = new SceneUpdateInfo(SceneUpdateCategory.Structure, "Undo/Redo caused Full Reload");

    // The following do not cause changes to History, but do require redrawing and possibly change to UI elements.
    public static SceneUpdateInfo Selection = new SceneUpdateInfo(SceneUpdateCategory.Selection, "Selection Changed"),
        EditHistory = new SceneUpdateInfo(SceneUpdateCategory.History, "History Entries Modified"), // History was cleared or a new entry was added that changed history. This may affect UI controls related to history, but the scene itself is unmodified.
        View = new SceneUpdateInfo(SceneUpdateCategory.View, "Changed View (Zoom etc)"), // change the view (by Zoom etc)
        Controls = new SceneUpdateInfo(SceneUpdateCategory.Controls, "UI Controls Updated") // UI elements changed in reaction to user mouse movements (or other interaction). The scene should be redrawn, even though no change occurred to the scene itself.

    ;
//    public SceneUpdateCategory category() { return category; }
//    public String description() { return description; }
//    public SceneUpdateInfo parent() { return parent; }
    public final SceneUpdateCategory category;
    public final String description;
    public final SceneUpdateInfo parent;

    public SceneUpdateInfo(final SceneUpdateCategory category, String description) { this(category, description, null); }
    protected SceneUpdateInfo(final SceneUpdateCategory category, String description, SceneUpdateInfo parentType) {
        parent = parentType==null?this:parentType;
        this.category = category;
        this.description = description;
    }
    public SceneUpdateInfo subType(String description) {
        return new SceneUpdateInfo(this.category, description, this.parent);
    }
}
