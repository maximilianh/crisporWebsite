package ur_rna.StructureEditor.services;

/**
 * Various types of changes that can occur to a Scene as a result of user interaction.
 */
public enum SceneUpdateCategory {
    /** Changes to the current nucleotide selection */
    Selection,
    /** Changes to Sequence or Bonds */
    Structure,
    /** Changes to Nucleotide styles (color etc) */
    Style,
    /** Changes to Nucleotide positions */
    Layout,
    /** Changes to Zoom, Offset, or ViewPort */
    View,
    /** Changes to the DrawHandles or other controls that need to be redrawn. */
    Controls,
    /** Changes to the Undo/Redo history */
    History
}
