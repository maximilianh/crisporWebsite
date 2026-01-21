package ur_rna.StructureEditor.models;

import java.awt.geom.AffineTransform;
import java.util.Collection;
import java.util.Iterator;

/**
 * Represents any grouping (collection) of nucleotides.
 */
public interface INucGroup extends Iterable<Nuc> {
    /**
     * Get a set of all nucleotides included in the group.
     * @return a Collection of nucleotides.
     */
    Collection<Nuc> getBases();

    /**
     * Returns an iterator over all Nucletotides in the group.
     * @return an Iterator.
     */
    default Iterator<Nuc> iterator() {
        return getBases().iterator();
    }

    default void transformAll(AffineTransform tr) {
        for(Nuc n : getBases())
            tr.transform(n.location, n.location);
    }
}
