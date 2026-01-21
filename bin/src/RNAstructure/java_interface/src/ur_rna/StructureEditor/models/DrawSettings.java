package ur_rna.StructureEditor.models;

import ur_rna.Utilities.Colors;

import java.awt.*;

/**
 * Container for settings related to drawing RNA scenes.
 *
 * The fields in this class are saved and loaded along with program settings,
 * but fields marked as transient are NOT saved or loaded.
 *
 */
public final class DrawSettings implements Cloneable {
    // **********************************************
    // Settings that affect the Scene as a whole.
    // **********************************************
    public Color sceneBgColor = Color.WHITE;
    /**
     * Spacing between nucs in a helix.
     * Does not include radius, which must be added in to calculate full point-to-point distance.
     * See {@link #calcOptimalNuc2NucDistance(Nuc, Nuc)} for details.
     */
    public float nucSpacing = 1.6f;
    public float nucRadius = 8f;
    /** Spacing between nucs in a loop or hairpin.
     * Does not include radius, which must be added in to calculate full point-to-point distance.
     * See {@link #calcOptimalLoopDistance(Nuc, Nuc, int)} for details.
     */
    public float loopSpacing = 4f;
    public double circularStartAngle = Math.PI /2; // start position when drawing circular. PI/2 is 90 degrees -- 12 o'clock

    // **********************************************
    // Settings that affect nucleotides and are overridden by NucStyles
    // **********************************************
    public Color nucFillColor = Colors.LightGray;
    public Color nucTextColor = Colors.Black;
    public Color nucOutlineColor = Colors.Black;
    public Font nucFont = new Font("Arial", Font.BOLD, 12);

    private transient NucStyle defaultStyle; // created based on the above.
    public float nucOutlineWidth = 1; // not currently overridable per-nuc

    // **********************************************
    // Settings that affect bonds
    // **********************************************
    public Color bondColor = Color.BLUE;
    public Color bondColorPseudo = Color.PINK;
    public float bondWidth = 2;
    public float pseudoBondWidth = 2.5f;
    public boolean drawBondsInHelix = true;
    public transient float bondLength; // calculated based on nuc radius.

    // **********************************************
    // Settings that affect the backbone
    // **********************************************
    public Color backboneColor = Color.GRAY;
    public float backboneWidth = 1;
    public boolean drawBackboneInLoops = true;
    public boolean drawBackboneInHelix = true;
    public boolean drawBackboneCurved = true;

    // **********************************************
    // Settings that affect interactive controls:  selection, bond editing, DrawHints etc.
    // **********************************************
    public Color nucSelColor = Color.BLUE;
    public Color nucFocusColor = Color.GREEN;
    public Color nucHitColor = new Color(0xFF004E);
    public Color selBandLineColor = new Color(0x2D95FF);
    public Color selBandFillColor = new Color(0x332D95FF, true);

    public Color hint = new Color(0x23B3FF);
    public Color hintEmphasis = new Color(0x28519D);

    public transient Stroke nucOutlineSelStroke = new BasicStroke(2);
    public transient Stroke editBondStroke = new BasicStroke(3, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10f, new float[] {4}, 0);
    public transient Stroke selBandOutline = new BasicStroke(2); //new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10f, new float[] { 5 }, 0f);

    // **********************************************
    // Settings that affect nucleotide numbering
    // **********************************************
    public Font numberFont = new Font("Arial", Font.BOLD, 14);
    public Color numberLineColor = Color.BLACK;
    public float numberLineMargin = 2;
    public float numberDistance = 30f;
    public float numberLineWidth = 1;
    public boolean drawNumbers = true;
    public int drawNumbersInterval = 10;
    public boolean drawNumbersAtHelix = false;

    // **********************************************
    // Strokes that are created based on *Width settings.
    // **********************************************
    public transient Stroke bondStroke;
    public transient Stroke bondStrokePseudo;
    public transient Stroke backboneStroke;
    public transient Stroke nucOutlineStroke;
    public transient Stroke numberLineStroke;

    public DrawSettings() {
        onSettingsChanged();
    }

    public NucStyle getStyle(Nuc n) {
        NucStyle s = new NucStyle();
        fillStyle(n, s);
        return s;
    }

    public void fillStyle(Nuc n, NucStyle fillProperties) {
        NucStyle s = n == null ? null : n.style;
        if (s == null)
            defaultStyle.copyTo(fillProperties);
        else {
            fillProperties.outlineColor = s.outlineColor == null ? defaultStyle.outlineColor : s.outlineColor;
            fillProperties.fillColor = s.fillColor == null ? defaultStyle.fillColor : s.fillColor;
            fillProperties.textColor = s.textColor == null ? defaultStyle.textColor : s.textColor;
            fillProperties.bondColor = s.bondColor == null ? defaultStyle.bondColor : s.bondColor;
            fillProperties.font = s.font == null ? defaultStyle.font : s.font;
        }
    }

    /**
     * Call this to update dependent settings whenever settings have been changed.
     * E.g. bondStroke needs to be re-created whenever bondWidth is changed.
     */
    public void onSettingsChanged() {
        bondLength = calcBondLength();

        bondStroke = new BasicStroke(bondWidth);
        bondStrokePseudo = new BasicStroke(pseudoBondWidth);
        backboneStroke = new BasicStroke(backboneWidth);
        nucOutlineStroke = new BasicStroke(nucOutlineWidth);
        numberLineStroke = new BasicStroke(numberLineWidth);

        defaultStyle = new NucStyle(nucFillColor, nucTextColor, nucOutlineColor, bondColor, nucFont, Float.NaN);
    }

    public float nucDiameter() { return 2*nucRadius; }
    //public float nucBackendDiameter() { return 2*nucRadius+nucSpacing; }
    public float calcBondLength() { return 0.9f*2*nucRadius; }

    /**
     * Calculate the optimal distance between nuc midpoints.
     * Does not include any extra space on the outside ends of the nucs.
     * e.g. For 2 adjacent nucs:      N1)-(N2   (two radii + 1 spacer)
     * e.g. For 3 consecutive nucs:   N1)-(N2)-(N3   (4 radii + 2 spacers)
     *
     * Let N be ABS(n1.index-n2.index), then distance = N*(2*radius+spacer)
     *
     * See {@link #calcOptimalOuterDistance(Nuc, Nuc, int)} for a calculation
     * that includes the end-to-end distance when the outer ends are included.
     */
    public float calcOptimalNuc2NucDistance(Nuc n1, Nuc n2) {
        return calcOptimalNuc2NucDistance(Math.abs(n2.indexInScene - n1.indexInScene)+1);
    }
    public float calcOptimalNuc2NucDistance(int nucCount) {
        return (nucCount-1)*(2 * nucRadius + nucSpacing);
    }

     /**
     * Calculate the optimal distance from the outside of one nuc to the outside of another.
     * This always includes the outer radii and can optionally include 0, 1 or 2 outside spacers.
     * e.g. For 2 adjacent nucs, endSpacers==0:        (N1)-(N2)     (4 radii + 1 total spacer)
     * e.g. For 2 adjacent nucs, endSpacers==1:       -(N1)-(N2)     (4 radii + 2 total spacers)
     * e.g. For 2 adjacent nucs, endSpacers==2:       -(N1)-(N2)-    (4 radii + 3 total spacers)
     * e.g. For 3 consecutive nucs, endSpacers==0:   (N1)-(N2)-(N3)  (6 radii + 2 total spacers)
     * e.g. For 3 consecutive nucs, endSpacers==1:  -(N1)-(N2)-(N3)  (6 radii + 3 total spacers)
     * e.g. For 3 consecutive nucs, endSpacers==2:  -(N1)-(N2)-(N3)- (6 radii + 4 total spacers)
     * */
    public float calcOptimalOuterDistance(Nuc n1, Nuc n2, int endSpacers) {
        return calcOptimalOuterDistance(Math.abs(n2.indexInScene-n1.indexInScene)+1, endSpacers);
    }
    public float calcOptimalOuterDistance(int nucCount, int endSpacers) {
        return nucCount*2*nucRadius + (nucCount+endSpacers-1)*nucSpacing;
    }

    /**
     * Calculate the optimal circumference for a circle containing nucs in the given region.
     * @param n1 the nucleotide at the start of the region (inclusive)
     * @param n2 the nucleotide at the end of the region (inclusive)
     * @param endSpacers the number of spacers required on the outside of the region.
     *                   0 means no space before the center of n1 or the after the center of n2
     *                   (e.g. the arc-length of a hairpin loop from the center of n1 to the center of n2)
     *                   1 means a single spacer (e.g. connecting n1 to n2 in a circle)
     *                   2 means a spacer before n1 and after n2 (e.g. a spacer connecting n1 to a previous nuc and another connecting n2 to a subsequent nuc)
     */
    public float calcOptimalLoopDistance(Nuc n1, Nuc n2, int endSpacers) {
        return calcOptimalLoopDistance(Math.abs(n2.indexInScene - n1.indexInScene)+1, endSpacers);
    }
    /**
     * Calculate the optimal circumference for a circle containing the given number of nucleotides.
     * @param nucCount the number of nucleotides in the circle
     * @param endSpacers the number of spacers required on the outside of the terminal nucs (n1=first and n2=last)
     *                   0 means no space before the center of n1 or the after the center of n2
     *                   (e.g. the arc-length of a hairpin loop from the center of n1 to the center of n2)
     *                   1 means a single spacer (e.g. connecting n1 to n2 in a circle)
     *                   2 means a spacer before n1 and after n2 (e.g. a spacer connecting n1 to a previous nuc and another connecting n2 to a subsequent nuc)
     */
    public float calcOptimalLoopDistance(int nucCount, int endSpacers) {
        return nucCount* 2 * nucRadius + (nucCount + endSpacers - 1) * loopSpacing;
    }

//    public NucSettings get(Nuc n) { return get(n, false); }
//    public NucSettings get(Nuc n, boolean createIfMissing) {
//        if (n == null) return defaultSettings;
//        NucSettings s = nucs.get(n);
//        if (s == null) {
//            if (createIfMissing)
//                nucs.put(n, s = defaultSettings.clone());
//            else
//                s = defaultSettings;
//        }
//        return s;
//    }

    //@Deprecated -- should deprecate this at some point and remove it.
    public float nucleotideRadius(final Nuc nuc) {
        return nucRadius;
    }

    public DrawSettings copy() {
        try {
            DrawSettings ds = (DrawSettings)super.clone();
            ds.onSettingsChanged();
            return ds;
        } catch (CloneNotSupportedException ex) {
            throw new InternalError(ex);
        }
    }
}
