package ur_rna.StructureEditor.models;

import java.awt.*;

/**
 * Contains formatting information, which affects how nucleotides are drawn.
 */
public class NucStyle implements Cloneable {
    //transient int index = -1; // an arbitrary number used to link each nucleotide with a Style when serialized/deserialized.
    //public int getIndex() { return index; }

    public Color fillColor;
    public Color textColor;
    public Color outlineColor;
    public Color bondColor;
    public Font font;
    public float numberOffset;

    public NucStyle() { }
    public NucStyle(final Color fillColor, final Color textColor, final Color outlineColor, final Color bondColor, final Font font, final float numberOffset) {
        this.fillColor = fillColor;
        this.textColor = textColor;
        this.outlineColor = outlineColor;
        this.bondColor = bondColor;
        this.font = font;
        this.numberOffset = numberOffset;
    }
    // public float outlineThickness;
    @Override
    protected NucStyle clone() {
        //try {
        NucStyle s = new NucStyle();
        copyTo(s);
        return s;
        //return new NucStyle(fillColor, textColor, outlineColor, font);
            //return (NucStyle)super.clone();
//        } catch (CloneNotSupportedException ex) {
//            throw new InternalError(ex);
//        }
    }
    public boolean isEmpty() {
        return font==null &&
                fillColor==null &&
                textColor==null &&
                outlineColor==null &&
                bondColor==null &&
                numberOffset==0;
    }
    public void copyTo(final NucStyle s) {
        s.fillColor = fillColor;
        s.textColor = textColor;
        s.outlineColor = outlineColor;
        s.bondColor = bondColor;
        s.font = font;
        s.numberOffset = numberOffset;
    }
}
