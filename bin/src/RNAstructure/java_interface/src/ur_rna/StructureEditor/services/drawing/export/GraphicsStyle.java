package ur_rna.StructureEditor.services.drawing.export;

import ur_rna.Utilities.annotation.NotNull;
import ur_rna.Utilities.annotation.Nullable;
import ur_rna.Utilities.geom.ShapeUtil;

import java.awt.*;
import java.awt.geom.AffineTransform;
/**
 * Tracks the state of various Graphics properties.
 */
public class GraphicsStyle implements Cloneable {
    public static final Color DEFAULT_BACKGROUND = Color.BLACK;
    public static final Color DEFAULT_COLOR = Color.WHITE;
    public static final Shape DEFAULT_CLIP = null;
    public static final Composite DEFAULT_COMPOSITE = AlphaComposite.SrcOver;
    public static final Font DEFAULT_FONT = Font.decode(null);
    public static final Color DEFAULT_PAINT = DEFAULT_COLOR;
    public static final Stroke DEFAULT_STROKE = new BasicStroke();
    public static final AffineTransform DEFAULT_TRANSFORM = new AffineTransform();
    public static final Color COLOR_EMPTY = new Color(0, true);
    public static final Color DEFAULT_XOR_MODE = COLOR_EMPTY;

    /** Current Rendering hints. */
    private RenderingHints hints;
    /** Current stroke color. */
    private Color color;
    /** Current background color. This only affects the clear operation; It is NOT the fill color for shapes. */
    private Color background;
    /** Method used for mixing colors. */
    private Composite composite;
    /** Stroke used for drawing shapes. */
    private Stroke stroke;
    /** Current font. */
    private Font font;
    /** Paint used to fill shapes. */
    private Paint paint;
    /** Shape used for clipping. */
    private Shape clip;
    /** XOR mode used for rendering. */
    private Color xorMode;
    /** Current transformation matrix. */
    private AffineTransform transform;

    public GraphicsStyle() {
        reset();
    }

    @Override
    public GraphicsStyle clone() {
        try {
            GraphicsStyle clone = (GraphicsStyle) super.clone();
            clone.hints = (RenderingHints) hints.clone();
            clone.clip = ShapeUtil.clone(clip);
            clone.transform = new AffineTransform(transform);
            return clone;
        } catch (CloneNotSupportedException ex) {
            throw  new InternalError(ex);
        }
    }

    public Shape transformShape(Shape shape) {
        return ShapeUtil.transform(shape, transform);
    }

    public Shape untransformShape(Shape shape) {
        return ShapeUtil.inverseTransform(shape, transform);
    }

    public RenderingHints getHints() {
        return hints;
    }

    /** Background color only affects the clear operation; It is NOT the fill color of shapes. */
    public Color getBackground() {
        return background;
    }

    /** Background color only affects the clear operation; It is NOT the fill color of shapes. */
    public void setBackground(Color background) {
        this.background = background;
    }

    public Color getColor() {
        return color;
    }

    public void setColor(Color color) {
        this.color = color;
        this.paint = color;
    }

    public Shape getClip() { return clip; }

    public void setClip(Shape clip) {
        this.clip = clip;
    }

    public Font getFont() {
        return font;
    }

    public void setFont(Font font) {
        this.font = font;
    }

    public Stroke getStroke() {
        return stroke;
    }

    public void setStroke(Stroke stroke) {
        this.stroke = stroke;
    }

    public Paint getPaint() {
        return paint;
    }

    public void setPaint(Paint paint) {
        this.paint = paint;
    }

    /**
     * Returns the current XOR background color if XOR Mode is in effect.
     * Otherwise returns null.
     */
    @NotNull
    public Color getXorMode() {
        return xorMode;
    }

    public void setXorMode(@Nullable Color xorMode) { this.xorMode = xorMode == null ? DEFAULT_XOR_MODE : xorMode; }

    public Composite getComposite() {
        return composite;
    }

    public void setComposite(Composite composite) {
        this.composite = composite;
        this.xorMode = COLOR_EMPTY; // indicate that XOR mode is no longer in effect.
    }


    public void setTransform(AffineTransform tx) { transform.setTransform(tx); }

    public AffineTransform getTransform() { return transform;  }

    @Override
    public boolean equals(Object obj) {
        return obj instanceof GraphicsStyle &&
                equals((GraphicsStyle)obj);
    }

    public boolean equals(GraphicsStyle gc) {
        return gc == this || gc != null &&
                hints.equals(gc.hints) &&
                background.equals(gc.background) &&
                color.equals(gc.color) &&
                composite.equals(gc.composite) &&
                font.equals(gc.font) &&
                paint.equals(gc.paint) &&
                stroke.equals(gc.stroke) &&
                transform.equals(gc.transform) &&
                xorMode.equals(gc.xorMode) &&
                (clip == gc.clip || clip != null && clip.equals(gc.clip));
    }

    public boolean isDefault() {
        return hints.isEmpty() &&
                DEFAULT_XOR_MODE.equals(xorMode) &&
                DEFAULT_BACKGROUND.equals(background) &&
                DEFAULT_COLOR.equals(color) &&
                DEFAULT_COMPOSITE.equals(composite) &&
                DEFAULT_FONT.equals(font) &&
                DEFAULT_PAINT.equals(paint) &&
                DEFAULT_STROKE.equals(stroke) &&
                clip == DEFAULT_CLIP &&
                DEFAULT_TRANSFORM.equals(transform);
    }
    /** Reset all properties to their default values. */
    public void reset() {
        hints = new RenderingHints(null);
        background = DEFAULT_BACKGROUND;
        color = DEFAULT_COLOR;
        clip = DEFAULT_CLIP;
        composite = DEFAULT_COMPOSITE;
        font = DEFAULT_FONT;
        paint = DEFAULT_PAINT;
        stroke = DEFAULT_STROKE;
        xorMode = DEFAULT_XOR_MODE;
        transform = new AffineTransform(DEFAULT_TRANSFORM);
    }
    public void copyFrom(final Graphics2D graphics) {
        transform = graphics.getTransform();
        clip = transformShape(graphics.getClip());
        hints = graphics.getRenderingHints();
        background = graphics.getBackground();
        color = graphics.getColor();
        paint = graphics.getPaint();
        stroke = graphics.getStroke();
        font = graphics.getFont();
//        Composite c = graphics.getComposite();
        composite = graphics.getComposite();
        if (composite == AlphaComposite.Xor) // In practice, I don't think AlphaComposite.Xor is used. The Composite used for setXORMode is sun.java2d.loops.XORComposite (at least on Windows, Oracle Java 1.8)
            xorMode = Color.GRAY;

        //TODO: XORMode is not supported. If graphics.setXORMode() was called prior to this, we would not be able to tell.
        // The reason is that the graphics interface has no method to getXORMode. Instead calling setXORMode affects the Composite.
        // But the XORComposite class (sun.java2d.loops.XORComposite) is proprietary and possibly OS-dependent.
        // For example, it does not exist in OpenJDK.
        // We could use reflection to get the class name and also to get XORComposite#getXorColor();
//        if (c instanceof XORComposite) {
//            xorMode = ((XORComposite) c).getXorColor();
//        } else if (c == AlphaComposite.Xor) {
//            xorMode = Color.GRAY;
//        }
    }
}
