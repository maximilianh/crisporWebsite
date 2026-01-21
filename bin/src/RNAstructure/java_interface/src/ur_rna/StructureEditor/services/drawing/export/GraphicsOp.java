package ur_rna.StructureEditor.services.drawing.export;

import ur_rna.Utilities.geom.ShapeUtil;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.util.Collections;

/**
 * Represents a graphics operation that can be applied to a Graphics2D object.
 */
public interface GraphicsOp {
    void replay(Graphics2D g);
    void apply(GraphicsStyle c);

    /** Represents an operation that does not change the state of the graphics context. */
    abstract class ActionOp implements GraphicsOp {
        @Override
        public abstract void replay(Graphics2D g);
        @Override
        public void apply(GraphicsStyle c) { /* This operation does not change the state of the graphics context. */ }
    }

    /**
     * Represents an operation that changes the state of the graphics context,
     * but does not perform any drawing on its own.
     */
    abstract class StateOp implements GraphicsOp {    }

    /**
     * Represents an operation that is for control purposes only.
     * It does not change the state of the graphics context or
     * perform drawing.
     */
    abstract class ControlOp implements GraphicsOp {
        public void replay(Graphics2D g) {}
        public void apply(GraphicsStyle c) { }
    }

    //<editor-fold desc="Control Operations">
    /**
     * Signifies the start of a new graphics record.
     * Clients should reset their Graphics context to the default state.
     */
    class InitOp extends ControlOp { }
    /**
     * Indicates the end of a graphics record.
     * Clients should revert to the state at which InitOp was encountered.
     */
    class EndOp extends ControlOp {}

    /**
     * Indicates the start of a new group of graphics objects (e.g. for SVG)
     */
    class GroupStartOp extends ControlOp {}
    /**
     * Indicates the end of group of graphics objects (e.g. for SVG)
     */
    class GroupEndOp extends ControlOp {}
    //</editor-fold>

    //<editor-fold desc="Action Operations">
    class DrawShapeOp extends ActionOp {
        private Shape shape;
        private boolean fill;

        public Shape getShape() {
            return shape;
        }
        public boolean isFill() {
            return fill;
        }

        public DrawShapeOp(final Shape shape, final boolean fill) {
            this.shape = shape;
            this.fill = fill;
        }
        @Override
        public void replay(final Graphics2D g) {
            if (fill)
                g.fill(shape);
            else
                g.draw(shape);
        }
    }
    class DrawStringOp extends ActionOp {
        public String text;
        public float x, y;

        public String getText() {
            return text;
        }
        public float getX() {
            return x;
        }
        public float getY() {
            return y;
        }
        public DrawStringOp(final String text, final float x, final float y) {
            this.text = text;
            this.x = x;
            this.y = y;
        }

        @Override
        public void replay(final Graphics2D g) {
            g.drawString(text, x, y);
        }
    }
    class DrawImageOp extends ActionOp {
        private Image img;
        private int width;
        private int height;
        private int originalWidth;
        private int originalHeight;
        private int x;
        private int y;

        public Image getImage() {
            return img;
        }
        public int getWidth() {
            return width;
        }
        public int getHeight() {
            return height;
        }
        public int getOriginalWidth() {
            return originalWidth;
        }
        public int getOriginalHeight() {
            return originalHeight;
        }
        public int getX() {
            return x;
        }
        public int getY() {
            return y;
        }
        public DrawImageOp(final Image img, final int originalWidth, final int originalHeight, final int x, final int y, final int width, final int height) {
            this.img = img;
            this.originalWidth = originalWidth;
            this.originalHeight = originalHeight;
            this.x = x;
            this.y = y;
            this.width = width;
            this.height = height;
        }

        @Override
        public void replay(final Graphics2D g) {
            g.drawImage(img, x, y, width, height, null);
        }
    }
    class ClearOp extends ActionOp {
        private Rectangle rc;
        public ClearOp(final Rectangle rc) {
            this.rc = rc;
        }
        public Rectangle getRect() {
            return rc;
        }
        @Override
        public void replay(final Graphics2D g) {
            g.clearRect(rc.x, rc.y, rc.width, rc.height);
        }
    }
    //</editor-fold>

    //<editor-fold desc="State Operations">
    class BackgroundOp extends StateOp {
        private Color color;
        public BackgroundOp(final Color color) {
            this.color = color;
        }
        public Color getColor() {
            return color;
        }
        @Override
        public void replay(final Graphics2D g) {
            g.setBackground(color);
        }
        @Override
        public void apply(final GraphicsStyle c) {
            c.setBackground(color);
        }
    }
    class FontOp extends StateOp {
        private Font font;
        public FontOp(final Font font) {
            this.font = font;
        }
        public Font getFont() {
            return font;
        }
        @Override
        public void replay(final Graphics2D g) {
            g.setFont(font);
        }
        @Override
        public void apply(final GraphicsStyle c) {
            c.setFont(font);
        }
    }
    class ColorOp extends StateOp {
        private Color color;
        public ColorOp(final Color color) {
            this.color = color;
        }
        public Color getColor() {
            return color;
        }
        @Override
        public void replay(final Graphics2D g) { g.setColor(color); }
        @Override
        public void apply(final GraphicsStyle c) { c.setColor(color); }
    }
    class XorModeOp extends StateOp {
        private Color color;
        public XorModeOp(final Color color) {
            this.color = color;
        }
        public Color getColor() {
            return color;
        }
        @Override
        public void replay(final Graphics2D g) {
            g.setXORMode(color);
        }
        @Override
        public void apply(final GraphicsStyle c) {
            c.setXorMode(color);
        }
    }
    class HintOp extends StateOp {
        private RenderingHints.Key key;
        private Object value;
        public HintOp(final RenderingHints.Key key, final Object value) {
            this.key = key;
            this.value = value;
        }
        public RenderingHints.Key getKey() {
            return key;
        }
        public Object getValue() {
            return value;
        }
        @Override
        public void replay(final Graphics2D g) {
            g.setRenderingHint(key, value);
        }
        @Override
        public void apply(final GraphicsStyle c) { c.getHints().put(key, value); }
    }
    class ClearHintsOp extends StateOp {
        @Override
        public void replay(final Graphics2D g) {
            g.setRenderingHints(Collections.emptyMap());
        }
        @Override
        public void apply(final GraphicsStyle c) {
            c.getHints().clear();
        }
    }
    class PaintOp extends StateOp {
        private Paint paint;
        public PaintOp(final Paint paint) {
            this.paint = paint;
        }
        public Paint getPaint() {
            return paint;
        }
        @Override
        public void replay(final Graphics2D g) {
            g.setPaint(paint);
        }
        @Override
        public void apply(final GraphicsStyle c) { c.setPaint(paint); }
    }
    class TransformOp extends StateOp {
        private AffineTransform tr;
        private boolean append;
        public TransformOp(final AffineTransform tr, boolean append) {
            this.tr = tr;
            this.append = append;
        }
        public AffineTransform getTransform() {
            return tr;
        }
        public boolean isAppend() {
            return append;
        }
        @Override
        public void replay(final Graphics2D g) {
            if (append)
                g.transform(tr);
            else
                g.setTransform(tr);
        }
        @Override
        public void apply(final GraphicsStyle c) {
            if (append) {
                c.getTransform().concatenate(tr);
            } else
                c.setTransform(tr);
        }
    }
    class CompositeOp extends StateOp {
        private Composite composite;
        public CompositeOp(final Composite composite) {
            this.composite = composite;
        }
        public Composite getComposite() {
            return composite;
        }
        @Override
        public void replay(final Graphics2D g) {
            g.setComposite(composite);
        }
        @Override
        public void apply(final GraphicsStyle c) { c.setComposite(composite); }
    }
    class ClipOp extends StateOp {
        private Shape clip;
        private boolean append;
        public ClipOp(final Shape clip, boolean append) {
            this.clip = clip;
            this.append = append;
        }
        public Shape getClip() {
            return clip;
        }
        public boolean isAppend() {
            return append;
        }
        @Override
        public void replay(final Graphics2D g) {
            if (append)
                g.clip(clip);
            else
                g.setClip(clip);
        }

        @Override
        public void apply(final GraphicsStyle c) {
            Shape newClip = clip == null ? null : c.transformShape(clip); // transform the clip, based on the current transformation.
            Shape oldClip = c.getClip();
            if (append && clip != null && oldClip != null)
                newClip = ShapeUtil.intersection(oldClip, newClip);
            c.setClip(newClip);
        }
    }
    class StrokeOp extends StateOp {
        Stroke stroke;
        public StrokeOp(final Stroke stroke) {
            this.stroke = stroke;
        }
        public Stroke getStroke() {
            return stroke;
        }
        @Override
        public void replay(final Graphics2D g) { g.setStroke(stroke); }
        public void apply(final GraphicsStyle c) { c.setStroke(stroke);  }
    }
    //</editor-fold>
}
