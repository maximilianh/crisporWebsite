package ur_rna.StructureEditor.services.drawing.export;

import ur_rna.StructureEditor.services.drawing.export.GraphicsOp.*;
import ur_rna.Utilities.annotation.NotNull;
import ur_rna.Utilities.geom.ShapeUtil;
import ur_rna.Utilities.swing.ImageUtil;

import java.awt.*;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.*;
import java.awt.image.*;
import java.awt.image.renderable.RenderableImage;
import java.text.AttributedCharacterIterator;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * A class that can be substituted for a Graphics2D object,
 * but instead of performing any actual drawing, it simply records
 * the drawing instructions. These instructions can later be replayed
 * onto a "real" Graphics2D object or converted into vector graphics
 * by an appropriate converter.
 */
public class Graphics2DRecorder extends Graphics2D implements Cloneable {
    private final List<GraphicsOp> ops;
    /** Device configuration settings. */
    private final GraphicsConfiguration deviceConfig;
    /** Context settings used to render fonts. */
    private final FontRenderContext fontRenderContext;
    /** Flag that tells whether this graphics object has been disposed. */
    private boolean disposed;

    private GraphicsStyle ctx;
    private int suspendRecordingCounter; // recording only occurs while this is 0.

    public Graphics2DRecorder() {
        this.ops = new ArrayList<>();
        GraphicsEnvironment graphicsEnvironment = GraphicsEnvironment.getLocalGraphicsEnvironment();
        GraphicsDevice graphicsDevice;
        if (graphicsEnvironment.isHeadlessInstance())
            deviceConfig = null;
        else {
            graphicsDevice = graphicsEnvironment.getDefaultScreenDevice();
            deviceConfig = graphicsDevice.getDefaultConfiguration();
        }
        fontRenderContext = new FontRenderContext(null, true, true);
        ctx = new GraphicsStyle();
//        setColor(Color.BLACK); // Required for EPS, PDF, and SVG
//        setStroke(new BasicStroke(1f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10f, null, 0f)); // EPS and PDF
        add(new InitOp());
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        try {
            Graphics2DRecorder clone = (Graphics2DRecorder) super.clone();
            clone.ctx = (GraphicsStyle) ctx.clone();
            return clone;
        } catch (CloneNotSupportedException e) {
            return null;
        }
    }

    @Override
    public void addRenderingHints(Map<?, ?> hints) {
        if (isDisposed())
            return;
        for (Map.Entry<?, ?> entry : hints.entrySet())
            setRenderingHint((RenderingHints.Key) entry.getKey(), entry.getValue());
    }

    @Override
    public void clip(Shape s) {
        add(new ClipOp(s, true));
    }

    @Override
    public void draw(Shape s) {
        draw(s, true);
    }

    /**
     * The the outline of the specified Shape.
     * @param s The shape to draw.
     * @param createCopy If true, a copy of the shape will be made so that if it changes
     *                   after the call, the recorded operation will still use the copy of the shape
     *                   at the time of the call to draw().
     *                   This should be set to true, unless it is known that the shape will never
     *                   be changed (i.e. if it is created internally) or if any changes are intentionally
     *                   meant to affect later application of the recorded operation.
     */
    public void draw(Shape s, boolean createCopy) {
        if (isDisposed() || s == null)
            return;
        if (createCopy) s = ShapeUtil.clone(s);
        add(new DrawShapeOp(s, false));
    }

    @Override
    public void drawGlyphVector(GlyphVector g, float x, float y) {
        draw(g.getOutline(x, y), false);
    }

    @Override
    public boolean drawImage(Image img, AffineTransform tr, ImageObserver obs) {
        BufferedImage b = getTransformedImage(img, tr);
        return drawImage(b, b.getMinX(), b.getMinY(), b.getWidth(), b.getHeight(), null, null);
    }

    /**
     * Returns a transformed version of an image.
     * @param image Image to be transformed
     * @param tr Affine transform to be applied
     * @return Image with transformed content
     */
    private BufferedImage getTransformedImage(Image image, AffineTransform tr) {
        AffineTransformOp op = new AffineTransformOp(tr, ctx.getHints());
        BufferedImage bufferedImage = ImageUtil.toBufferedImage(image);
        return op.filter(bufferedImage, null);
    }

    @Override
    public void drawImage(BufferedImage img, BufferedImageOp op, int x, int y) {
        if (op != null) {
            img = op.filter(img, null);
        }
        drawImage(img, x, y, img.getWidth(), img.getHeight(), null, null);
    }

    @Override
    public void drawRenderableImage(RenderableImage img, AffineTransform tr) {
        drawRenderedImage(img.createDefaultRendering(), tr);
    }

    @Override
    public void drawRenderedImage(RenderedImage img, AffineTransform tr) {
        drawImage(ImageUtil.toBufferedImage(img), tr, null);
    }

    @Override
    public void drawString(String str, int x, int y) {
        drawString(str, (float) x, (float) y);
    }

    @Override
    public void drawString(String str, float x, float y) {
        if (isDisposed() || str == null || str.trim().length() == 0) {
            return;
        }
//        if (isTextAsVectors) {
//            TextLayout layout = new TextLayout(str, getFont(), getFontRenderContext());
//            Shape s = layout.getOutline(AffineTransform.getTranslateInstance(x, y));
//            fill(s);
//        } else {
            add(new DrawStringOp(str, x, y));
//        }
    }

    @Override
    public void drawString(AttributedCharacterIterator iterator, int x, int y) {
        drawString(iterator, (float) x, (float) y);
    }

    @Override
    public void drawString(AttributedCharacterIterator iterator, float x,
            float y) {
        // TODO Draw styled text
        StringBuilder buf = new StringBuilder();
        for (char c = iterator.first(); c != AttributedCharacterIterator.DONE;
             c = iterator.next()) {
            buf.append(c);
        }
        drawString(buf.toString(), x, y);
    }

    @Override
    public void fill(Shape s) {
        fill(s, true);
    }
    /**
     * Fill the specified Shape.
     * @param s The shape to fill.
     * @param createCopy If true, a copy of the shape will be made so that if it changes
     *                   after the call, the recorded operation will still use the copy of the shape
     *                   at the time of the call to fill().
     *                   This should be set to true, unless it is known that the shape will never
     *                   be changed (i.e. if it is created internally) or if any changes are intentionally
     *                   meant to affect later application of the recorded operation.
     */
    public void fill(Shape s, boolean createCopy)  {
        if (isDisposed() || s == null) return;
        if (createCopy)
            s = ShapeUtil.clone(s);
        add(new DrawShapeOp(s, true));
    }

    @Override
    public Color getBackground() {
        return ctx.getBackground();
    }

    @Override
    public Composite getComposite() {
        return ctx.getComposite();
    }

    @Override
    public GraphicsConfiguration getDeviceConfiguration() {
        return deviceConfig;
    }

    @Override
    public FontRenderContext getFontRenderContext() {
        return fontRenderContext;
    }

    @Override
    public Paint getPaint() {
        return ctx.getPaint();
    }

    @Override
    public Object getRenderingHint(RenderingHints.Key hintKey) {
        return ctx.getHints().get(hintKey);
    }

    @Override
    public RenderingHints getRenderingHints() {
        return ctx.getHints();
    }

    @Override
    public Stroke getStroke() {
        return ctx.getStroke();
    }

    @Override
    public boolean hit(Rectangle rect, Shape s, boolean onStroke) {
        Shape hitShape = s;
        if (onStroke) {
            hitShape = getStroke().createStrokedShape(hitShape);
        }
        hitShape = ctx.transformShape(hitShape);
        return hitShape.intersects(rect);
    }

    @Override
    public void setBackground(Color color) {
        if (isDisposed() || color == null || getColor().equals(color)) {
            return;
        }
        add(new BackgroundOp(color));
    }

    @Override
    public void setComposite(Composite comp) {
        if (isDisposed()) return;
        if (comp == null)
            throw new IllegalArgumentException("Cannot set a null composite.");
        add(new CompositeOp(comp));
    }

    @Override
    public void setPaint(Paint paint) {
        if (isDisposed() || paint == null)
            return;
        if (paint instanceof Color) {
            setColor((Color) paint);
            return;
        }
        if (getPaint().equals(paint))
            return;
        add(new PaintOp(paint));
    }

    @Override
    public void setRenderingHint(RenderingHints.Key hintKey, Object hintValue) {
        if (isDisposed()) return;
        add(new HintOp(hintKey, hintValue));
    }

    @Override
    public void setRenderingHints(Map<?, ?> hints) {
        if (isDisposed()) return;
        add(new ClearHintsOp());
        for (Map.Entry<?, ?> hint : hints.entrySet())
            setRenderingHint((RenderingHints.Key) hint.getKey(), hint.getValue());
    }

    @Override
    public void setStroke(Stroke s) {
        if (s == null)
            throw new IllegalArgumentException("Cannot set a null stroke.");
        if (isDisposed()) return;
        add(new StrokeOp(s));
    }

    @Override
    public AffineTransform getTransform() {
        return new AffineTransform(ctx.getTransform());
    }

    @Override
    public void setTransform(AffineTransform tx) {
        if (isDisposed() || tx == null || ctx.getTransform().equals(tx))
            return;
        add(new TransformOp(tx, false));
    }

    @Override
    public void shear(double shx, double shy) {
        if (shx != 1.0 || shy != 1.0)
            transform(AffineTransform.getShearInstance(shx, shy));
    }

    @Override
    public void transform(AffineTransform Tx) {
        if (!Tx.isIdentity())
            add(new TransformOp(Tx, true));
    }

    @Override
    public void translate(int x, int y) {
        translate((double) x, (double) y);
    }

    @Override
    public void translate(double tx, double ty) {
        if (tx != 0.0 || ty != 0.0)
            transform(AffineTransform.getTranslateInstance(tx, ty));
    }

    @Override
    public void rotate(double theta) {
        rotate(theta, 0.0, 0.0);
    }

    @Override
    public void rotate(double theta, double x, double y) {
        if (theta != 0.0)
            transform(AffineTransform.getRotateInstance(theta, x, y));
    }

    @Override
    public void scale(double sx, double sy) {
        if (sx != 1.0 || sy != 1.0)
            transform(AffineTransform.getScaleInstance(sx, sy));
    }

    @Override
    public void clearRect(int x, int y, int width, int height) {
        add(new ClearOp(new Rectangle(x, y, width, height)));
//        Color colorOld = getColor();
//        setColor(getBackground());
//        fillRect(x, y, width, height);
//        setColor(colorOld);
    }

    @Override
    public void clipRect(int x, int y, int width, int height) {
        clip(new Rectangle(x, y, width, height));
    }

    @Override
    public void copyArea(int x, int y, int width, int height, int dx, int dy) {
        throw notImpl();
    }

    @Override
    public Graphics create() {
        if (isDisposed()) return null;
        try {
            return (Graphics2DRecorder) this.clone();
            //add(new CreateCommand(clone));
        } catch (CloneNotSupportedException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public void dispose() {
        if (isDisposed()) return;
        add(new EndOp());
        disposed = true;
    }

    @Override
    public void drawArc(int x, int y, int width, int height, int startAngle, int arcAngle) {
        draw(new Arc2D.Double(x, y, width, height, startAngle, arcAngle, Arc2D.OPEN), false);
    }

    @Override
    public boolean drawImage(Image img, int x, int y, ImageObserver observer) {
        return drawImage(img, x, y, img.getWidth(observer), img.getHeight(observer), null, observer);
    }

    @Override
    public boolean drawImage(Image img, int x, int y, Color bgcolor, ImageObserver observer) {
        return drawImage(img, x, y, img.getWidth(observer), img.getHeight(observer), bgcolor, observer);
    }

    @Override
    public boolean drawImage(Image img, int x, int y, int width, int height, ImageObserver observer) {
        return drawImage(img, x, y, width, height, null, observer);
    }

    @Override
    public boolean drawImage(Image img, int x, int y, int width, int height, Color bgcolor, ImageObserver observer) {
        if (isDisposed() || img == null)
            return true;

        int imageWidth = img.getWidth(observer);
        int imageHeight = img.getHeight(observer);
        Rectangle bounds = new Rectangle(x, y, width, height);

        if (bgcolor != null) {
            // Fill rectangle with bgcolor
            Color bgcolorOld = getColor();
            setColor(bgcolor);
            fill(bounds, false);
            setColor(bgcolorOld);
        }
        add(new DrawImageOp(img, imageWidth, imageHeight, x, y, width, height));
        return true;
    }

    @Override
    public boolean drawImage(Image img, int dx1, int dy1, int dx2, int dy2,
            int sx1, int sy1, int sx2, int sy2, ImageObserver observer) {
        return drawImage(img, dx1, dy1, dx2, dy2, sx1, sy1, sx2, sy2, null, observer);
    }

    @Override
    public boolean drawImage(Image img, int dx1, int dy1, int dx2, int dy2,
            int sx1, int sy1, int sx2, int sy2, Color bgcolor, ImageObserver observer) {
        if (img == null)
            return true;
        int sx = Math.min(sx1, sx2);
        int sy = Math.min(sy1, sy2);
        int dx = Math.min(dx1, dx2);
        int dy = Math.min(dy1, dy2);
        int dw = Math.abs(dx2 - dx1);
        int dh = Math.abs(dy2 - dy1);

        // Draw image on rectangle
        Image cropped = ImageUtil.toBufferedImage(img).getSubimage(sx, sy, Math.abs(sx2 - sx1), Math.abs(sy2 - sy1));
        return drawImage(cropped, dx, dy, dw, dh, bgcolor, observer);
    }

    @Override
    public void drawLine(int x1, int y1, int x2, int y2) {
        draw(new Line2D.Double(x1, y1, x2, y2), false);
    }

    @Override
    public void drawOval(int x, int y, int width, int height) {
        draw(new Ellipse2D.Double(x, y, width, height), false);
    }

    @Override
    public void drawPolygon(Polygon p) {
        draw(p, true);
    }

    @Override
    public void drawPolygon(int[] xPoints, int[] yPoints, int nPoints) {
        draw(new Polygon(xPoints, yPoints, nPoints), false);
    }

    @Override
    public void drawPolyline(int[] xPoints, int[] yPoints, int nPoints) {
        Path2D p = new Path2D.Float();
        for (int i = 0; i < nPoints; i++) {
            if (i > 0) {
                p.lineTo(xPoints[i], yPoints[i]);
            } else {
                p.moveTo(xPoints[i], yPoints[i]);
            }
        }
        draw(p, false);
    }

    @Override
    public void drawRect(int x, int y, int width, int height) {
        draw(new Rectangle(x, y, width, height), false);
    }

    @Override
    public void drawRoundRect(int x, int y, int width, int height,
            int arcWidth, int arcHeight) {
        draw(new RoundRectangle2D.Double(x, y, width, height,
                arcWidth, arcHeight), false);
    }

    @Override
    public void fillArc(int x, int y, int width, int height,
            int startAngle, int arcAngle) {
        fill(new Arc2D.Double(x, y, width, height,
                startAngle, arcAngle, Arc2D.PIE), false);
    }

    @Override
    public void fillOval(int x, int y, int width, int height) {
        fill(new Ellipse2D.Double(x, y, width, height), false);
    }

    @Override
    public void fillPolygon(Polygon p) {
        fill(p, true);
    }

    @Override
    public void fillPolygon(int[] xPoints, int[] yPoints, int nPoints) {
        fill(new Polygon(xPoints, yPoints, nPoints), false);
    }

    @Override
    public void fillRect(int x, int y, int width, int height) {
        fill(new Rectangle(x, y, width, height), false);
    }

    @Override
    public void fillRoundRect(int x, int y, int width, int height,
            int arcWidth, int arcHeight) {
        fill(new RoundRectangle2D.Double(x, y, width, height,
                arcWidth, arcHeight), false);
    }

    @Override
    public Shape getClip() {
        return ctx.untransformShape(ctx.getClip());
    }

    @Override
    public Rectangle getClipBounds() {
        if (getClip() == null) {
            return null;
        }
        return getClip().getBounds();
    }

    @Override
    public Color getColor() {
        return ctx.getColor();
    }

    @Override
    public Font getFont() {
        return ctx.getFont();
    }

    @Override
    public FontMetrics getFontMetrics(Font f) {
        BufferedImage bi = new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB_PRE);
        Graphics g = bi.getGraphics();
        FontMetrics fontMetrics = g.getFontMetrics(getFont());
        g.dispose();
        return fontMetrics;
    }

    @Override
    public void setClip(Shape clip) {
        if (isDisposed()) return;
        add(new ClipOp(clip, false));
    }

    @Override
    public void setClip(int x, int y, int width, int height) {
        setClip(new Rectangle(x, y, width, height));
    }

    @Override
    public void setColor(Color c) {
        if (isDisposed() || c == null || getColor().equals(c)) {
            return;
        }
        add(new ColorOp(c));
    }

    @Override
    public void setFont(Font font) {
        if (isDisposed() || (font != null && getFont().equals(font))) {
            return;
        }
        add(new FontOp(font));
        ctx.setFont(font);
    }

    @Override
    public void setPaintMode() {
        setComposite(AlphaComposite.SrcOver);
    }

    /**
     * Returns the current XOR background color.
     * If the value is {@link GraphicsStyle#COLOR_EMPTY } then
     * XOR Mode is not in effect.
     */
    @NotNull
    public Color getXORMode() {
        return ctx.getXorMode();
    }

    @Override
    public void setXORMode(Color backColor) {
        if (isDisposed() || backColor == null)
            return;
        add(new XorModeOp(backColor));
    }

    /** Add the operation and apply it to the current graphics context. */
    private void add(GraphicsOp op) {
        op.apply(ctx);
        if (suspendRecordingCounter==0) ops.add(op);
    }
//    private void add(GraphicsOp op, boolean preventApply) {
//        op.apply(ctx);
//        ops.add(op);
//    }

    protected boolean isDisposed() {
        return disposed;
    }

    private RuntimeException notImpl() {
        return new UnsupportedOperationException("Graphics operation not yet implemented.");
    }

    /**
     * Returns a {@code CommandSequence} representing all calls that were issued to this {@code VectorGraphics2D} object.
     * @return Sequence of commands since.
     */
    public List<GraphicsOp> getGraphicsOps() {
        return ops;
    }
    public void replay(final Graphics2D g) {
        for (GraphicsOp o : ops)
            o.replay(g);
    }
    public void suspendRecording() {
        suspendRecordingCounter++;
    }
    public void resumeRecording() {
        if (suspendRecordingCounter == 0)
            throw new IllegalStateException("Cannot resume. Already recording. This probably means you have unbalanced calls to suspendRecording and resumeRecording.");
        suspendRecordingCounter--;
    }
    public GraphicsStyle getContext() { return ctx; }
    public void startGroup() {
        add(new GroupStartOp());
    }
    public void endGroup() {
        add(new GroupEndOp());
    }
}
