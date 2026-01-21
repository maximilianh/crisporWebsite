package ur_rna.StructureEditor.services.drawing;

import ur_rna.StructureEditor.models.DrawSettings;
import ur_rna.Utilities.geom.GraphicsUtil;

import java.awt.*;
import java.awt.geom.Point2D;

import static ur_rna.Utilities.Colors.setAlpha;

/**
 * Provides a way to specify the drawing of shapes and text outside of the actual rendering process.
 * This allows events like mouse-drags to use the event information to create visual cues (drawing hints)
 * to the user.
 */
public abstract class DrawHint {
    public static BasicStroke HintLineDefault = new BasicStroke(1);
    public static BasicStroke HintLineThick = new BasicStroke(3, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
    public static final Color HintColorDefault = new Color(0, true);
    public static final Color HintColorReference = new Color(0, true);

//    public static Color hintFill = RnaDrawController.Colors.setAlpha(RnaDrawController.Colors.cool, 0x99);

    public Color colorFill = null;
    public Color colorLine = HintColorDefault;
    public boolean isModelCoords = false;
    public Stroke lineStyle = HintLineDefault;

    public void draw(Graphics2D g, View2D view, final DrawSettings settings) {
//        if (!isDeviceCoords) {
//            AffineTransform prevTr = g.getTransform();
//            g.transform(view.trToScreen);
//            drawHint(g, view, prevTr, settings);
//            g.setTransform(prevTr);
//        } else
//            drawHint(g, view, null, settings);
        drawHint(g, view, settings);
    }

    protected abstract void drawHint(Graphics2D g, final View2D view, final DrawSettings settings);
    protected Color getDrawColor(Color c, DrawSettings s) {
        if (colorLine == HintColorDefault)
            return s.hint;
        else if (colorLine == HintColorReference)
            return s.hintEmphasis;
        else
            return c;
    }

    public DrawHint fromModel() { isModelCoords = true; return this; }

    public DrawHint color(Color line) { colorLine = line; return this; }
    public DrawHint color(Color line, Color fill) { colorLine = line; colorFill = fill; return this; }
    public DrawHint color(Color line, int fillAlpha) { colorLine = line; colorFill = setAlpha(line, fillAlpha); return this; }
    public DrawHint fill(Color fill) { colorFill = fill; return this; }
    public DrawHint fillAlpha(int fillAlpha) { colorFill = setAlpha(colorLine, fillAlpha); return this; }
    public DrawHint noFill() { colorFill = null; return this; }
    public DrawHint noLine() { colorLine = null; return this; }
    public DrawHint line(float width, float... dash) {
        BasicStroke copy = lineStyle instanceof BasicStroke ? (BasicStroke)lineStyle : HintLineDefault;
        lineStyle = new BasicStroke(width, copy.getEndCap(), copy.getLineJoin(), copy.getMiterLimit(), dash == null ? copy.getDashArray() : dash, copy.getDashPhase());
        return this;
    }
    public DrawHint line(float width) { return line(width, (float[]) null); }
    public DrawHint line(Stroke s) {
        this.lineStyle = s;
        return this;
    }

    /** Do something to emphasize this hint. */
    public abstract DrawHint bold();
    public DrawHint dashed() { return line(1, 2, 4); }

    public static class DrawShape extends DrawHint {
        public final Shape shape;
        public DrawShape(final Shape shape) { this.shape = shape; }
        @Override
        protected void drawHint(final Graphics2D g, final View2D view, final DrawSettings settings) {
            Shape screenShape = isModelCoords ? view.trToScreen.createTransformedShape(shape) : shape;
            if (colorFill != null) {
                g.setColor(getDrawColor(colorFill, settings));
                g.fill(screenShape);
            }
            if (lineStyle != null && colorLine != null) {
                g.setStroke(lineStyle);
                g.setColor(getDrawColor(colorLine, settings));
                g.draw(screenShape);
            }
        }

        @Override
        public DrawHint bold() {
            lineStyle = HintLineThick;
            return this;
        }
    }
    public static class DrawText extends DrawHint {
        public final String text;
        public final Point2D textPos;
        public final Point2D.Float screenPos = new Point2D.Float();
        private boolean isBold;
        //public int textAlign;
        public DrawText(final String text, final Point2D textPos) {
            this.text = text;
            this.textPos = textPos;
        }
        @Override
        protected void drawHint(final Graphics2D g, final View2D view, final DrawSettings settings) {
            if (colorLine != null && text != null) {
                g.setColor(getDrawColor(colorLine, settings));
                if (isModelCoords)
                    view.toScreen(textPos, screenPos);
                else
                    screenPos.setLocation(textPos);
                Font tmpFont = null;

                if (isBold) {
                    tmpFont = g.getFont();
                    g.setFont(tmpFont.deriveFont(tmpFont.getStyle() | Font.BOLD));
                }
                GraphicsUtil.drawTextCentered(g, text, screenPos);
                if (isBold)
                    g.setFont(tmpFont);
            }
        }
        @Override
        public DrawHint bold() {
            isBold = true;
            return this;
        }
    }

    public static DrawHint shape(final Shape shape) {
        return new DrawShape(shape);
    }
    public static DrawHint text(final String text, final Point2D textPos) {
        return new DrawText(text, textPos);
    }

//        public DrawHint(final Shape shape, final Color colorLine, final Color colorFill, final boolean isDeviceCoords) {
//            this.shape = shape;
//            this.colorLine = colorLine;
//            this.colorFill = colorFill;
//            this.isDeviceCoords = isDeviceCoords;
//            text = null;
//            textPos = null;
//            textAlign = 0;
//        }
//        public DrawHint(final Shape shape, final Color color, final boolean isDeviceCoords) { this(shape, color, moreTransp(color), isDeviceCoords); }
//        public DrawHint(final Shape shape, final boolean isDeviceCoords) { this(shape, cool, isDeviceCoords); }
//        public DrawHint(final Shape shape) { this(shape, cool, false); }
//        public DrawHint(final String text, final Point2D textPos, final Color color, final boolean isDeviceCoords) {
//            this.shape = null;
//            this.colorFill = null;
//            this.colorLine = color;
//            this.isDeviceCoords = isDeviceCoords;
//            this.text = text;
//            this.textPos = textPos;
//            textAlign = 0;
//        }
}
