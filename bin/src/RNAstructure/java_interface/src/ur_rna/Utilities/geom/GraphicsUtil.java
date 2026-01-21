package ur_rna.Utilities.geom;

import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.annotation.Nullable;

import java.awt.*;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

/**
 * Utilities to help with {@link java.awt.Graphics2D } drawing and graphics.
 */
public class GraphicsUtil {
    public static Rectangle2D.Float calcTextBox(final Graphics2D g, final String text, final double centerX, final double centerY) {
        return calcTextBox(g, text, centerX, centerY, false, null);
    }
    public static Rectangle2D.Float calcTextBox(final Graphics2D g, final String text, final double centerX, final double centerY, boolean addAscentForDraw) {
        return calcTextBox(g, text, centerX, centerY, addAscentForDraw, null);
    }
    public static Rectangle2D.Float calcTextBox(final Graphics2D g, final String text, final double centerX, final double centerY, boolean addAscentForDraw, @Nullable ObjTools.RefInt ascent) {
        FontMetrics fm = g.getFontMetrics();
        Rectangle2D rcs = fm.getStringBounds(text, g);

        //
        //

        /*
        Java text metrics are complicated with respect to vertical measurements.
        e.g. with drawString(x, y) the y component is NOT the top of the text, but rather the baseline:
        Fonts are rendered on a baseline, running along the bottom of the text. The vertical space is allocated as follows:

        --- top
         ^
         |  leading
         |
         --
         ^              Y     Y
         |               Y   Y
         |                Y Y
         |  ascent         Y     y     y
         |                 Y      y   y
         |                 Y       y y
         __ baseline ______Y________y_________
         |                         y
         v  descent              yy
         -- bottom

        The leading is simply the font's recommended space between lines.
        For the sake of centering vertically between two points, you should ignore leading (it's ledding, BTW, not leeding; in general typography it is/was the lead spacing inserted between lines in a printing plate).

        So for centering the text ascenders and descenders, you want:

        baseline=(top+((bottom+1-top)/2) - ((ascent + descent)/2) + ascent;
        Without the final "+ ascent", you have the position for the top of the font; therefore adding the ascent goes from the top to the baseline.

        Also, note that the font height should include leading, but some fonts don't include it, and due to rounding differences, the font height may not exactly equal (leading + ascent + descent).
         */
        Rectangle2D.Float rc = new Rectangle2D.Float(
                (float)(centerX - rcs.getWidth() / 2), // calculate horizontal center
                (float)(centerY) - (1 + fm.getAscent() + fm.getDescent()) / 2f + (addAscentForDraw?fm.getAscent():0), // calculate vertical center (does NOT include ascent. it must be added before a call to drawString)
                (float)rcs.getWidth(), (float)rcs.getHeight());
        if (ascent != null) ascent.value = fm.getAscent();
        return rc;
    }
    public static Rectangle2D.Float drawTextCentered(final Graphics2D g, final String text, final Rectangle2D boundingBox) {
        // ----- Debug text layout -------------
//        g.setColor(Color.GREEN);
//        g.draw(boundingBox);
        // ----- ^ Debug text layout ^ -------------
        return drawTextCentered(g, text, boundingBox.getCenterX(), boundingBox.getCenterY());
    }
    public static Rectangle2D.Float drawTextCentered(final Graphics2D g, final String text, final Point2D centerPos) { return drawTextCentered(g, text, centerPos.getX(), centerPos.getY()); }
    public static Rectangle2D.Float drawTextCentered(final Graphics2D g, final String text, final double centerX, final double centerY) {
//        FontMetrics fm = g.getFontMetrics();
//
//        //g.setFont(symbol);
//        Rectangle2D rct = fm.getStringBounds(text, g);
//
//        float textX, textY, textH;
//        textX = (float) (centerX - rct.getWidth() / 2);
//        textH = (1 + fm.getAscent() + fm.getDescent());
//        textY = (float) (centerY - textH / 2f  + fm.getAscent());

        // ----- Debug text layout -------------
//        g.setColor(Color.RED);
//        g.setStroke(new BasicStroke(1));
//        //g.draw(rcs);
//        g.draw(new Rectangle2D.Double(centerX - rct.getWidth() / 2, centerY - rct.getHeight() / 2, rct.getWidth(), rct.getHeight()));
//        g.setStroke(new BasicStroke(0.25f));
//        g.setColor(Color.YELLOW);
//        for (int i = 0; i < fm.getAscent() + fm.getDescent() + 1; i++) {
//            g.draw(new Rectangle2D.Double(textX, textY + fm.getDescent() - i, rct.getWidth(), 1));
//        }
//        g.setColor(Color.BLUE);
//        g.draw(new Rectangle2D.Double(textX, textY, rct.getWidth(), 1));
//        g.draw(new Rectangle2D.Double(textX, textY - fm.getAscent(), rct.getWidth(), 1));
//        g.draw(new Rectangle2D.Double(textX, textY + fm.getDescent(), rct.getWidth(), 1));
//        g.setColor(Color.BLACK);
//        g.drawString("g", textX, textY);
        // ----- ^ Debug text layout ^ -------------

        ObjTools.RefInt ascent = new ObjTools.RefInt();
        Rectangle2D.Float rct = calcTextBox(g, text, centerX, centerY, false, ascent);

//        Color c = g.getColor();
//        g.setColor(Color.GREEN);
//        g.draw(rct);
//        g.setColor(c);

        g.drawString(text, rct.x, rct.y+ascent.value);
        return rct;
    }

}
