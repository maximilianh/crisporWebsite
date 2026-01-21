package ur_rna.Utilities.geom;

import java.awt.*;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

public abstract class Rectangles {
    public static Rectangle fromPoints(Point p1, Point p2) {
        int w = Math.abs(p1.x - p2.x);
        int h = Math.abs(p1.y - p2.y);
        int x = Math.min(p1.x, p2.x);
        int y = Math.min(p1.y, p2.y);
        return new Rectangle(x, y, w, h);
    }
    public static Rectangle2D.Double fromPoints(Point2D p1, Point2D p2) {
        double w = Math.abs(p1.getX() - p2.getX()),
               h = Math.abs(p1.getY() - p2.getY()),
               x = Math.min(p1.getX(), p2.getX()),
               y = Math.min(p1.getY(), p2.getY());
        return new Rectangle2D.Double(x, y, w, h);
    }
    public static Rectangle2D.Double fromPoints(Vec2D p1, Vec2D p2) {
        return new Rectangle2D.Double(Math.min(p1.x, p2.x), Math.min(p1.y, p2.y),  Math.abs(p1.x - p2.x), Math.abs(p1.y - p2.y));
    }
    public static Rectangle2D.Double fromCenter(Point2D p, double w, double h) {
        return  new Rectangle2D.Double(p.getX() - w /2, p.getY() - h / 2, w, h);
    }
    public static Rectangle2D.Double fromCenter(Vec2D p, double w, double h) {
        return  new Rectangle2D.Double(p.x - w /2, p.y - h / 2, w, h);
    }
    public static Rectangle fromCenter(Point p, int w, int h) {
        return  new Rectangle(p.x - w /2, p.y - h / 2, w, h);
    }
    public static Rectangle2D.Double fromCenter(Point2D p, double squareSize) {
        return  new Rectangle2D.Double(p.getX() - squareSize /2, p.getY() - squareSize / 2, squareSize, squareSize);
    }
    public static Rectangle2D.Double fromCenter(Vec2D p, double squareSize) {
        return  new Rectangle2D.Double(p.x - squareSize /2, p.y - squareSize / 2, squareSize, squareSize);
    }
    public static Rectangle fromCenter(Point p, int squareSize) {
        return  new Rectangle(p.x - squareSize /2, p.y - squareSize / 2, squareSize, squareSize);
    }

    public static void grow(Rectangle2D rc, double growBothDims) {grow(rc, growBothDims, growBothDims); }
    public static void grow(Rectangle2D rc, double growX, double growY) {
        rc.setFrame(rc.getX() - growX, rc.getY() - growY, rc.getWidth() + 2 * growX, rc.getHeight() + 2 * growY);
    }
}
