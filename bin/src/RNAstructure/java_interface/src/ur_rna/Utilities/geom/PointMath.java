package ur_rna.Utilities.geom;

import java.awt.geom.Point2D;

/**
 * @author Richard M. Watson
 */
public class PointMath {
    /** sets result to p1 + p2.  If result is null, a new Point2D (cloned from p1) will be returned. */
    public static Point2D sum(Point2D p1, Point2D p2, Point2D result) {
        return sum(p1, p2.getX(), p2.getY(), result);
    }
    public static Point2D sum(Point2D p1, Point2D p2) { return sum(p1, p2, null); }

    /** sets result to (p.x + dx, p.y + dy)  If result is null, a new Point2D (cloned from p1) will be returned. */
    public static Point2D sum(Point2D p, double dx, double dy, Point2D result) {
        if (result == null)
            result = (Point2D) p.clone();
        result.setLocation(p.getX() + dx, p.getY() + dy);
        return result;
    }
    public static Point2D sum(Point2D p, double dx, double dy) { return sum(p, dx, dy, null); }

    /** sets result to p1 - p2.  If result is null, a new Point2D (cloned from p1) will be returned. */
    public static Point2D diff(Point2D p1, Point2D p2, Point2D result) {
        return sum(p1, -p2.getX(), -p2.getY(), result);
    }
    public static Point2D diff(Point2D p1, Point2D p2) { return diff(p1, p2, null); }

    /** Returns the mid-point between two points. */
    public static Point2D midpoint(Point2D p1, Point2D p2, Point2D result) {
        if (result == null)
            result = (Point2D)p1.clone();
        result.setLocation((p1.getX() + p2.getX()) / 2, (p1.getY() + p2.getY()) / 2);
        return result;
    }
    public static Point2D midpoint(Point2D p1, Point2D p2) { return midpoint(p1, p2, null); }

    /** sets result to -p.  If result is null, a new Point2D (cloned from p) will be returned. */
    public static Point2D negative(Point2D p, Point2D result) {
        if (result == null)
            result = (Point2D) p.clone();
        result.setLocation(-p.getX(), -p.getY());
        return result;
    }
    public static Point2D negative(Point2D p) { return negative(p, null); }
    /**
     * Assumes p is a vector and sets result to the corresponding normalized unit vector (i.e. divides x and y by the magnitude of vector (x, y).
     * If result is null, a new Point2D (cloned from p) will be returned.
     */
    public static Point2D normalized(Point2D p, Point2D result) {
        if (result == null)
            result = (Point2D) p.clone();
        double x = p.getX(), y = p.getY();
        double magnitude = Math.sqrt(x * x + y * y);
        result.setLocation(x / magnitude, y / magnitude);
        return result;
    }
    /** Returns the absolute length of the vector pointing to p */
    public static double magnitude(Point2D p) {
        double x = p.getX(), y = p.getY();
        return Math.sqrt(x * x + y * y);
    }
    public static Point2D normalized(Point2D p) { return normalized(p, null); }

    /** sets result to scale * p.  If result is null, a new Point2D (cloned from p) will be returned. */
    public static Point2D scaled(Point2D p, double scale, Point2D result) {
        return scaled(p, scale, scale, result);
    }
    public static Point2D scaled(Point2D p, double scale) { return scaled(p, scale, null); }

    public static Point2D scaled(Point2D p, double sx, double sy, Point2D result) {
        if (result == null)
            result = (Point2D) p.clone();
        result.setLocation(p.getX() * sx, p.getY() * sy);
        return result;
    }
    public static Point2D scaled(Point2D p, double sx, double sy) { return scaled(p, sx, sy, null); }

    // ************************************************************************************
    // The following all mutate the given point instead of setting the result in a new point:
    // ************************************************************************************

    /** translates p by (p.x, p.y). Returns p itself (for chaining). */
    public static Point2D translate(Point2D p, Point2D delta) {
        translate(p, delta.getX(), delta.getY());
        return p;
    }

    /** ranslates p by (dx, dy). Returns p itself (for chaining). */
    public static Point2D translate(Point2D p, double dx, double dy) {
        p.setLocation(p.getX() + dx, p.getY() + dy);
        return p;
    }

    /** sets result to -p.  Returns p itself (for chaining). */
    public static Point2D negate(Point2D p) {
        p.setLocation(-p.getX(), -p.getY());
        return p;
    }
    /**
     * Assumes p represents a vector and mormalizes it to a unit vector -- i.e. divides x and y by the magnitude of vector (x, y).
     */
    public static Point2D normalize(Point2D p) {
        normalized(p, p);
        return p;
    }

    /** scales p by the given factor.  Returns p itself (for chaining). */
    public static Point2D scale(Point2D p, double factor) {
        return scale(p, factor, factor);
    }

    public static Point2D scale(Point2D p, double sx, double sy) {
        p.setLocation(p.getX() * sx, p.getY() * sy);
        return p;
    }
    public static double dist(final Point2D p1, final Point2D p2) {
        return Math.sqrt(distSqr(p1, p2));
    }
    public static double distSqr(final Point2D p1, final Point2D p2) {
        return Point2D.distanceSq(p1.getX(), p1.getY(), p2.getX(), p2.getY());
    }
}
