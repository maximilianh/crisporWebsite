package ur_rna.Utilities.geom;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;

public class Vec2D extends Point2D.Double implements Cloneable {
    public static final double TwoPI = 2 * Math.PI;

    //public double x, y;
    public Vec2D(final double x, final double y) {
        this.x = x;
        this.y = y;
    }
    public Vec2D() { }
    public Vec2D(Point2D p) { x = p.getX(); y = p.getY(); }
    public Vec2D(Point2D from, Point2D to) { x = to.getX() - from.getX(); y = to.getY() - from.getY(); }
    public static Vec2D fromAngle(double theta) {
        return new Vec2D(Math.cos(theta), Math.sin(theta));
    }
    /** The magnitude (aka Euclidean length) of the vector.  */
    public double length() { return Math.sqrt(x*x + y*y); }
    public double lengthSqr() { return x*x + y*y; }
    //public Vec2D add(Vec2D v) { x += v.x; y += v.y; return this;}
    public Vec2D add(Point2D p) { x += p.getX(); y += p.getY(); return this;}
    public Vec2D add(double dx, double dy) { x += dx; y += dy; return this;}
    public Vec2D addScaled(Point2D p, double scale) { x += p.getX() * scale; y += p.getY()* scale; return this;}
    //public Vec2D subtr(Vec2D v) { x -= v.x; y -= v.y; return this;}
    public Vec2D subtr(Point2D p) { x -= p.getX(); y -= p.getY(); return this;}
    //public Vec2D multiply(Vec2D v) { x *= v.x; y *= v.y; return this;}
    public Vec2D multiply(Point2D p) { x *= p.getX(); y *= p.getY(); return this;}
    public double dot(Vec2D v) { return x * v.x + y * v.y; }
    public double dot(Point2D endPoint, Point2D referencePoint) { return x * (endPoint.getX() - referencePoint.getX()) + y * (endPoint.getY() - referencePoint.getY()); }
    public Vec2D negate() { x = -x; y = -y; return this;}
    public Vec2D normalize() { double mag = length(); x/=mag; y/=mag; return this;}
    public Vec2D normalize(double x_if_Zero, double y_if_Zero) {
        double mag = length();
        if (mag < 0.000000000001) {
            x = x_if_Zero;
            y = y_if_Zero;
        } else {
            x /= mag;
            y /= mag;
        }
        return this;
    }
    public double getAngle() { return Math.atan2(y, x); }
    // public Point2D toPoint() { return new Point2D.Double(x, y);}
    // public Vec2D setPoint(Point2D p) { p.setLocation(x, y); return this; }
    /**
     * Returns the location of result to be the point that would be reached by starting
     * at startPoint and traveling along the direction of this vector the distance specified by lengthAlongVector.
     * (The current length of this vector is irrelevant. I.e. It acts as if it were a unit vector.)
     * This is equivalent to {@code this.clone().setLength(lengthAlongVector).add(startPoint) }
     */
    public Vec2D getPointAlong(Point2D startPoint, double lengthAlongVector) { return this.clone().setLength(lengthAlongVector).add(startPoint); }
    /**
     * Sets the location of result to be the point that would be reached by starting
     * at startPoint and traveling along the direction of this vector the specified multiple of this vector's length.
     * This is equivalent to {@code this.clone().scale(lengthMultiple).add(startPoint) }
     * For example  @{code new Vec2D(p1, p2).setPointAlongScaled(p1, 0.5, p3) } would set the
     * point p3 to be the midpoint between p1 and p2.
     */
    public Vec2D getPointAlongScaled(Point2D startPoint, double lengthMultiple) {
        return this.clone().scale(lengthMultiple).add(startPoint);
    }
    public Point toPoint() { return new Point((int)x, (int)y);}
    public Point2D.Float toPointF() { return new Point2D.Float((float)x, (float)y);}
    /** Returns a line starting at (0,0) and ending at the tip of the vector -- (x, y). */
    public Line2D toLine() { return new Line2D.Double(0, 0, x, y);}
    /** Returns a line starting at the given point and ending at the starting point + the vector -- i.e. (start.x + x, start.y + y). */
    public Line2D toLine(Point2D start) { return new Line2D.Double(start.getX(), start.getY(), start.getX() + x, start.getY() + y);}
    /** Returns a line starting at the given point, in the direction of the vector, and with the specified distance. */
    public Line2D toLine(Point2D start, double distance) {
        double mag = length();
        return new Line2D.Double(start.getX(), start.getY(), start.getX() + x * distance / mag, start.getY() + y * distance / mag);
    }
    public Vec2D scale(double scale) { x *= scale; y *= scale; return this;}
    public Vec2D scale(double scaleX, double scaleY) { x *= scaleX; y *= scaleY; return this;}
    /**
     * Sets the length (magnitude) of this vector.
     * If newLength is negative, the direction of the vector is also reversed.
     * This is equivalent to {@link #normalize()} followed by {@link #scale(double)}
     */
    public Vec2D setLength(double newLength) { scale(newLength / length()); return this;}
    public Vec2D transform(AffineTransform tr) {
        double px = x * tr.getScaleX() + y * tr.getShearX() + tr.getTranslateX();
        double py = y * tr.getScaleY() + x * tr.getShearY() + tr.getTranslateY();
        x = px; y = py;
        return this;
    }
    /**
     * Rotates this vector 90 degrees counter-clockwise.  (x is replaced with -y and y is replaced with x).
     * The rotated vector is perpendicular (orthogonal) to its previous state.
     */
    public Vec2D rotate90() { return setTo(-y, x); }
    public Vec2D rotate(final double theta) {
        double cosT = Math.cos(theta);
        double sinT = Math.sin(theta);
        // x' = x cos theta - y sin theta
        // y' = x sin theta + y cos theta
        double xp = x;
        x = x * cosT - y * sinT;
        y = xp * sinT + y * cosT;
        return this;
    }


    /** Alias of {@link #setLocation(double, double)} */
    public Vec2D setTo(double newX, double newY) {x = newX; y = newY; return this; }
    /** Alias of {@link #setLocation(Point2D)} */
    public Vec2D setTo(Point2D p) {x = p.getX(); y = p.getY(); return this; }
    public Vec2D setTo(Vec2D p) {x = p.x; y = p.y; return this; }

    public Vec2D getMidpoint(Point2D other) { return new Vec2D((x + other.getX())/2, (y + other.getY())/2); }
    public static Vec2D getMidpoint(Point2D start, Point2D end) { return new Vec2D((start.getX() + end.getX())/2, (start.getY() + end.getY())/2); }
    public static Vec2D getPointBetween(Point2D start, Point2D end, double fractionOfDistance) { return new Vec2D(start.getX() + (end.getX() - start.getX())*fractionOfDistance, start.getY() + (end.getY() - start.getY())*fractionOfDistance); }
    /**
     * Returns the distance between two points.
     * @param p1 one of the points
     * @param p2 one of the points
     * @return the distance between the two points
     */
    public static double distance(Point2D p1, Point2D p2) {
        double x = p1.getX() - p2.getX();
        double y = p1.getY() - p2.getY();
        return Math.sqrt(x * x + y * y);
    }
    /**
     * Returns the square of the distance between two points, i.e.  {@code (p1.x - p2.x2)^2  + (p1.y - p2.y)^2 }
     * @param p1 one of the points
     * @param p2 one of the points
     * @return the square of he distance between the two points
     */
    public static double distanceSq(Point2D p1, Point2D p2) {
        double x = p1.getX() - p2.getX();
        double y = p1.getY() - p2.getY();
        return x * x + y * y;
    }

    //public Vec2D getAdded(Vec2D v) { return this.clone().add(v); }
    public Vec2D getAdded(Point2D p) { return this.clone().add(p); }
    //public Vec2D getSubtr(Vec2D v) { return this.clone().subtr(v); }
    public Vec2D getSubtr(Point2D p) { return this.clone().subtr(p); }
    //public Vec2D getMultiplied(Vec2D v) { return this.clone().multiply(v); }
    public Vec2D getMultiplied(Point2D p){ return this.clone().multiply(p); }
    public Vec2D getNegated() { return this.clone().negate();}
    public Vec2D getNormalized() { return this.clone().normalize(); }
    public Vec2D getScaled(double scale) { return this.clone().scale(scale); }
    public Vec2D getScaled(double scaleX, double scaleY) { return this.clone().scale(scaleX, scaleY); }
    public Vec2D getTransformed(AffineTransform tr) { return this.clone().transform(tr); }
    /** Get a copy of this vector rotated counter-clockwise by 90 degrees */
    public Vec2D getRotated90() { return new Vec2D(-y, x); }
    /** Get a copy of this vector rotated counter-clockwise by the specified angle (in radians) */
    public Vec2D getRotated(double theta) { return this.clone().rotate(theta); }

    public Vec2D clone() {
//        try {
            return (Vec2D)super.clone();
//        }catch (CloneNotSupportedException ex) {
//            throw new InternalError(ex);
//        }
    }
    /** Gets the smallest angle between two vectors. The result is always between 0 and PI (inclusive). Based on the dot product {@code A · B = |A| * |B| * cos(θ)} */
    public double angleBetween(final Vec2D v2) {
        return Math.acos(dot(v2) / Math.sqrt(lengthSqr() * v2.lengthSqr()));
    }
    /** Gets the full angle from this vector to v2, sweeping counter-clockwise. The result is from -PI to PI. The order of the vectors is important: {@code v1.angleTo(v2) == -v2.angleTo(v1)}  */
    public double angleTo(final Vec2D v2) {
        // The dot product returns the smallest angle between the two vectors,
        // so the result is always between 0 and PI.

        // But if we want to determine the full angle moving counter clockwise from v1 to v2,
        // we have to compare the angle to each vector.
        double dot = x*v2.x+y*v2.y; // dot product = cos(a) * |v1||v2|
        double det = x*v2.y-y*v2.x; // determinant = sin(a) * |v1||v2|  (same as length of cross product)
        return Math.atan2(det, dot); // i.e. atan(r sin(a), r cos(a)) where r=|v1||v2| is the magnitude an doesn't affect the angle.
    }
    /**
     * This calculates the angle from one vector to another, similar to {@link #angleTo(Vec2D)} except that the
     * output is intended for applications in which the direction of travel is indicated by the sign of the result
     * and the extent of travel is indicated by its absolute value. Thus instead of returning a value from -PI to PI, this
     * function returns a value from -2PI to 2PI.
     *
     * For example, want to draw an arc from 0 degrees to 270 degrees counter-clockwise. (I.e. the arc should span 3 quadrants.)
     * So v1 = (1, 0) and v2 = (0, -1)
     * But if {@link #angleTo(Vec2D)} were used, it would return -90 degrees, which would result in the arc being drawn CLOCKWISE through only the bottom-right quadrant.
     * This function would return 270 degrees (reverse = false) or -90 (reverse = true);
     *
     * As another example, assume v1 points to (-1, 0) and v2 points to (1, 0). This function would return 180 (reverse = false) or -180 (reverse = true);
     * As another example, assume v1 points to (1, 1) and v2 points to (-1, 1). This function would return 90 (reverse = false) or -270 (reverse = true);
     *
      * @param v2 the other vector
     * @param reverse whether the arc should be drawn counter-clockwise (reverse=false) or clockwise (reverse=true) from v1 to v2.
     *                         Alternatively reverse=true could be interpreted as "draw the arc counter-clockwise from v2 to v1 instead of from v1 to v2"
     *                If reverse is true, the output will be in the range (-360, 0]
     *                If reverse is false, the output will be in the range [0, 360)
     * @return
     */
    public double arcAngleTo(final Vec2D v2, boolean reverse) {
        double angle = angleTo(v2);
        if (reverse && angle > 0)
            angle -= TwoPI;
        else if (!reverse && angle < 0)
            angle += TwoPI;
        return angle;
    }

    /**
     * Given any angle in radians, this returns an equivalent angle in the range [-PI, PI]
     * In general, odd multiples of PI map to -PI if they are negative and to PI if they are positive.
     * (Even multiples of PI, e.g. 2PI, of course, map to 0)
     */
    public static double wrapToPI(double angle) {
        return angle < 0 ?
            (angle+Math.PI) % TwoPI - Math.PI : // returns -PI if angle is an exact multiple of PI
            (angle-Math.PI) % TwoPI + Math.PI;  // returns PI if angle is an exact multiple of PI
    }

    /**
     * Given any angle in radians, this returns an equivalent angle in the range [0, 2*PI]
     * In general, multiples of 2*PI map to 0 if they are negative and to 2*PI if they are positive.
     */
    public static double wrapTo2PI(double angle) {
        return wrapToPI(angle)+Math.PI;
    }

    /**
     * Given any angle in radians, this returns an equivalent angle in the range [0, 2*PI )
     * Note that 2*PI is an exclusive upper limit (unlike wrapTo2PI).
     * In general, all multiples of 2*PI map to 0.
     */
    public static double mod2PI(double angle) {
        return (angle = angle % TwoPI) < 0 ? angle + TwoPI : angle;
    }

    /**
     * Calculates the projection of this vector onto the vector specified by direction.
     * If direction is a unit vector this is equal to {@code direction * (this·direction)}
     * More generally, this is equal to {@code direction.normalize().scale(this·direction/direction.length())}
     */
    public Vec2D projectedOn(final Vec2D direction) {

        // Skip normalization and just divide by direction.lengthSqr()
        // This is equivalent to direction.clone().normalize().scale(this.dot(direction)/direction.length())
        return direction.clone().scale(dot(direction) / direction.lengthSqr());
    }

    private static final double toStringRoundFactor = 10000000;

    @Override
    public String toString() {
        return "[" + this.getClass().getSimpleName() +
                '(' + x + ", " + y + ')' +
                " L: " + Math.round(this.length() * toStringRoundFactor) / toStringRoundFactor +
                " θ: " + Math.round(this.getAngle() * 180 / Math.PI * toStringRoundFactor) / toStringRoundFactor +
                ']';
    }
    public static Vec2D reflect(Point2D pt, Vec2D reflectionVector, Point2D pointOnReflectionLine) {
        // Define R = reflectionVector, p0 = pointOnReflectionLine, p = pt, pR = result (reflected)
        // Let P = p0 - p (the ray from p to p0) which points from p towards R (at an unspecified angle)
        // Then P has two components: one parallel to R (L) and one perpendicular to it (K), (i.e. P=L+K)
        // The result we want is p + 2K = p + 2*(P - L)
        // L is simply the projection of d onto R, i.e. R*(R dot d)
        // So pR = p + 2*(P - R*(R dot d))    =    p + 2(p0 - p) - 2R*(R dot d)   =   -p + 2*p0 + 2*R*(d.R)
        Vec2D P = new Vec2D(pt, pointOnReflectionLine); // ray from arbitrary point on the line to the subject point.
        Vec2D K = P.projectedOn(reflectionVector).scale(-2); // get projection of P onto R. Scale to -2K
        P.scale(2).add(pt); // this gives P = p + 2(p0 - p)
        return P.add(K);
    }
}
