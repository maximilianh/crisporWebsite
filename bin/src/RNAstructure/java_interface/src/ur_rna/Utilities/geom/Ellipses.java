package ur_rna.Utilities.geom;

import java.awt.geom.Arc2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

public class Ellipses {
    public static Ellipse2D.Double fromRect(Rectangle2D rc) {
        Ellipse2D.Double e = new Ellipse2D.Double();
        e.setFrame(rc);
        return e;
    }
    public static Ellipse2D.Double fromCenter(Point2D p, double w, double h) {
        return fromCenter(p.getX(), p.getY(), w, h);
    }
    public static Ellipse2D.Double fromCenter(Vec2D p, double w, double h) {
        return fromCenter(p.x, p.y, w, h);
    }
    public static Ellipse2D.Double fromCenter(double x, double y, double w, double h) {
        return new Ellipse2D.Double(x - w / 2, y - h / 2, w, h);
    }
    public static Ellipse2D.Double fromCenter(Point2D p, double radius) {
        return fromCenter(p.getX(), p.getY(), radius);
    }
    public static Ellipse2D.Double fromCenter(Vec2D p, double radius) {
        return fromCenter(p.x, p.y, radius);
    }
    public static Ellipse2D.Double fromCenter(double x, double y, double radius) {
        return new Ellipse2D.Double(x - radius, y - radius, 2 * radius, 2 * radius);
    }

    public static Arc2D.Double arc(Point2D center, double radius, double startAngle, double arcAngle, int arcType) {
        return new Arc2D.Double(center.getX() - radius, center.getY() - radius, radius * 2, radius * 2, startAngle * toScreenDegrees, arcAngle * toScreenDegrees, arcType);
    }
    /**
     * The normal Arc2D constructor expects the angles in degrees.
     * Additionally it
     *
     * @param center     the center of the arc's circle
     * @param radius     the radius of the arc's circle
     * @param startAngle the starting angle in radians. This will be multiplied by {@link #toScreenDegrees} to account for the upside-down screen coordinates and convert from radians to degrees.
     * @param endAngle   the end angle in radians. This will be multiplied by {@link #toScreenDegrees} to account for the upside-down screen coordinates and convert from radians to degrees.
     * @param reversed   The direction in which the arc should be drawn clockwise(reversed=true) or counter-clockwise(reversed=false).
     * @param arcType    The type of arc to draw: {@link Arc2D#CHORD}, {@link Arc2D#PIE}, or {@link Arc2D#OPEN}
     */
    public static Arc2D.Double arc(Point2D center, double radius, double startAngle, double endAngle, boolean reversed, int arcType) {
        return arc(center, radius, startAngle, calcArcExtent(startAngle, endAngle, reversed), arcType);
    }

    // Multiply an angle in model-space radians by this to get degrees in screen-space (where the y-coordinate is upside-down)
    public final static double toScreenDegrees = -180 / Math.PI;

    public static double calcArcExtent(double startAngle, double endAngle, boolean reversed) {
        double extent = (endAngle - startAngle) % Math.PI;
        if (!reversed && extent < 0)
            extent += 2 * Math.PI;
        else if (reversed && extent > 0)
            extent -= 2 * Math.PI;
        return extent;
    }

    public static class Circle extends Ellipse2D.Float {
        private float _radius = 1;
        private Point2D.Float _center = new Point2D.Float();

        public Circle(Point2D center, double radius) { this(center.getX(), center.getY(), radius); }
        public Circle(double radius) { this(0, 0, radius); }
        public Circle(double centerX, double centerY, double radius) {
            this.radius(radius);
            this.center(centerX, centerY);
        }
        /**
         * This calculates the circumference of a circle based on an arc that is cut by a chord.
         * The arc-length and the chord length are known, but the full circumference is not (and neither is the sweep angle or the radius).
         * There is no analytical solution to the problem of how to determine the full circumference of a circle from a partial arclength and the chord length, so just estimate it roughly here.
         */
        public static double getRadiusFromArcAndChord(double arcLength, double chordLength) {
            // Let `theta` be HALF of the sweep angle of the chord.
            // (Note that valid values are 0 <= theta < 180 degrees, because if theta is 180, then the full sweep angle would be 360, which means arcLength is 0 and the circumference is therefore infinite)
            // Then:
            // S = chord length = 2 r sin theta
            // K = arc-length = full_circumference - r * 2 theta = 2 PI r - r * 2 theta  = 2r (PI - theta)
            // If we can determine theta, then we can obtain the radius from:
            // r = S / (2 sin theta)
            // so C = full_circumference = K + r * 2 theta
            double ratio = chordLength / arcLength; // calculate S/K, which is equal to (sin theta)/(PI-theta) ... independent of the radius of the circle.

    //        if (ratio < 2 / Math.PI) // ratio = 2/pi = 0.636619772 corresponds to theta = 90. (The full sweep angle is 2*theta=180) The chord BISECTS the circle.  (e.g. assume the radius is 1. The arc length (K) is PI and the chord length (S) is 2. S/K = 2/PI
    //            return (chordLength * Math.cosh(ratio * Math.PI / 2) + arcLength) / 2 / Math.PI; // cosh is just a rough estimate (which I discovered) that approximates the correct result with the error in the circumference less than 1% for all S/K <= 2/PI

            // -0.1549x2 + 1.0574x + 0.0975
            double sinHalfTheta = -0.1549*ratio*ratio + 1.0574*ratio + 0.0975; // determined by graph of [sin(theta/2) vs S/K] fit with an polynomial (order 2) curve. Theta/2 was used instead of Theta to make the curve closer to linear.
            double theta = 2 * Math.asin(sinHalfTheta);
            double r = chordLength / 2 / Math.sin(theta);
            // The next part seems redundant, because we already have r.
            // But in fact the equation is NOT accurate for the intermediate r (or intermediate theta),
            // but it IS VERY accurate (err < 0.16 %) for the resulting circumference
            // (and the corresponding radius and theta derived from it)
            // The reason is simply that even though neither r nor theta are accurate, the product r*2* Math.asin(sinTheta) is accurate for the missing segment of the circumference.
            double circ = arcLength + r * 2 * theta; // full circumference
            return circ / (2 * Math.PI);

    //        if (ratio < 0.99468492) { // this is just an arbitrary number, chosen such that the maximum error in the following formula is less than 1%. It corresponds to theta = 169.76 (full sweep angle = 339.52)
    //            // The accuracy of this equation (which I discovered) is essentially 100% until theta = 150 degrees.
    //            double sinTheta =-0.1549*ratio*ratio + 1.0574*ratio + 0.0975; // determined by graph of [sin(theta) vs S/K] fit with an polynomial (order 2) curve.
    //            return chordLength / 2 / sinTheta;
    //        }
            // when theta is greater than 339.52, the above formula quickly loses accuracy, and fitting becomes difficult. Very small changes in theta result in huge changes to the radius, because as the chord moves closer to the "top" of the circle, the circumference has to increase dramatically to accommodate the chord.
            // at this point, the most accurate method is to simply search for a radius that works. This is equivalent mathematically to searching for an angle, however due to machine precision, it is best to use the radius (or circumference) as the independent variable and back calculate the angle from it.

    //        // do a "binary" search to discover the radius.
    //        double prev = chordLength / 1.99201995; // this number is 2*sin(theta), where theta is 169.76 degrees, which corresponds to the
    //        while(true) {
    //            double next = prev * 2; // radius
    //
    //            double foundLen = arcLength * Math.sin(theta)/(Math.PI  - theta);
    //            if (Math.abs(chordLength - foundLen)/chordLength <= 0.1)
    //                return arcLength / (2 * Math.PI - 2 * theta);
    //
    //        } while ();
        }

        //public Ellipse2D ellipse() { return Ellipses.fromCenter(center, radius); }
        public Arc2D arc(final double startAngle, final double arcAngle) {
            return Ellipses.arc(_center, _radius, startAngle, arcAngle, Arc2D.OPEN);
        }

        /** get the x-coordinate of the circle at the specified angle. */
        public double getCircleX(final double angle) {
            return _center.x + _radius * Math.cos(angle);
        }

        /** get the y-coordinate of the circle at the specified angle. */
        public double getCircleY(final double angle) {
            return _center.y + _radius * Math.sin(angle);
        }

        public void getPoint(final double angle, Point2D result) {
            result.setLocation(getCircleX(angle), getCircleY(angle));
        }
        public Point2D getPoint(final double angle) {
            return new Point2D.Double(getCircleX(angle), getCircleY(angle));
        }

//        @Override
//        public void setFrame(final float x, final float y, final float w, final float h) {
//            super.setFrame(x, y, w, h);
//            update();
//        }
//        /**
//         * {@inheritDoc}
//         *
//         * @param x
//         * @param y
//         * @param w
//         * @param h
//         * @since 1.2
//         */
//        @Override
//        public void setFrame(final double x, final double y, final double w, final double h) {
//            super.setFrame(x, y, w, h);
//            update();
//        }

//        private void update() {
//            //ptCenter.setLocation(getCenterX(), getCenterY());
//            //radius = (float)getWidth() / 2;
//        }

        public void translate(float dx, float dy) {
            this.x += dx;
            this.y += dy;
            _center.x += dx;
            _center.y += dy;
        }
        public float radius() {
            return _radius;
        }
        public void radius(double radius) {
            this._radius = (float) radius;
            this.width = this.height = (float) (radius * 2);
            this.x = _center.x - this._radius;
            this.y = _center.y - this._radius;
        }
        public Point2D.Float center() {
            return _center;
        }
        public void center(Point2D center) { this.center(center.getX(), center.getY()); }
        public void center(double x, double y) {
            this._center.setLocation(x, y);
            this.x = (float) (x - this._radius);
            this.y = (float) (y - this._radius);
        }
        public void setCircle(final Point2D center, final float radius) {
            this.center(center);
            this.radius(radius);
        }
        /**
         * Calculate the center and radius of a circle from from 3 points on the circle.
         *
         * @param p1 A point on the circle
         * @param p2 A point on the circle
         * @param p3 A point on the circle
         * @return A Circle that contains the three given points or NULL if the circle cannot be calculated
         * (e.g. because two or more points are identical or because all three points are co-linear).
         */
        public static Circle calcFromPoints(Point2D p1, Point2D p2, Point2D p3) {
            // Define equation of circle from radius (r) and center (xc, yc)
            //    (x-xc)^2 + (y-yc)^2 = r^2
            // Expand equation and rearrange to formulate into a LINEAR equation in 3 variables (ri, xi, yi)
            //    x^2 + y^2 = (2 xc)x + (2 yc)y + (r^2 - xc^2 - yc^2);
            //        ri    =   a*xi  +   b*yi   +        c
            // c = r^2 - xc^2 - yc^2   ==>  r = sqrt(c + xc^2 + yc^2)
            double x1 = p1.getX(), y1 = p1.getY();
            double x2 = p2.getX(), y2 = p2.getY();
            double x3 = p3.getX(), y3 = p3.getY();
            Matrix3D m = new Matrix3D(
                    x1, y1, 1,
                    x2, y2, 1,
                    x3, y3, 1
            );
            // determine the result vector: [ri] = [r1, r2, r3] = [ x1^2+y1^2, x2^2+y2^2, x3^2+y3^2 ]
            double[] results = new double[]{x1 * x1 + y1 * y1, x2 * x2 + y2 * y2, x3 * x3 + y3 * y3}; // xi^2 + yi^2
            if (!m.invert()) return null;
            m.transform(results, 0, 0);
            double xc = results[0] / 2, yc = results[1] / 2;
            return new Circle(xc, yc, Math.sqrt(results[2] + xc * xc + yc * yc));
        }

        public static Circle calcLeastSquaresFit(Point2D... points) {
            /*
             *   If you have (x, y) data distributed in a ring-shape on the xy-plane, you can use least squares
             *   regression to find the equation of the circle that best fits the data. That is, you determine the
             *   values of the radius (r) and center (h, k) such that the curve
             *
             *   (x - h)^2 + (y - k)^2 = r^2
             *
             *   provides a good fit around the data points.
             *   With least squares, "best fit" means that you minimize the equation
             *
             *   F(h, k, r) = ∑[(xi - h)^2 + (yi - k)^2 - r^2]^2
             *
             *   by solving the system ∂F/∂h = 0, ∂F/∂k = 0, and ∂F/∂r = 0.
             *   The equation of the circle can be linearized as follows:
             *
             *   (x - h)^2 + (y - k)^2 = r^2
             *   x^2 - 2hx + h^2 + y^2 - 2ky + k^2 = r^2
             *   x^2 + y^2 = 2hx + 2ky + r^2 - h^2 - k^2
             *   x^2 + y^2 =  Ax +  By +       C
             *
             *   Where A=2h, B=2k and C=(r^2 - h^2 - k^2)
             *   r = sqrt(C+h^2+k^2)
             *
             *   The matrix equation for circular regression is
             *
             *   [ ∑(xi^2),  ∑(xi.yi), ∑(xi) ]   [A]    [ ∑(xi.(xi^2+yi^2) ]
             *   [ ∑(xi.yi), ∑(yi^2),  ∑(yi) ] * [B] =  [ ∑(yi.(xi^2+yi^2) ]
             *   [ ∑(xi),    ∑(yi),      N   ]   [C]    [ ∑(xi^2+yi^2)     ]   (where N is the number of data points)
             *
             *   This will be represented by the following variable names:
             *           mx                      vec
             *   -----------------          -------------
             *   [ sx2, sxy, sx  ]   [A]    [ sr1       ]
             *   [ sxy, sy2, sy  ] * [B] =  [ sr2       ]
             *   [ sx,  sy,  p.L ]   [C]    [ sx2 + sy2 ]  ( because ∑(xi^2+yi^2) = ∑(xi^2) + ∑(yi^2) )
             *
             *   If the 3-by-3 matrix on the left is invertible, then there is a unique set of values for A, B, and C
             *   that generates the circle of best fit.
             *   i.e. mx.invert() * vec = [ A,B,C ]
             *   Note that multiplication of a vector by a matrix is represented by mx.transform(vec).
             *
             *   Once you have values for A, B, and C, you can compute h, k, and r with these transformations:
             *   h = A/2
             *   k = B/2
             *   r = sqrt(4C + A^2 + B^2) / 2   or more simply,  sqrt(C + h^2 + k^2)
             */
            if (points.length < 3)
                throw new IllegalArgumentException("At least 3 points are required to determine a circle.");
            if (points.length == 3)
                return calcFromPoints(points[0], points[1], points[2]);
            //http://www.had2know.com/academics/best-fit-circle-least-squares.html

            //  x^2 + y^2 = (2 xc)x + (2 yc)y + (r^2 - xc^2 - yc^2);

            double sx2 = 0, sy2 = 0, sxy = 0, sx = 0, sy = 0, sr1 = 0, sr2 = 0;
            for (Point2D p : points) {
                double x = p.getX(), y = p.getY(); // xi and yi
                double x2 = x * x,   y2 = y * y;   // xi^2 and yi^2
                sx2 += x2;     // ∑(xi^2)
                sy2 += y2;     // ∑(yi^2)
                sxy += x * y;  // ∑(xi.yi)
                sx += x;       // ∑(xi)
                sy += y;       // ∑(yi)
                sr1 += x * (x2 + y2);  // ∑(xi.(xi^2+yi^2)
                sr2 += y * (x2 + y2);  // ∑(yi.(xi^2+yi^2)
            }
            Matrix3D mx = new Matrix3D(
                    sx2, sxy, sx,
                    sxy, sy2, sy,
                    sx, sy, points.length);

            if (!mx.invert()) return null; // a NON-invertible matrix, so just return null. no circle can be calc'd
            double[] vec = new double[]{sr1, sr2, sx2 + sy2};
            mx.transform(vec);
            double cx = vec[0] / 2;
            double cy = vec[1] / 2;
            double r = Math.sqrt(vec[2] + cx * cx + cy * cy);
            return new Circle(cx, cy, r);
        }
    }
}
