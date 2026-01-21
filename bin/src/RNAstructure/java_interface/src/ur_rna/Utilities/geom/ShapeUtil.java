package ur_rna.Utilities.geom;

import java.awt.*;
import java.awt.geom.*;

/**
 * Shape Utilities.
 */
public class ShapeUtil {
    public static Shape intersection(Shape s1, Shape s2) {
        if (s1 instanceof Rectangle2D && s2 instanceof Rectangle2D) {
            return ((Rectangle2D) s1).createIntersection((Rectangle2D) s2);
        } else {
            Area a1 = new Area(s1);
            a1.intersect(new Area(s2));
            return a1;
        }
    }

    public static Shape transform(Shape s, AffineTransform tx) {
        if (s == null)
            return null;

        if (tx == null || tx.isIdentity())
            return clone(s);

        boolean isRectlinearTr = (tx.getType() & (AffineTransform.TYPE_GENERAL_TRANSFORM | AffineTransform.TYPE_GENERAL_ROTATION)) == 0;
        if (s instanceof Rectangle2D && isRectlinearTr) {
            Rectangle2D rect = (Rectangle2D) s;
            double corners[] = new double[] { rect.getMinX(), rect.getMinY(), rect.getMaxX(), rect.getMaxY() };
            tx.transform(corners, 0, corners, 0, 2);
            rect = new Rectangle2D.Double();
            rect.setFrameFromDiagonal(corners[0], corners[1], corners[2], corners[3]);
            return rect;
        }
        return tx.createTransformedShape(s);
    }

    public static Shape inverseTransform(Shape s, AffineTransform tx) {
        if (s == null)
            return null;
        try {
            AffineTransform inverse = tx.createInverse();
            return transform(s, inverse);
        } catch (NoninvertibleTransformException e) {
            return null;
        }
    }

    @SuppressWarnings("unchecked") // clone() will always return the same type, so type cast error should never occur.
    public static Shape clone(Shape shape) {
        if (shape == null)
            return null;
        if (shape instanceof Line2D)
            return (Shape)((Line2D)shape).clone();
        if (shape instanceof Rectangle2D)
            return (Shape)((Rectangle2D)shape).clone();
        if (shape instanceof RoundRectangle2D)
            return (Shape)((RoundRectangle2D)shape).clone();
        if (shape instanceof Ellipse2D)
            return (Shape)((Ellipse2D)shape).clone();
        if (shape instanceof Arc2D)
            return (Shape)((Arc2D)shape).clone();
        if (shape instanceof Polygon) {
            Polygon p = (Polygon) shape;
            return new Polygon(p.xpoints, p.ypoints, p.npoints);
        }
        if (shape instanceof CubicCurve2D)
            return (Shape)((CubicCurve2D)shape).clone();
        if (shape instanceof QuadCurve2D)
            return (Shape)((QuadCurve2D)shape).clone();
        if (shape instanceof Path2D.Float)
            return (Shape)((Path2D)shape).clone();
        return new Path2D.Double(shape);
    }

//    /**
//     * Returns a hashCode for the state of the Shape.
//     * This can be used to detect changes to the shape.
//      * @param shape
//     * @return
//     */
//    public int stateHash(Shape shape) {
//
//    }

    public static boolean equals(Shape shapeA, Shape shapeB) {
        PathIterator pathAIterator = shapeA.getPathIterator(null);
        PathIterator pathBIterator = shapeB.getPathIterator(null);

        if (pathAIterator.getWindingRule() != pathBIterator.getWindingRule())
            return false;

        double[] pathASegment = new double[6];
        double[] pathBSegment = new double[6];
        while (!pathAIterator.isDone()) {
            int pathASegmentType = pathAIterator.currentSegment(pathASegment);
            int pathBSegmentType = pathBIterator.currentSegment(pathBSegment);
            if (pathASegmentType != pathBSegmentType)
                return false;
            for (int segmentIndex = 0; segmentIndex < pathASegment.length; segmentIndex++)
                if (pathASegment[segmentIndex] != pathBSegment[segmentIndex])
                    return false;
            pathAIterator.next();
            pathBIterator.next();
        }

        // if the shapes are equal then when shapeA is done shapeB is too.
        return pathBIterator.isDone();
    }

}
