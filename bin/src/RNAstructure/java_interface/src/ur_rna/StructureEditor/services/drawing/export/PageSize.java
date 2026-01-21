package ur_rna.StructureEditor.services.drawing.export;

import java.awt.geom.Dimension2D;

/**
 * Represents the dimensions of a fixed-size page (e.g. for SVG or PS etc)
 */
public class PageSize extends Dimension2D {
    public static final double MM_PER_INCH = 25.4;

    public static final PageSize LETTER = new PageSize(8.5*MM_PER_INCH, 11*MM_PER_INCH);
    public static final PageSize LEGAL = new PageSize(8.5*MM_PER_INCH, 14*MM_PER_INCH);
    public static final PageSize A0 = new PageSize(841 , 1189);
    public static final PageSize A1 = new PageSize(594 , 841);
    public static final PageSize A2 = new PageSize(420 , 594);
    public static final PageSize A3 = new PageSize(297 , 420);
    public static final PageSize A4 = new PageSize(210 , 297);
    public static final PageSize A5 = new PageSize(148 , 210);
    public static final PageSize A6 = new PageSize(105 , 148);
    public static final PageSize TABLOID = new PageSize(11*MM_PER_INCH, 17*MM_PER_INCH);
    public static final PageSize LEDGER = TABLOID.getLandscape();

    /** Width of the page in mm */
    public final double width;
    /** Height of the page in mm */
    public final double height;

    public PageSize(double width, double height) {
        this.width = width;
        this.height = height;
    }

    public PageSize(Dimension2D size) {
        this(size.getWidth(), size.getHeight());
    }

    public PageSize getPortrait() {
        if (width <= height) {
            return this;
        }
        //return new PageSize(x, y, height, width);
        return new PageSize(height, width); //swap height & width
    }

    public PageSize getLandscape() {
        if (width >= height) {
            return this;
        }
        //return new PageSize(x, y, height, width);
        return new PageSize(height, width); //swap height & width
    }

    @Override
    public double getWidth() {
        return width;
    }

    @Override
    public double getHeight() {
        return height;
    }

    @Override
    public void setSize(final double width, final double height) {
        throw new UnsupportedOperationException("PageSize dimensions are fixed.");
    }
}
