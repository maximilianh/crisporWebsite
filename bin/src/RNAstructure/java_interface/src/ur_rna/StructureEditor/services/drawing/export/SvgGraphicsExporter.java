package ur_rna.StructureEditor.services.drawing.export;

import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;
import org.w3c.dom.DocumentType;
import org.w3c.dom.Element;
import ur_rna.StructureEditor.services.drawing.export.GraphicsOp.*;
import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.Strings;
import ur_rna.Utilities.swing.FontUtil;
import ur_rna.Utilities.swing.ImageUtil;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.awt.*;
import java.awt.color.ColorSpace;
import java.awt.geom.*;
import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.Writer;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;
import java.util.function.Function;

/**
 * Saves Graphics operations (recorded by a Graphics2DRecorder) to SVG format.
 *
 * @author Richard Watson on 2/8/2017.
 */
public class SvgGraphicsExporter {
    private static final String SVG_DOCTYPE_NAME = "svg",SVG_DOCTYPE_PUBLIC_ID = "-//W3C//DTD SVG 1.1//EN",
            SVG_DOCTYPE_SYSTEM_ID = "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd",
            SVG_NAMESPACE_URI = "http://www.w3.org/2000/svg";
    private static final String XLINK_NAMESPACE = "xlink", XLINK_NAMESPACE_URI = "http://www.w3.org/1999/xlink";
    public static int ScreenDPI = 72, PrintDPI = 300, TrueScreenDPI = Toolkit.getDefaultToolkit().getScreenResolution();

    private final double DPI, PX_PER_MM;
    private static final String CHARSET_UTF8 = "UTF-8";

    private final Map<Integer, Element> clippingPathElements = new HashMap<>();
    private final Stack<GraphicsStyle> graphicsStyles = new Stack<>();
    private final Map<Element,GraphicsStyle> groupStyles = new HashMap<>(); // the styles in effect when the group was created.

    private final PageSize pageSize;
    private final Document doc;
    private Element root, group, defs;
    private GraphicsStyle currentStyle, groupStyle;

    public SvgGraphicsExporter(final PageSize pageSize, final double dpiResolution) {
        this.pageSize = pageSize;
        DPI =  dpiResolution;
        PX_PER_MM = DPI/PageSize.MM_PER_INCH;

        DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
        docFactory.setValidating(false);
        DocumentBuilder docBuilder;
        try {
            docBuilder = docFactory.newDocumentBuilder();
        } catch (ParserConfigurationException e) {
            throw new IllegalStateException("Could not create an XML builder.");
        }

        // Create a new XML document and set its DOCTYPE
        DOMImplementation domImpl = docBuilder.getDOMImplementation();
        DocumentType docType = domImpl.createDocumentType(SVG_DOCTYPE_NAME, SVG_DOCTYPE_PUBLIC_ID, SVG_DOCTYPE_SYSTEM_ID);
        doc = domImpl.createDocument(SVG_NAMESPACE_URI, "svg", docType);

        try {
            doc.setXmlStandalone(false);
        } catch (AbstractMethodError e) {
            System.err.println("Standalone XML documents not supported.");
        }

        root = doc.getDocumentElement();
        initDocRoot();
        group = root;
        currentStyle =new GraphicsStyle();
        graphicsStyles.push(currentStyle);
        groupStyles.put(group, groupStyle=currentStyle.clone());
    }

    private void initDocRoot() {
        double x = 0;
        double y = 0;
        double width = pageSize.width;
        double height = pageSize.height;

        // Add svg element
        root.setAttribute("xmlns:" + XLINK_NAMESPACE, XLINK_NAMESPACE_URI);
        root.setAttribute("version", "1.1");
        root.setAttribute("x", fmtNumber(x/PX_PER_MM) + "mm");
        root.setAttribute("y", fmtNumber(y/PX_PER_MM) + "mm");
        root.setAttribute("width", fmtNumber(width/PX_PER_MM) + "mm");
        root.setAttribute("height", fmtNumber(height/PX_PER_MM) + "mm");
        root.setAttribute("viewBox", join(" ", x, y, width, height));
    }

    public void write(List<GraphicsOp> ops, String path) throws IOException{
        try(Writer writer = newUTF8Writer(path, false)) {
            write(ops, writer);
        }
    }
    public void write(List<GraphicsOp> ops,  Writer output) throws IOException{
        for (GraphicsOp op : ops)
            applyOperation(op);
        writeTo(new StreamResult(output));
    }
    /**
     * Simplified version of {@link Files#newBufferedWriter(Path, Charset, OpenOption...)} that
     * uses {@link StandardCharsets#UTF_8} as the charset and accepts a boolean "append" instead of
     * an {@link OpenOption} argument list.
     */
    private static Writer newUTF8Writer(String path, boolean append) throws IOException {
        OpenOption[] options = append ? new OpenOption[] { StandardOpenOption.APPEND } : new OpenOption[0];
        return Files.newBufferedWriter(Paths.get(path), StandardCharsets.UTF_8, options);
    }

    /** Mapping of stroke endcap values from Java to SVG. */
    private static final Map<Integer, String> STROKE_ENDCAPS = ObjTools.toMap(
            new Integer[] { BasicStroke.CAP_BUTT, BasicStroke.CAP_ROUND, BasicStroke.CAP_SQUARE },
            new String[] { "butt", "round", "square" }
    );

    /** Mapping of line join values for path drawing from Java to SVG. */
    private static final Map<Integer, String> STROKE_LINEJOIN = ObjTools.toMap(
            new Integer[] { BasicStroke.JOIN_MITER, BasicStroke.JOIN_ROUND, BasicStroke.JOIN_BEVEL },
            new String[] { "miter", "round", "bevel" }
    );

    private static final DecimalFormat _fmtSimpleDouble = new DecimalFormat("#.######");
    private static String fmtNumber(double number) {
        return _fmtSimpleDouble.format(number);
    }
    private static String fmtObj(Object value) {
        if (value instanceof Number)
            return _fmtSimpleDouble.format(((Number)value).doubleValue());
        if (value == null)
            return "null";
        return value.toString();
    }

    public void writeTo(StreamResult out) throws IOException {
        TransformerFactory transformerFactory = TransformerFactory.newInstance();
        try {
            Transformer transformer = transformerFactory.newTransformer();
            transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            transformer.setOutputProperty(OutputKeys.STANDALONE, "no");
            transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
            transformer.setOutputProperty(OutputKeys.ENCODING, CHARSET_UTF8);
            transformer.setOutputProperty(OutputKeys.DOCTYPE_PUBLIC,
                    doc.getDoctype().getPublicId());
            transformer.setOutputProperty(OutputKeys.DOCTYPE_SYSTEM,
                    doc.getDoctype().getSystemId());
            transformer.transform(new DOMSource(doc), out);
        } catch (TransformerException e) {
            throw new IOException(e.getMessage());
        }
    }

    @Override
    public String toString() {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        try {
            writeTo(new StreamResult(out));
            return out.toString(CHARSET_UTF8);
        } catch (IOException e) {
            return "";
        }
    }

    private Element getClipElement(Shape clip) {
        // Look for existing entries
        Element path = clippingPathElements.get(clip.hashCode());
        if (path != null) {
            return path;
        }

        // Make sure <defs> exists
        if (defs == null) {
            defs = doc.createElement("defs");
            root.insertBefore(defs, root.getFirstChild());
        }

        // Store clipping path in <defs> without styling information
        path = doc.createElement("clipPath");
        path.setAttribute("id", "clip" + clip.hashCode());
        Element shape = getElement(clip);
        shape.removeAttribute("style");
        path.appendChild(shape);
        defs.appendChild(path);

        // Register path
        clippingPathElements.put(clip.hashCode(), path);

        return path;
    }

    private void newGroup() {
        Element prevGroup = group;
        group = doc.createElement("g");
        prevGroup.appendChild(group);
        initGroup(group);
    }

    private void initGroup(Element g) {
        GraphicsStyle prev = groupStyles.get((Element)g.getParentNode());
        Shape clip = currentStyle.getClip();
        if (clip != null && !Objects.equals(prev.getClip(), clip)) {
            Element clipElem = getClipElement(clip);
            String ref = "url(#" + clipElem.getAttribute("id") + ")";
            group.setAttribute("clip-path", ref);
        }

        AffineTransform tx = currentStyle.getTransform();
        if (!Objects.equals(prev.getTransform(), tx)) {
            group.setAttribute("transform", fmtTransform(tx));
        }

        Font f = currentStyle.getFont();
        if (!Objects.equals(prev.getFont(), f))
            group.setAttribute("style", getFontStyle(f));

        groupStyles.put(g, groupStyle=currentStyle.clone());
    }

    private void addToGroup(Element e) { group.appendChild(e); }

    private void exitGroup() {
        if (group == root)
            throw new IllegalStateException("Cannot exit root group.");
        group = (Element)group.getParentNode();
        groupStyle = groupStyles.get((Element)group);
    }

    public void applyOperation(GraphicsOp op) {
//        if (op instanceof GroupStartOp) {
////            GroupStartOp c = (GroupStartOp) command;
////            applyStateCommands(c.getValue());
////            if (containsGroupCommand(c.getValue())) {
////                newGroup();
////            }
//        } else

        if (op instanceof InitOp || op instanceof EndOp) {
            // do nothing for now.
        } else if (op instanceof GroupStartOp) {
            newGroup();
        } else if (op instanceof GroupEndOp) {
            exitGroup();
        } else if (op instanceof ActionOp) {
            if (group == root)
                newGroup(); // start a new group for these drawing commands
            if (op instanceof DrawImageOp) {
                DrawImageOp c = (DrawImageOp) op;
                Element e = getElement(c.getImage(), c.getX(), c.getY(), c.getWidth(), c.getHeight());
                addToGroup(e);
            } else if (op instanceof DrawShapeOp) {
                DrawShapeOp c = (DrawShapeOp) op;
                Element e = getElement(c.getShape());
                e.setAttribute("style", getDrawStyle(c.isFill()));
                addToGroup(e);
            } else if (op instanceof ClearOp) {
                ClearOp c = (ClearOp) op;
                Element e = getElement(c.getRect());
                GraphicsStyle ctx = pushState();
                ctx.setStroke(null);
                ctx.setColor(ctx.getBackground());
                ctx.setComposite(AlphaComposite.Src);
                e.setAttribute("style", getDrawStyle(true));
                popState(); // return to previous color and stroke
                addToGroup(e);
            } else if (op instanceof DrawStringOp) {
                DrawStringOp c = (DrawStringOp) op;
                Element e = getElement(c.getText(), c.getX(), c.getY());
                e.setAttribute("style", getDrawStyle(true) + getFontStyle(currentStyle.getFont()));
                addToGroup(e);
            } else
                throw new IllegalArgumentException("Unknown Graphics Operation: " + op.getClass().getSimpleName());
        } else if (op instanceof GraphicsOp.StateOp) {
            op.apply(currentStyle);
        } else
            throw new IllegalArgumentException("Unknown Graphics Operation: " + op.getClass().getSimpleName());
    }

//    private void applyStateCommands(List<GraphicsOp> commands) {
//        for (GraphicsOp command : commands) {
//            GraphicsStyle state = currentStyle;
//            if (command instanceof BackgroundOp) {
//                BackgroundOp c = (BackgroundOp) command;
//                state.setBackground(c.getValue());
//            } else if (command instanceof ClipOp) {
//                ClipOp c = (ClipOp) command;
//                state.setClip(c.getValue());
//            } else if (command instanceof ColorOp) {
//                ColorOp c = (ColorOp) command;
//                state.setColor(c.getValue());
//            } else if (command instanceof CompositeOp) {
//                CompositeOp c = (CompositeOp) command;
//                state.setComposite(c.getComposite());
//            } else if (command instanceof FontOp) {
//                FontOp c = (FontOp) command;
//                state.setFont(c.getValue());
//            } else if (command instanceof PaintOp) {
//                PaintOp c = (PaintOp) command;
//                state.setPaint(c.getValue());
//            } else if (command instanceof StrokeOp) {
//                StrokeOp c = (StrokeOp) command;
//                state.setStroke(c.getValue());
//            } else if (command instanceof TransformOp) {
//                TransformOp c = (TransformOp) command;
//                state.setTransform(c.getValue());
//            } else if (command instanceof AffineTransformOp) {
//                AffineTransformOp c = (AffineTransformOp) command;
//                AffineTransform stateTransform = state.getTransform();
//                AffineTransform transformToBeApplied = c.getValue();
//                stateTransform.concatenate(transformToBeApplied);
//                state.setTransform(stateTransform);
//            } else if (command instanceof SetHintOp) {
//                SetHintOp c = (SetHintOp) command;
//                state.getHints().put(c.getKey(), c.getValue());
//            } else if (command instanceof CreateOp) {
//                try {
//                    graphicsStyles.push((GraphicsStyle) currentStyle.clone());
//                } catch (CloneNotSupportedException e) {
//                    e.printStackTrace();
//                }
//            } else if (command instanceof DisposeOp) {
//                graphicsStyles.pop();
//            }
//        }
//    }

//    private boolean containsGroupCommand(List<Command<?>> commands) {
//        for (Command<?> command : commands) {
//            if ((command instanceof SetClipCommand) ||
//                    (command instanceof SetTransformCommand) ||
//                    (command instanceof AffineTransformCommand)) {
//                return true;
//            }
//        }
//        return false;
//    }

    /** Pushes a new GraphicsStyle (cloned from the current one) to the top of the stack and returns it. */
    private GraphicsStyle pushState() {
        return graphicsStyles.push(currentStyle = currentStyle.clone());
    }
    /** Pops the GraphicsStyle from the top of the stack and returns it. */
    private GraphicsStyle popState() {
        GraphicsStyle prev = graphicsStyles.pop();
        currentStyle = graphicsStyles.peek();
        return prev;
    }
    private String getDrawStyle(boolean filled) {
        StringBuilder style = new StringBuilder();

        Color color = currentStyle.getColor();
        String colorOutput = getOutput(color);
        double opacity = color.getAlpha()/255.0;

        if (filled) {
            appendStyle(style, "fill", colorOutput);
            if (color.getAlpha() < 255) {
                appendStyle(style, "fill-opacity", opacity);
            }
        } else {
            appendStyle(style, "fill", "none");
        }

        if (!filled) {
            appendStyle(style, "stroke", colorOutput);
            if (color.getAlpha() < 255) {
                appendStyle(style, "stroke-opacity", opacity);
            }
            Stroke stroke = currentStyle.getStroke();
            if (stroke instanceof BasicStroke) {
                BasicStroke bs = (BasicStroke) stroke;
                if (bs.getLineWidth() != 1f) {
                    appendStyle(style, "stroke-width", bs.getLineWidth());
                }
                if (bs.getMiterLimit() != 4f) {
                    appendStyle(style, "stroke-miterlimit", bs.getMiterLimit());
                }
                if (bs.getEndCap() != BasicStroke.CAP_BUTT) {
                    appendStyle(style, "stroke-linecap", STROKE_ENDCAPS.get(bs.getEndCap()));
                }
                if (bs.getLineJoin() != BasicStroke.JOIN_MITER) {
                    appendStyle(style, "stroke-linejoin", STROKE_LINEJOIN.get(bs.getLineJoin()));
                }
                if (bs.getDashArray() != null) {
                    appendStyle(style, "stroke-dasharray", join(",", bs.getDashArray()));
                    if (bs.getDashPhase() != 0f) {
                        appendStyle(style, "stroke-dashoffset", bs.getDashPhase());
                    }
                }
            }
        } else {
            appendStyle(style, "stroke", "none");
        }

        return style.toString();
    }

    private String getFontStyle(Font font) {
        if (Objects.equals(groupStyle.getFont(), font)) return "";

        StringBuilder out = new StringBuilder();
        //if (!GraphicsStyle.DEFAULT_FONT.getFamily().equals(font.getFamily())) {
        String physicalFamily = FontUtil.getPhysicalFont(font).getFamily();
        out.append("font-family:\"").append(physicalFamily).append("\";");
        //}
        //if (font.getSize2D() != GraphicsStyle.DEFAULT_FONT.getSize2D()) {
            // Java Font sizes are in units of "points", assuming a fixed-screen resolution of 72 PPI
            // So the conversion to pixels is simply  SIZE (points) / 72 (points/inch) * 72 (pixels/inch) = SIZE (pixels)
            // Note that most screen resolutions are NOT 72 PPI, but rather 96 PPI, so it is more accurate to
            // say that Java's fonts are in pixel units, NOT true points. Modern web browsers, in contrast,
            // often do use the true screen resolution, so if a font is specified as 12pt (at 96 PPI) it would have
            // a pixel size of 12pt / 72 (pt/in) * 96 (px/in) = 16px.
            out.append("font-size:").append(fmtNumber(font.getSize2D())).append("px;");
        //}
        if ((font.getStyle() & Font.ITALIC) != 0) {
            out.append("font-style:italic;");
        }
        if ((font.getStyle() & Font.BOLD) != 0) {
            out.append("font-weight:bold;");
        }
        return out.toString();
    }

    private static void appendStyle(StringBuilder style, String attribute, Object value) {
        style.append(attribute).append(":")
                .append(fmtObj(value)).append(";");
    }

    private static String fmtTransform(AffineTransform tx) {
        StringBuilder out = new StringBuilder();

        int type = tx.getType();

        if (type == AffineTransform.TYPE_IDENTITY) return "";


        if (type == AffineTransform.TYPE_GENERAL_TRANSFORM ||
                0!=(type & AffineTransform.TYPE_MASK_ROTATION)) {
            double[] matrix = new double[6];
            tx.getMatrix(matrix);
            out.append("matrix(").append(join(" ", matrix)).append(")");
            return out.toString();
        }

        // TODO: We can actually calculate the rotation from the AffineTransform, and specify it explicitly ...
        // potentially along with a translation and scale.
//        private static double extractAngle(AffineTransform at) {
//            Point2D p0 = new Point2D.Double();
//            Point2D p1 = new Point2D.Double(1,0);
//            Point2D pp0 = at.transform(p0, new Point2D.Double());
//            Point2D pp1 = at.transform(p1, new Point2D.Double());
//            double dx = pp1.getX() - pp0.getX();
//            double dy = pp1.getY() - pp0.getY();
//            double angle = Math.atan2(dy, dx);
//            return angle;
//        }


        if ((type & AffineTransform.TYPE_TRANSLATION)!=0) {
                out.append("translate(")
                        .append(fmtNumber(tx.getTranslateX())).append(" ")
                        .append(fmtNumber(tx.getTranslateY())).append(") ");
        }

        if ((type & AffineTransform.TYPE_MASK_SCALE)!=0) {
            out.append("scale(")
                    .append(fmtNumber(tx.getScaleX())).append(" ")
                    .append(fmtNumber(tx.getScaleY())).append(") ");
        }

        return out.toString();
    }

    private static String join(String delimiter, double... values) {
        return Strings.join(delimiter, ObjTools.toList(values), SvgGraphicsExporter::fmtNumber);
    }
    private static String join(String delimiter, float... values) {
        return Strings.join(delimiter, ObjTools.toList(values), (Function<Float,String>)SvgGraphicsExporter::fmtNumber);
    }

    private static String getOutput(Color color) {
        if (color.getColorSpace().getType() == ColorSpace.TYPE_CMYK) {
            float[] cmyk = color.getComponents(null);
            return String.format((Locale) null,
                    "rgb(%d,%d,%d) icc-color(Generic-CMYK-profile,%f,%f,%f,%f)",
                    color.getRed(), color.getGreen(), color.getBlue(),
                    cmyk[0], cmyk[1], cmyk[2], cmyk[3]);
        } else {
            return String.format((Locale) null, "rgb(%d,%d,%d)",
                    color.getRed(), color.getGreen(), color.getBlue());
        }
    }

    private static String getOutput(Shape shape) {
        StringBuilder out = new StringBuilder();
        PathIterator segments = shape.getPathIterator(null);
        double[] coords = new double[6];
        for (int i = 0; !segments.isDone(); i++, segments.next()) {
            if (i > 0) {
                out.append(" ");
            }
            int segmentType = segments.currentSegment(coords);
            switch (segmentType) {
                case PathIterator.SEG_MOVETO:
                    out.append("M").append(coords[0]).append(",").append(coords[1]);
                    break;
                case PathIterator.SEG_LINETO:
                    out.append("L").append(coords[0]).append(",").append(coords[1]);
                    break;
                case PathIterator.SEG_CUBICTO:
                    out.append("C")
                            .append(coords[0]).append(",").append(coords[1]).append(" ")
                            .append(coords[2]).append(",").append(coords[3]).append(" ")
                            .append(coords[4]).append(",").append(coords[5]);
                    break;
                case PathIterator.SEG_QUADTO:
                    out.append("Q")
                            .append(coords[0]).append(",").append(coords[1]).append(" ")
                            .append(coords[2]).append(",").append(coords[3]);
                    break;
                case PathIterator.SEG_CLOSE:
                    out.append("Z");
                    break;
                default:
                    throw new IllegalStateException("Unknown path operation.");
            }
        }
        return out.toString();
    }

    private static String getImageOutput(Image image, boolean lossyAllowed) {
        BufferedImage bufferedImage = ImageUtil.toBufferedImage(image);

        String encoded = encodeImage(bufferedImage, "png");
        if (!ImageUtil.usesAlpha(bufferedImage) && lossyAllowed) {
            String encodedLossy = encodeImage(bufferedImage, "jpeg");
            if (encodedLossy.length() > 0 && encodedLossy.length() < encoded.length()) {
                encoded = encodedLossy;
            }
        }

        return encoded;
    }

    private static String encodeImage(BufferedImage bufferedImage, String format) {
        throw new RuntimeException("Not implemented.");
//        ByteArrayOutputStream byteStream = new ByteArrayOutputStream();
//        Base64EncodeStream encodeStream = new Base64EncodeStream(byteStream);
//        try {
//            ImageIO.write(bufferedImage, format, encodeStream);
//            encodeStream.close();
//            String encoded = byteStream.toString("ISO-8859-1");
//            return String.format("data:image/%s;base64,%s", format, encoded);
//        } catch (IOException e) {
//            return "";
//        }
    }

    private Element getElement(Shape shape) {
        Element elem;
        if (shape instanceof Line2D) {
            Line2D s = (Line2D) shape;
            elem = doc.createElement("line");
            elem.setAttribute("x1", fmtNumber(s.getX1()));
            elem.setAttribute("y1", fmtNumber(s.getY1()));
            elem.setAttribute("x2", fmtNumber(s.getX2()));
            elem.setAttribute("y2", fmtNumber(s.getY2()));
        } else if (shape instanceof Rectangle2D) {
            Rectangle2D s = (Rectangle2D) shape;
            elem = doc.createElement("rect");
            elem.setAttribute("x", fmtNumber(s.getX()));
            elem.setAttribute("y", fmtNumber(s.getY()));
            elem.setAttribute("width", fmtNumber(s.getWidth()));
            elem.setAttribute("height", fmtNumber(s.getHeight()));
        } else if (shape instanceof RoundRectangle2D) {
            RoundRectangle2D s = (RoundRectangle2D) shape;
            elem = doc.createElement("rect");
            elem.setAttribute("x", fmtNumber(s.getX()));
            elem.setAttribute("y", fmtNumber(s.getY()));
            elem.setAttribute("width", fmtNumber(s.getWidth()));
            elem.setAttribute("height", fmtNumber(s.getHeight()));
            elem.setAttribute("rx", fmtNumber(s.getArcWidth()/2.0));
            elem.setAttribute("ry", fmtNumber(s.getArcHeight()/2.0));
        } else if (shape instanceof Ellipse2D) {
            Ellipse2D s = (Ellipse2D) shape;
            elem = doc.createElement("ellipse");
            elem.setAttribute("cx", fmtNumber(s.getCenterX()));
            elem.setAttribute("cy", fmtNumber(s.getCenterY()));
            elem.setAttribute("rx", fmtNumber(s.getWidth()/2.0));
            elem.setAttribute("ry", fmtNumber(s.getHeight()/2.0));
        } else {
            elem = doc.createElement("path");
            elem.setAttribute("d", getOutput(shape));
        }
        return elem;
    }

    private Element getElement(String text, double x, double y) {
        Element elem = doc.createElement("text");
        elem.appendChild(doc.createTextNode(text));
        elem.setAttribute("x", fmtNumber(x));
        elem.setAttribute("y", fmtNumber(y));
        return elem;
    }

    private Element getElement(Image image, double x, double y, double width, double height) {
        Element elem = doc.createElement("image");
        elem.setAttribute("x", fmtNumber(x));
        elem.setAttribute("y", fmtNumber(y));
        elem.setAttribute("width", fmtNumber(width));
        elem.setAttribute("height", fmtNumber(height));
        elem.setAttribute("preserveAspectRatio", "none");
//        boolean lossyAllowed = currentStyle.getHints().get(VectorHints.KEY_EXPORT) ==
//                VectorHints.VALUE_EXPORT_SIZE;
        elem.setAttribute("xlink:href", getImageOutput(image, true /*lossyAllowed*/));
        return elem;
    }
}
