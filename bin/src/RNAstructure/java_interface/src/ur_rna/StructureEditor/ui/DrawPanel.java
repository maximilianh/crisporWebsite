package ur_rna.StructureEditor.ui;

import ur_rna.StructureEditor.models.HistoryUpdateEvent;
import ur_rna.StructureEditor.models.IScreenObject;
import ur_rna.StructureEditor.services.RnaDrawController;
import ur_rna.StructureEditor.services.SceneRenderer;
import ur_rna.StructureEditor.services.SceneUpdateCategory;
import ur_rna.StructureEditor.services.drawing.ICanvas;
import ur_rna.StructureEditor.services.drawing.View2D;
import ur_rna.Utilities.EventSource;
import ur_rna.Utilities.annotation.Nullable;

import javax.swing.*;
import java.awt.*;
import java.awt.event.InputEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Dimension2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.Iterator;

/**
 * A double-buffered canvas on which the RNA scene can be drawn, implemented on top of a JPanel
 */
public class DrawPanel extends JPanel implements ICanvas {
    private float scale = 1f;
    private Point2D.Float offset = new Point2D.Float();
    private AffineTransform tr = new AffineTransform();

    private final int MARGIN = 5;
    private View2D view = new View2D(), drawnView = new View2D(); //, scrolledView = new View2D(); // view is the "target" view that will become drawnView once the paint operation is invoked.

    private SceneRenderer renderer;
    private static final int ZOOM_MASK = InputEvent.CTRL_DOWN_MASK;

    public final EventSource.NoArg viewChanged = new EventSource.NoArg();

    public DrawPanel() {
        setDoubleBuffered(true);
        addMouseWheelListener(e -> {
            int mod = e.getModifiersEx();
            if ((mod & ZOOM_MASK)==ZOOM_MASK) {
                float newScale = (float) (scale - 0.1 * scale * e.getPreciseWheelRotation());
                //System.out.println("Scale when wheeled: " + scale + " -> " + newScale);
                // get the point where the mouse currently is (on the full content panel)
                setZoomScroll(e.getPoint());
//
//                // Adjust the point to account for the change in scale
//                loc.x = (int)(loc.x*(newScale/drawnScale - 1)-this.getX());
//                loc.y = (int)(loc.y*(newScale/drawnScale - 1)-this.getY());
//
//                float tmp = drawnScale;
//
//                // add an event that will fire when the viewport actually changes size (after revalidate events etc)
//                SwingUtilities.invokeLater(()->verifyScrollToPoint(tmp, loc) );
                setScale(newScale);

                e.consume();
            } else {
                delegateToParent(e);
            }
        });


//        AffineTransform tt = new AffineTransform();
//        tt.setToTranslation(20, 20);
//        tt.scale(10, 10);
//        //tt.scale(100,100);
//        Point2D[] pts = new Point2D[]{
//                new Point2D.Float(-20, -20),
//                new Point2D.Float(0, 0),
//                new Point2D.Float(10, 10),
//                new Point2D.Float(10, 0)
//        };
//        for (int i = 0; i < pts.length; i++) {
//            Point2D pt = tt.transform(pts[i], null);
//            System.out.println("pt" + i + ": [" + pt.getX() + ", " + pt.getY() + "]");
//        }
    }

    private void setZoomScroll(Point screenLocation) {
        //log("mouse point: %s", screenLocation);
        // Convert the location from panel coords to model space using the most recent drawn coordinates.
        Point2D.Float model = new Point2D.Float();
        drawnView.toModel(screenLocation, model);
        //log("model point: %s", model);

        // Translate to view coordinates
        screenLocation.x += getX();
        screenLocation.y += getY();

        //log("view point: %s", screenLocation);

        scrollTargetModel = model;
        scrollTargetView = screenLocation;
    }

    private static void log(String format, Object... args) { log(false, format, args); }
    private static void loge(String format, Object... args) { log(true, format, args); }
    private static void log(boolean useError, String format, Object... args) {
        for (int i = 0; i < args.length; i++) {
            if (args[i] instanceof Rectangle)
                args[i] = fmt((Rectangle)args[i]);
            else if (args[i] instanceof Point)
                args[i] = fmt((Point)args[i]);
            else if (args[i] instanceof Dimension)
                args[i] = fmt((Dimension)args[i]);
        }
        (useError?System.err:System.out).println(args.length==0?format:String.format(format, args));
    }

    Point2D scrollTargetModel; Point scrollTargetView; // Dimension targetSize;
//    protected void processComponentEvent(final ComponentEvent e) {
//        super.processComponentEvent(e);
//        if (e.getID() == ComponentEvent.COMPONENT_RESIZED) {
//            log("Resized to %s", getSize());
//            if (getParent() != null)
//            log("HBar: %s, %s", ((JScrollPane) getParent().getParent()).getHorizontalScrollBar().getMaximum(),
//                    ((JViewport)getParent()).getViewSize().width - ((JViewport)getParent()).getViewRect().width);
//
//            if (getSize().equals(targetSize)) {
//                scrolledView.setView(drawnView);
//                handleZoomScroll();
//                scrollTargetModel = null;
//            } else {
//                loge("not sized yet");
//            }
//
////            Point2D p = scrollTargets.poll();
////            if (p != null)
////                scrollToScenePoint(p);
//        }
//    }

//    protected void processComponentEvent(final ComponentEvent e) {
//        super.processComponentEvent(e);
//        if (e.getID() == ComponentEvent.COMPONENT_MOVED) {
//            log("Moved to %s", getLocation());
//        }
//    }

    private Point getZoomScroll(Dimension contentBounds, Dimension viewRect) {
        if (scrollTargetModel==null) return null;
        Point target = new Point();
        drawnView.toScreen(scrollTargetModel, target); // use drawnView here instead of scrollView because we want to anticipate the next view.

        //log("target point: %s", target);
        target.x -= scrollTargetView.x;
        target.y -= scrollTargetView.y;
        //log("target point on view: %s", target);

        if (viewRect == null) {
            if (getParent() instanceof JViewport)
                viewRect = ((JViewport)getParent()).getExtentSize();
            else
                viewRect = getVisibleRect().getSize();
        }
        target.x = Math.max(0, Math.min(target.x, contentBounds.width - viewRect.width));
        target.y = Math.max(0, Math.min(target.y, contentBounds.height - viewRect.height));

        return target;
    }

//    private void verifyScrollToPoint(final float requiredScale, final Point loc) {
//        if (requiredScale == drawnScale) {
//            scrollToScenePoint(loc);
//            System.out.printf("verifyScrollToPoint r:%s, d:%s, cur:%s loc:%s\n", requiredScale, drawnScale, this.getScale(), fmt(loc));
//        } else
//            System.err.printf("verifyScrollToPoint r:%s, d:%s, cur:%s loc:%s\n", requiredScale, drawnScale, this.getScale(), fmt(loc));
//    }
    public void scrollToPoint(final Point p) {
        if (p==null|| !(getParent() instanceof JViewport)) return;
        JViewport port = (JViewport)getParent();

        Rectangle viewRect = port.getViewRect(); // actual size
        //System.out.printf("P: %s,  frame:%s, content:%s, max:%s, \n", fmt(p),  fmt(frame), fmt(content), fmt(new Dimension(content.width - frame.width, content.height - frame.height)));

        //System.out.printf("P: %s\n", fmt(p));
        //p.x -= target.x;
        //p.y -= target.y;
        //System.out.printf("P.new: %s\n", fmt(p));
        //p.x = Math.max(0, Math.min(p.x, content.width - frame.width));
        //p.y = Math.max(0, Math.min(p.y, content.height - frame.height));
        //System.out.printf("P.ajd: %s\n", fmt(p));

        Dimension contentBounds = getSize();
        p.x = Math.max(0, Math.min(p.x, contentBounds.width - viewRect.width));
        p.y = Math.max(0, Math.min(p.y, contentBounds.height - viewRect.height));

//        p.x = Math.max(0, Math.min(p.x, maxx));
//        p.y = Math.max(0, Math.min(p.y, maxy));

        if (getX()!=-p.x||getY()!=-p.y)
            setLocation(-p.x, -p.y);
        port.setViewPosition(p);

//        if (!p.equals(port.getViewPosition()))
//            port.setViewPosition(p);
    }

//    @Override
//    public void reshape(final int x, final int y, final int w, final int h) {
//        log("reshape %s, %s, %s, %s", x, y, w, h);
//        super.reshape(x, y, w, h);
//    }
    /** bubble the mouse wheel event up to the parent JScrollPane */
    private void delegateToParent(MouseWheelEvent e) {
        // even with scroll bar set to never the event doesn't reach the parent scroll frame
        JScrollPane ancestor = (JScrollPane) SwingUtilities.getAncestorOfClass(
                JScrollPane.class, this);
        if (ancestor != null) {
            MouseWheelEvent converted = null;
            for (MouseWheelListener listener : ancestor
                    .getMouseWheelListeners()) {
                listener.mouseWheelMoved(converted != null ? converted
                        : (converted = (MouseWheelEvent) SwingUtilities
                        .convertMouseEvent(this, e, ancestor)));
            }
        }
    }

    private void updateTransform() {
        tr.setToScale(scale, scale);
        tr.translate(offset.x, offset.y);
        updateView();
        viewChanged.invoke();
        repaint();
    }

    public void updateView() {
        view.setView(tr, getVisibleRect());
        view.screenMargin.setSize(MARGIN, MARGIN);
    }

    public void setRenderer(final SceneRenderer renderer) {
        this.renderer = renderer;
        if (renderer instanceof RnaDrawController) {
            ((RnaDrawController) renderer).historyEvent.add(this::controllerHistoryUpdate);
        }
    }
    private void controllerHistoryUpdate(final HistoryUpdateEvent event) {
        if (event.storing) {
            event.state.put("view-scale", scale);
            event.state.put("view-offset", new Point2D.Float(offset.x, offset.y));
        } else {
            if (event.state.update.category == SceneUpdateCategory.Layout)
                setView(event.state.get("view-scale", scale), event.state.get("view-offset", offset));
        }
    }
    public SceneRenderer getRenderer() {
        return renderer;
    }

    @Override
    public void repaintRequired() {
        repaint();
    }

    @Override
    public Component getComponent() {
        return this;
    }
    public View2D getView() {  return view;  }
//    @Override
//    public void handleEvent(final SceneUpdateEvent eventInfo) {
//        repaint();
//    }
    @Override
    protected void paintComponent(final Graphics gg) {
        //log("Painting");
        if (renderer != null) {
            Graphics2D g = (Graphics2D) gg;
            double scale =  view.trToScreen.getScaleX() / drawnView.trToScreen.getScaleX();
            drawnView.setView(view);
            if (!Double.isNaN(scale) && (scale > 1.000001 || scale < .999999)) {
                if (prevBounds==null) prevBounds = getBounds();
                // log("Current bounds: %s", getBounds());

                // Estimate new bounds based on previous bounds and the change in scale.
                prevBounds.width = (int)((prevBounds.width-MARGIN) * scale + MARGIN);
                prevBounds.height = (int)((prevBounds.height-MARGIN) * scale + MARGIN);

                if (getParent() instanceof JViewport) {
                    Dimension port = ((JViewport) getParent()).getExtentSize();
                    Point p = getZoomScroll(prevBounds.getSize(), port);
                    if (p != null) {
                        Point loc = getLocation();
                        g.translate(-(loc.x + p.x), -(loc.y + p.y));
                        prevBounds.setLocation(-p.x, -p.y);
                    }
                    if (prevBounds.width < port.width) prevBounds.width = port.width;
                    if (prevBounds.height < port.height) prevBounds.height = port.height;
                }
                this.setBounds(prevBounds);
                // log("New bounds: %s == %s at %s", getBounds(), prevBounds, getLocation());
                //drawnView.setView(tr,  getVisibleRect());
            }
            // Try to guess the best size. if we are wrong, it will redraw.
            updateContentSize(getBounds(renderer.render(g, drawnView, true)));
        }
    }
    private Rectangle prevBounds;

    private void updateContentSize(Rectangle bounds) {
        boolean allowOffsetChange = true;
        if (bounds == null) return;

        if (renderer instanceof RnaDrawController) {
            RnaDrawController rdc = (RnaDrawController)renderer;
            allowOffsetChange = !rdc.isUserActionInProgress();
        }

        bounds.grow(MARGIN, MARGIN);

        float dx = bounds.x < 0 ? -bounds.x : 0;
        float dy = bounds.y < 0 ? -bounds.y : 0;

        if (dx != 0 || dy != 0) {
            if (allowOffsetChange) {
                // System.out.println(fmt("Scale: %s, Current: (%s, %s) New: (%s, %s) Target: (%s, %s)",  scale, offset.x, offset.y, offset.x+ dx/scale, offset.y+ dy/scale, offset.x + dx, offset.y + dy));
                setOffset(offset.x + dx / scale, offset.y + dy / scale);
                return;
            }
            // for view calculation
            bounds.translate((int)Math.ceil(dx), (int)Math.ceil(dy));
        }

        // Calculate new view rectangle.
        Dimension size = new Dimension(bounds.x + bounds.width, bounds.y + bounds.height);
        if (!this.getPreferredSize().equals(size)) {
            //log("Scale when revalidated: %s", drawnScale);
//            log("HBar: %s, %s", ((JScrollPane)((JViewport)this.getParent()).getParent()).getHorizontalScrollBar().getMaximum(),
//                    ((JViewport)this.getParent()).getViewSize().width - ((JViewport)this.getParent()).getViewRect().width
//                    );
            // targetSize = size;
            this.setPreferredSize(size);

            if (prevBounds==null)
                prevBounds = getBounds();
            else
                prevBounds.setSize(size);

            Dimension port = ((JViewport)this.getParent()).getExtentSize();
            if (size.width < port.width) size.width = port.width;
            if (size.height < port.height) size.height = port.height;
            this.setSize(size);
            scrollToPoint(getZoomScroll(size, port));
            scrollTargetModel = null;
            revalidate();

//            log("HBar: %s, %s", ((JScrollPane)((JViewport)this.getParent()).getParent()).getHorizontalScrollBar().getMaximum(),
//                    ((JViewport)this.getParent()).getViewSize().width - ((JViewport)this.getParent()).getViewRect().width
//            );
            // this.revalidate(); // forces owner (ScrollPane) to re-evaluate size and scrollbars.
            // log("Updated bounds: %s", getBounds());
        }
    }
//    public void showViewInfo() {
//        JViewport view = ((JViewport)this.getParent());
//        System.out.printf("Frame: %s, content: %s, pref-size: %s size: %s\n", fmt(view.getViewRect()), fmt(view.getViewSize()), fmt(this.getPreferredSize()), fmt(this.getSize()));
//    }
    private static String fmt(Rectangle2D rc) {
        if (rc==null) return "NULL";
        return String.format("[%5.0fx,%5.0fy,%5.0fw,%5.0fh]", rc.getX(), rc.getY(), rc.getWidth(), rc.getHeight());
    }
    private static String fmt(Point2D p) {
        if (p==null) return "NULL";
        return String.format("[%5.0fx,%5.0fy]", p.getX(), p.getY());
    }
    private static String fmt(Dimension2D d) {
        if (d==null) return "NULL";
        return String.format("[%5.0fw,%5.0fh]", d.getWidth(), d.getHeight());
    }

    public void setView(final float scale, final Point2D offset) {
        setView(scale, (float)offset.getX(), (float)offset.getY());
    }
    public void setView(final float scale, final float offsetX, final float offsetY) {
        offset.x = offsetX;
        offset.y = offsetY;
        this.scale = scale;
        updateTransform();
    }
    public Point2D getOffset() {
        return new Point2D.Float(offset.x, offset.y);
    }
    public void setOffset(final float dx, final float dy) { setView(this.scale, dx, dy); }
    public void setOffset(final Point2D offset) { setView(this.scale, offset); }
    public float getScale() {
        return scale;
    }
    public void setScale(final float newScale) { setView(newScale, offset.x, offset.y); }

    final View2D IDENTITY_VIEW = new View2D();
    public void autoScale(Dimension preferredSize) {
        Rectangle rc = calcPureRender();
        if (rc == null) return;
        float w = (float)(preferredSize.width - 2 * MARGIN) / (rc.width + 2 * MARGIN);
        float h = (float)(preferredSize.height - 2 * MARGIN) / (rc.height + 2 * MARGIN);
        float newScale = Math.min(w, h);
        setView(newScale, -rc.x + MARGIN / newScale, -rc.y + MARGIN / newScale);
    }

    @Nullable
    private Rectangle getBounds(Iterable<IScreenObject> objs) {
        Iterator<IScreenObject> i = objs.iterator();
        if (!i.hasNext()) return null;
        Rectangle bounds = i.next().getBounds();
        while(i.hasNext())
            bounds = bounds.union(i.next().getBounds());
        return bounds;
    }

    public void zoom(final int steps) {
        float newScale = (float) (scale + 0.1 * scale * steps);
        // attempt to preserve center of view
        if (getParent() instanceof JViewport) {
            JViewport port = (JViewport)getParent();
            Dimension size = port.getExtentSize();
            Point center = new Point(size.width/2- getX(), size.height/2- getY());
            setZoomScroll(center);
        }
        setScale(newScale);
    }

    /**
     * Renders the scene with no applied transform and no output to the device (because it is clipped).
     */
    private Rectangle calcPureRender() {
        if (renderer == null) return null;
        Graphics2D g = (Graphics2D) this.getGraphics();
        if (g == null) return null;
        g.setClip(0,0,0,0);
        Rectangle rc = renderer.calcBounds(g, IDENTITY_VIEW, true);
        g.dispose();
        return rc;
    }
    public void resetOffset() {
        // The offset automatically grows when the coordinates of the scene are negative, but the excess is not trimmed unless this method is called.
        Rectangle rc = calcPureRender();
        if (rc == null) return;
        rc.grow(MARGIN, MARGIN); // 5 px margin on all sides
        setOffset(-rc.x, -rc.y);
    }
}
