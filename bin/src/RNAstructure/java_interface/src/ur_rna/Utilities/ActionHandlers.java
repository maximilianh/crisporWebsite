package ur_rna.Utilities;

import java.awt.event.ActionListener;
import java.util.concurrent.Callable;
import java.util.function.Consumer;

/**
 * @author Richard M. Watson
 */
public abstract class ActionHandlers {
    public static ActionHandler fromRunnable(final Runnable r) {
        return new ActionHandler() {
            @Override
            public void handle(Object source, EventArgs e) {
                r.run();
            }
        };
    }
    public static ActionHandler fromCallable(final Callable c) {
        return new ActionHandler() {
            @Override
            public void handle(Object source, EventArgs e) {
                try {
                    c.call();
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
            }
        };
    }
    public static ActionHandler fromConsumer(final Consumer<EventArgs> consumer) {
        return new ActionHandler() {
            @Override
            public void handle(Object source, EventArgs e) {
                consumer.accept(e);
            }
        };
    }

    public static <T> ActionHandler<EventArgs.ObjEventArgs<T>> fromConsumerT(final Consumer<T> consumer) {
        return new ActionHandler<EventArgs.ObjEventArgs<T>>() {
            public void handle(Object source, EventArgs.ObjEventArgs<T> e) {
                consumer.accept(e.value);
            }
        };
    }

    public static ActionHandler<EventArgs.ActionEventArgs> fromEventListner(final ActionListener ev) {
        return new ActionHandler<EventArgs.ActionEventArgs>() {
            @Override
            public void handle(Object source, EventArgs.ActionEventArgs e) {
                ev.actionPerformed(e.event);
            }
        };
    }
}
