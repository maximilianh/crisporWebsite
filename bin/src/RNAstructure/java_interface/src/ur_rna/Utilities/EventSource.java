package ur_rna.Utilities;

import ur_rna.Utilities.annotation.NotNull;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Consumer;

/**
 *
 */
public abstract class EventSource<TSubscriber> {
    public static class Manager {
        private static class Subscriber {
            public final Object handler, event;
            public <T> Subscriber(final EventSource<T> event, final T handler) {
                this.handler = handler;
                this.event = event;
            }
        }
        private ArrayList<Subscriber> subscribers = null;

        public <T> boolean add(EventSource<T> source, T handler) {
            synchronized (this) {
                if (subscribers == null) subscribers = new ArrayList<>();
                return subscribers.add(new Subscriber(source, handler));
            }
        }
        public <T> boolean rem(EventSource<T> source, T handler) {
            synchronized (this) {
                if (subscribers == null) return false;
                int length = subscribers.size();
                for (int i = 0; i < length; i++) {
                    if (subscribers.get(i).event == source && subscribers.get(i).handler == handler) {
                        subscribers.remove(i);
                        return true;
                    }
                }
            }
            return false;
        }
        @SuppressWarnings("unchecked")
        public <T> List<T> get(EventSource<T> source) {
            synchronized (this) {
                List<T> list = null;
                if (subscribers != null) {
                    for (Subscriber subscriber : subscribers) {
                        if (subscriber.event == source) {
                            if (list == null) list = new ArrayList<T>();
                            list.add((T) subscriber.handler);
                        }
                    }
                }
                return list == null ? Collections.emptyList() : list;
            }
        }
    }

    public EventSource() { this (new Manager()); }
    public EventSource(Manager m) { manager = m; }

    private final Manager manager;
    public Manager getManager() { return manager; }

    protected void invoke(Consumer<TSubscriber> call) {
        for (TSubscriber o : manager.get(this))
            call.accept(o);
    }
    public boolean add(@NotNull TSubscriber handler) {
        return manager.add(this, handler);
    }

    public boolean remove(TSubscriber handler) {
        return manager.rem(this, handler);
    }

    public static class NoArg extends EventSource<Runnable> {
        public void invoke() { super.invoke(Runnable::run); }
    }

    public static class OneArg<TArg> extends EventSource<Consumer<TArg>> {
        public void invoke(TArg arg) { super.invoke(c -> c.accept(arg)); }
    }
    public static class TwoArgs<TArg1, TArg2> extends EventSource<BiConsumer<TArg1, TArg2>> {
        public void invoke(TArg1 arg1, TArg2 arg2) { super.invoke(o -> o.accept(arg1, arg2)); }

    }
}


//public abstract class EventSource<TSubscriber> {
//    private ArrayList<TSubscriber> subscribers = null;
//
//    protected void invoke(Consumer<TSubscriber> call) {
//        if (subscribers != null)
//            for (TSubscriber o : subscribers)
//                call.accept(o);
//    }
//    public void add(@NotNull TSubscriber handler) {
//        if (subscribers == null)
//            subscribers = new ArrayList<>();
//        subscribers.add(handler);
//    }
//    public boolean remove(TSubscriber handler) {
//        return subscribers != null && subscribers.remove(handler);
//    }
//
//    public static class OneArg<TArg> extends EventSource<Consumer<TArg>> {
//        public void invoke(TArg arg) { super.invoke(c -> c.accept(arg)); }
//    }
//    public static class TwoArgs<TArg1, TArg2> extends EventSource<BiConsumer<TArg1, TArg2>> {
//        public void invoke(TArg1 arg1, TArg2 arg2) { super.invoke(o -> o.accept(arg1, arg2)); }
//    }
//}
