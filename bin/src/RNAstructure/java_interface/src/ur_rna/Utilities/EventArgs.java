package ur_rna.Utilities;

import java.awt.event.ActionEvent;
import java.util.EventObject;
import java.util.HashMap;

/**
 * @author Richard M. Watson
 */
public class EventArgs {
    public static final EventArgs Default = new EventArgs();

    public static class EventObjectArgs extends EventArgs {
        public final EventObject event;
        public EventObjectArgs(final EventObject event) {
            this.event = event;
        }
    }

    public static class ActionEventArgs extends EventArgs {
        public final ActionEvent event;
        public ActionEventArgs(final ActionEvent event) {
            this.event = event;
        }
    }

    public static class StringEventArgs extends EventArgs {
        public final String value;
        public StringEventArgs(final String value) {
            this.value = value;
        }
    }

    public static class IntEventArgs extends EventArgs {
        public final int value;
        public IntEventArgs(final int value) {
            this.value = value;
        }
    }

    public static class BoolEventArgs extends EventArgs {
        public final boolean value;
        public BoolEventArgs(final boolean value) {
            this.value = value;
        }
    }

    public static class ObjEventArgs<T> extends EventArgs {
        public final T value;
        public ObjEventArgs(final T value) {
            this.value = value;
        }
    }

    public static class InfoBagEventArgs extends EventArgs {
        public final HashMap<String, Object> info = new HashMap<>();
        public InfoBagEventArgs() {}
        public InfoBagEventArgs(String name, Object arg) {
            info.put(name, arg);
        }
        public InfoBagEventArgs(String name1, Object arg1, String name2, Object arg2) {
            info.put(name1, arg1);
            info.put(name2, arg2);
        }
        public InfoBagEventArgs(String name1, Object arg1, String name2, Object arg2, String name3, Object arg3) {
            info.put(name1, arg1);
            info.put(name2, arg2);
            info.put(name3, arg3);
        }
        public InfoBagEventArgs(Object... nameArgPairs) {
            if ((nameArgPairs.length & 1)==1)
                throw new IllegalArgumentException("The list of name-value argument pairs must have an even number of items.");
            for (int i = 0; i < nameArgPairs.length; i+=2)
                info.put(ObjTools.toStr(nameArgPairs[i], ""), nameArgPairs[i+1]);
        }
    }

}
