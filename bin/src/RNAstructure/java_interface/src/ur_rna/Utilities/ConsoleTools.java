package ur_rna.Utilities;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

public class ConsoleTools {
    private static OutputStream lastStream = null;
    private static boolean      isFixed    = false;

    private static class FixedStream extends OutputStream {

        private final OutputStream target;

        public FixedStream(OutputStream originalStream) {
            target = originalStream;
        }

        @Override
        public void write(int b) throws IOException {
            if (lastStream!=this) swap();
            target.write(b);
        }

        @Override
        public void write(byte[] b) throws IOException {
            if (lastStream!=this) swap();
            target.write(b);
        }

        @Override
        public void write(byte[] b, int off, int len) throws IOException {
            if (lastStream!=this) swap();
            target.write(b, off, len);
        }

        private void swap() throws IOException {
            if (lastStream!=null) {
                lastStream.flush();
                try { Thread.sleep(200); } catch (InterruptedException e) {}
            }
            lastStream = this;
        }

        @Override public void close() throws IOException { target.close(); }
        @Override public void flush() throws IOException { target.flush(); }
    }

    /**
     * Inserts a 200ms delay into the System.err or System.out OutputStreams
     * every time the output switches from one to the other. This prevents
     * the Eclipse console from showing the output of the two streams out of
     * order. This function only needs to be called once.
     */
    public static void fixConsole() {
        if (isFixed) return;
        isFixed = true;
        System.setErr(new PrintStream(new FixedStream(System.err)));
        System.setOut(new PrintStream(new FixedStream(System.out)));
    }
}