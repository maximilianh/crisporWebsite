package ur_rna.Utilities.swing;

import org.junit.Test;
import ur_rna.Utilities.StopWatch;

import javax.swing.*;
import java.awt.*;

/**
 * @author Richard M. Watson
 */
public class WaitDialogTest {
    @Test
    public void testStopWatch() {
        StopWatch w = new StopWatch(true);
        sleep(10);
        w.println("Test 10");
        sleep(100);
        w.println("Test 110").restart();
        sleep(50);
        w.println("Test 50").restart();
    }

    private boolean sleep(int ms) {
        try {
            Thread.sleep(ms);
            return true;
        } catch (InterruptedException ex) {
            return false;
        }

    }

    @Test
    public void testWaitDialog() {
        WaitDialog wd = new WaitDialog(null, "Banana Title", "Wait for me!", true, true);
        wd.setVisible(true);
        Thread t = new Thread(()->{
            try { Thread.sleep(2000); } catch (InterruptedException ex) { return; }
            for (int i = 0; i < 100; i++) {
                try {
                    final int value = i;
                    Thread.sleep(100);
                    SwingUtilities.invokeLater(() -> wd.setProgressValue(value));
                    SwingUtilities.invokeLater(() -> wd.setMessage("Message #" + value));
                    if (wd.isCanceled()) break;
                } catch (InterruptedException ex) {
                    break;
                }
            }
            SwingUtilities.invokeLater(() -> wd.close());
        });
        t.start();

        JDialog d = new JDialog((Frame)null,"", true);
        d.setSize(200,100);
        d.setVisible(true);
    }
}