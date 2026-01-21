import javax.swing.*;
import java.awt.*;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.image.BufferStrategy;

/**
 *
 * @author Tom
 */
public class Main extends JFrame implements KeyListener
{
    int x = 0;
    int y = 0;
    int x1 = 0;
    int y1 = 0;
    int x2 = 0;
    int y2 = 0;
    int x3 = 0;
    int y3 = 0;
    int x4 = 0;
    int y4 = 0;
    int x5 = 0;
    int y5 = 0;


    BufferStrategy bs;
    DrawPanel panel = new DrawPanel();

    public Main()
    {
        setIgnoreRepaint(true);
        setTitle("Active Rendering on a JPanel");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setSize(1024,720);
        setResizable(false);
        setVisible(true);
        createBufferStrategy(2);
        bs = getBufferStrategy();
        getContentPane().add(panel);
        panel.setIgnoreRepaint(true);
        addKeyListener(this);
    }

    public void startNow()
    {
        panel.drawStuff();
    }
    public void keyTyped(KeyEvent e)
    {

    }

    public void keyPressed(KeyEvent e)
    {
        System.exit(0);
    }

    public void keyReleased(KeyEvent e)
    {

    }

    public class DrawPanel extends JPanel
    {
        public void drawStuff()
        {
            while(true)
            {
                try
                {
                    x = (int)(Math.round(Math.random()*920));
                    x1 = (int)(Math.round(Math.random()*920));
                    x2 = (int)(Math.round(Math.random()*920));
                    x3 = (int)(Math.round(Math.random()*920));
                    x4 = (int)(Math.round(Math.random()*920));
                    x5 = (int)(Math.round(Math.random()*920));
                    y = (int)(Math.round(Math.random()*600));
                    y1 = (int)(Math.round(Math.random()*600));
                    y2 = (int)(Math.round(Math.random()*600));
                    y3 = (int)(Math.round(Math.random()*600));
                    y4 = (int)(Math.round(Math.random()*600));
                    y5 = (int)(Math.round(Math.random()*600));

                    Graphics2D g = (Graphics2D)bs.getDrawGraphics();
                    g.setColor(Color.BLACK);
                    g.fillRect(0,0,1024,768);
                    g.setColor(Color.red);
                    g.fillRect(x,y,100,100);
                    g.setColor(Color.blue);
                    g.fillRect(x1,y1,100,100);
                    g.setColor(Color.yellow);
                    g.fillRect(x2,y2,100,100);
                    g.setColor(Color.green);
                    g.fillRect(x3,y3,100,100);
                    g.setColor(Color.magenta);
                    g.fillRect(x4,y4,100,100);
                    g.setColor(Color.ORANGE);
                    g.fillRect(x5,y5,100,100);
                    bs.show();
                    Toolkit.getDefaultToolkit().sync();
                    g.dispose();
                    Thread.sleep(20);
                }
                catch (Exception e)
                {
                    System.exit(0);
                }
            }
        }
    }

    public static void main(String[] args)
    {
        Main main = new Main();
        main.startNow();
    }

}