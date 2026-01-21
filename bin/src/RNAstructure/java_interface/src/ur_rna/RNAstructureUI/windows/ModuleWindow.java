/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI.windows;

import ur_rna.RNAstructureUI.RNAstructure;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.ui.ProgressPanel;
import ur_rna.RNAstructureUI.utilities.UserOptions;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.PathTools;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyVetoException;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import static ur_rna.Utilities.PathTools.isDir;

/**
 * A class responsible for creating and displaying a window that enables a
 * user to input data for one of the RNAstructure modules.
 *
 * @author Jessica S. Reuter
 */
public abstract class ModuleWindow
	extends InternalWindow {
	private static final long serialVersionUID = 20120802;
    protected final ModuleWindow self = this;  // for use in nested and anonymous classes instead of ModuleWindow.this
    protected static final AppLog log = RNAstructure.log;

	/**
	 * The GridBagConstraints for this window.
	 */
	private GridBagConstraints constraints;

	/**
	 * The layout for this window.
	 */
	private GridBagLayout layout;

	/**
	 * A boolean, true if the calculation handles RNA, false if not.
	 */
	protected boolean isRNA;

	/**
	 * The progress bar dialog for this window.
	 */
	protected volatile Component progressIndicator;
	protected Thread mainCalcThread;
	protected JButton startButton;
	private ActionListener calcCancelListener = new ActionListener() {
		@Override
		public void actionPerformed(final ActionEvent e) {
			cancelMainCalculation();
		}
	};
    private long calcCompleted; // always set to 0 before runMainCalculation.
    private ActionListener calcProgressListener = new ActionListener() {
        @Override
        public void actionPerformed(final ActionEvent e) {
            Component c = progressIndicator;
            if (c instanceof ProgressPanel)
                updateCalcProgress((ProgressPanel) c);
        }
    };

    long lastUpdate;
    protected void updateCalcProgress(ProgressPanel panel) {
        int p = backend.getProgressNumber();
        if (p >= 100) {
            // Some calculations complete (100%) BEFORE they are actually done.
            // e.g. AllSub reports 100% before any tracebacks are performed,
            // so the calculation is really only about 30% done.
            // To deal with this, switch to indeterminate mode at about 1 second after completion.
            if (calcCompleted == 0) {
                calcCompleted = System.currentTimeMillis();
                panel.setProgress(100);
                panel.setAllowCancel(false);
            } else if (System.currentTimeMillis() - calcCompleted > 1000) {
                panel.setIndeterminate(true); // the calculation has completed, but some operations are still being peformed, delaying the result.
                panel.setMessage("Performing final steps...");
            }
        } else
            panel.setProgress(p);
    }

    private boolean calcCanceled;
    protected boolean calcWasCanceled() {
        return calcCanceled;
    }
	protected void cancelMainCalculation() {
        calcCanceled = true;
        backend.cancelOperation();
//		if (mainCalcThread != null)
//			mainCalcThread.interrupt();
//		closeProgress();
//		if (startButton != null)
//			startButton.setEnabled(true);
	}

    protected void closeProgress() { closeProgress(true); }
	protected void closeProgress(final boolean resetUI) {
        final Component c;
        synchronized (this) {
            // only one thread can capture progressIndicator. All subsequent calls will be null.
            c = progressIndicator;
            progressIndicator = null;
        }
        if (c != null)
            runOnUiThread(new Runnable() {
                @Override
                public void run() {
                    enableControls(self, true);
                    if (c instanceof ProgressPanel)
                        ((ProgressPanel) c).stop();
                    self.remove(c);
                    if (resetUI) self.pack();
                }
            });
    }

    void runOnUiThread(Runnable r) {
        if (SwingUtilities.isEventDispatchThread())
            r.run();
        else
        try {
            SwingUtilities.invokeAndWait(r);  // throws InterruptedException, InvocationTargetException
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private HashMap<Component, Boolean> enabledState = new HashMap<>();
    private void enableControls(Container parent, final boolean enable) {
        for (Component c : parent.getComponents()) {
            if (!enable) {
                enabledState.put(c, c.isEnabled());
                c.setEnabled(false);
            } else {
                if (enabledState.get(c) == Boolean.TRUE) // value can be TRUE, FALSE, or null. Only enable if identically TRUE.
                    c.setEnabled(true);
            }
            if (c instanceof Container)
                enableControls((Container) c, enable);
        }
    }


    /**
	 * Constructor.
	 *
	 * @param acid    The nucleic acid type, as a string.
	 * @param title   The title of the window.
	 */
	protected ModuleWindow( String acid, String title ) {

		// Initialize whether this window handles RNA or DNA.
		isRNA = ("RNA".equalsIgnoreCase(acid));

		// Initialize the layout and constraints.
		layout = new GridBagLayout();
		setLayout( layout );
		constraints = new GridBagConstraints();

		// Build the window.
		setCaption( title );
		makeInputControls();
	}

	/**
	 * Combine two file names into a new file name.
	 *
	 * @param file1   The first file name.
	 * @param file2   The second file name.
	 * @param ext     The extension of the new file name.
	 * @return        The combined file name.
	 */
	protected final String combineFileNames(
		String file1, String file2, String ext ) {
		String separator = File.separator;
		String name1 = new File( file1 ).getAbsolutePath();
		String name2 = new File( file2 ).getAbsolutePath();
		String dir = name1.substring( 0, name1.lastIndexOf( separator ) + 1 );
		name1 = name1.substring(
			name1.lastIndexOf( separator ) + 1,
			name1.lastIndexOf( "." ) );
		name2 = name2.substring(
			name2.lastIndexOf( separator ) + 1,
			name2.lastIndexOf( "." ) );
		return dir + name1 + "_" + name2 + "." + ext;
	}


	/**
	 * Perform an action specific to this window.
	 * @param ci        Information about the command.
	 */
	@Override
	protected void processCommand(final CommandInfo ci) {
		// If the command comes from the "START" button, run the main
		// calculation in a new background thread.
        if( ci.command.equals( "START" ) ) {
            calcCompleted = 0;
            calcCanceled = false;
			mainCalcThread = new Thread() {
				public void run() {
                    if (runMainCalculation())
                        self.close();
                    else
                        closeProgress();
                }
			};
			mainCalcThread.start();
		} else
			super.processCommand(ci);
	}

	/**
	 * Handle finishing of a module.
	 * <br><br>
	 * Finishing of a module involves closing of the window and showing any
	 * error that may have occurred.
	 *
	 * @param text   The result of the calculation spawned by the module.
	 */
	protected final boolean displayCalcError( String text ) {
        closeProgress();
        if(text == null || text.isEmpty())
            return true;

        if (backend.wasCanceled())
            Dialogs.showMessage( "Calculation cancelled." );
        else
            Dialogs.showError( text );
        return false;
	}

    protected boolean drawStructures(final String file) {
//        if (!SwingUtilities.isEventDispatchThread()) {
//            RunFunc r = new RunFunc(file) { public Object run(Object param) { return drawStructures((String)param); } };
//            try {
//                SwingUtilities.invokeAndWait(r);
//            } catch (Exception ex) {
//                ex.printStackTrace();
//            }
//            return r.getBoolResult();
//        }
        if (promptDraw()) {
            showDrawProgress();
            DrawingWindow drawing = new DrawingWindow( file );
            if(drawing.isError())
                return false;
            drawing.showWindow();
        }
        return true;
    }

	/**
	 * Get a particular input control from the content pane.
	 *
	 * @param index   The index of the control to get, one-indexed.
	 * @return        The input control.
	 */
	protected final JComponent getInputControl( int index ) {
		return getInputControl( getContentPane(), index );
	}

	/**
	 * Get a particular input control from a given container.
	 *
	 * @param hold    The container that holds the input control.
	 * @param index   The index of the control to get, one-indexed.
	 * @return        The input control.
	 */
	protected final JComponent getInputControl( Container hold, int index ) {
		return (JComponent)hold.getComponent( index - 1 );
	}

	/**
	 * Make any input controls necessary to get data from this window for a
	 * calculation.
	 */
	protected abstract void makeInputControls();

	/**
	 * Create the "START" button.
	 * <br><br>
	 * Note that this method simply makes the button, adds the action
	 * listener, applies standard padding, then places the button at the
	 * specified (X,Y) layout coordinate.
	 * <br><br>
	 * If more customization of the button itself or its positioning is
	 * needed, it should be done in subclasses.
	 *
	 * @param x   The X index at which to place the button in the layout.
	 * @param y   The Y index at which to place the button in the layout.
	 */
	protected final void makeStartButton( int x, int y ) {
		startButton = new JButton( "START" );
		startButton.addActionListener( this );
		setPad( 15, 15 );
		setInsets( 10, 10, 10, 10 );
		setFillCenter();
		placeComponent( x, y, startButton );
	}

	/**
	 * Check to see if the module was initialized correctly.
	 * <br><br>
	 * Check what the result string retrieved from the module was, and return
	 * true only if it's the empty string. Otherwise, show the result string
	 * as an error and return false.
	 *
	 * @param result   The result of the data structure construction.
	 * @return         True if the data structure was initialized correctly,
	 *                 false if not.
	 */
	protected final boolean verifyBackendResult(String result, String details, Object... detailArgs) {
		if( result.equals( "" ) ) { return true; }
		else {
			if (details != null && details.length() != 0) {
				result += "\n\nAdditional Details:\n" + String.format(details, detailArgs);
			}
			Dialogs.showError( result );
			return false;
		}
	}

	/**.
	 * Set a component in its place.
	 *
	 * @param x           The X coordinate of the grid.
	 * @param y           The Y coordinate of the grid.
	 * @param component   The component to set in place.
	 */
	protected final void placeComponent( int x, int y, JComponent component ) {
		constraints.gridx = x;
		constraints.gridy = y;
		add( component, constraints );
	}

	private static HashMap<String, Boolean> knownDirs = new HashMap<>();
    /**
	 * Replace the extension on a file name.
	 * @param file        The file name whose extension to replace.
	 * @param extension   The new extension.
	 * @return            The new file name with extension replaced.
	 */
    public static String getOutputFile( String file, String extension ) {
        return getOutputFile(PathTools.changeExtension(file, extension));
    }
	protected static String protectedExamplesDir;
    public static String getOutputFile( String file ) {
		// Verify the target directory to ensure that it is writeable.
		String dir;
		File f = new File(file);
		try {
			if (protectedExamplesDir == null) {
				File examples = PathTools.findLocalPath("examples");
				if (examples == null)
					protectedExamplesDir = "<>"; //invalid.
				else
					protectedExamplesDir = examples.getCanonicalFile().toString();
			}

			dir = f.getCanonicalFile().getParent();
			if (dir != null) {
				Boolean write = knownDirs.get(dir);
				if (write == null) {
					File parent = f.getParentFile();
					if (parent == null) parent = new File(".");
					write = PathTools.verifyWritable(parent);
					knownDirs.put(dir, write);
				}
				if (!write) dir = null;
			}
		}catch (IOException ex) {
			dir = null;
		}

		if (dir == null || dir.equals(protectedExamplesDir)) {
			File home = PathTools.getHomeDir();
			// The home directory is something like C:\Users\<USERNAME> on windows  or /home/<USERNAME> on Linux or /Users/<USERNAME> on Mac
			// But we'd like to put the file inside the Documents sub-folder if it exists.
			if (!"documents".equalsIgnoreCase(home.getName())) {
				if (isDir(home + "/Documents"))
					home = new File(home + "/Documents");
				else if (isDir(home + "/documents"))
					home = new File(home + "/documents");
			}
			dir = home.toString();
		} else
			dir = f.getParent(); //  previously dir was canonical

		return dir + File.separator + f.getName();
	}

	/**
	 * Run the module calculation in the back end.
	 * <br><br>
	 * Data must be pulled from the module window before the back end
	 * calculation is run, so this data capture should be done in this method.
	 */
	protected abstract boolean runMainCalculation();

	/**
	 * Set the anchoring area for the GridBagConstraints to CENTER.
	 */
	protected final void setAnchorCenter() {
		constraints.anchor = GridBagConstraints.CENTER;
	}

	/**
	 * Set the anchoring area for the GridBagConstraints to NORTH.
	 */
	protected final void setAnchorNorth() {
		constraints.anchor = GridBagConstraints.NORTH;
	}

	/**
	 * Set the fill for the GridBagConstraints as CENTER.
	 */
	protected final void setFillCenter() {
		constraints.fill = GridBagConstraints.CENTER;
	}

	/**
	 * Set the fill for the GridBagConstraints as HORIZONTAL.
	 */
	protected final void setFillHorizontal() {
		constraints.fill = GridBagConstraints.HORIZONTAL;
	}

	/**
	 * Set the amount of space a component takes up in the grid.
	 *
	 * @param width    The width of the component.
	 * @param height   The height of the component.
	 */
	protected final void setGrid( int width, int height ) {
		constraints.gridwidth = width;
		constraints.gridheight = height;
	}

	/**
	 * Set the external padding of the constraints.
	 *
	 * @param top      The top padding.
	 * @param left     The left padding.
	 * @param bottom   The bottom padding.
	 * @param right    The right padding.
	 */
	protected final void setInsets(
		int top, int left, int bottom, int right ) {
		constraints.insets = new Insets( top, left, bottom, right );
	}

	/**
	 * Set the internal padding of the constraints.
	 *
	 * @param xPad   The padding in the X direction.
	 * @param yPad   The padding in the Y direction.
	 */
	protected final void setPad( int xPad, int yPad ) {
		constraints.ipadx = xPad;
		constraints.ipady = yPad;
	}

	/**
     * Create a determinate progress bar.
     */
    protected void showProgress() { showProgress( true ); }
	/**
	 * Create a progress bar.
	 *
	 * @param determinate   True if the progress bar is determinate,
	 *                      false if not.
	 */
    protected void showProgress( boolean determinate ) { showProgress( determinate, false ); }
    protected void showProgress( boolean determinate, boolean canCancel) {  showProgress(determinate, canCancel, null); }
    protected void showProgress( boolean determinate, boolean canCancel, String message) {
        closeProgress(false); //in case a progress indicator is already showing
        enableControls(this, false);
        // prevent any other thread from calling closeProgress while this is being created.

        if (message == null) message = determinate ? "Calculation progress: " : "Calculating...";

        // Add the progress bar panel to a dialog.
        ProgressPanel p = new ProgressPanel(message,
                canCancel ? calcCancelListener : null,
                calcProgressListener
        ); //new JPanel()
        p.setIndeterminate(!determinate);
        p.setProgress(0);
        p.start();

        // Add the ProgressPanel to the bottom of the current window.
        // Get the number of layout rows and columns
        int[][] dims = layout.getLayoutDimensions();
        //	Note: 	dim[0] = new int[layoutInfo.width];
        //	Note: 	dim[1] = new int[layoutInfo.height];

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(2, 2, 2, 2);
        gbc.ipadx = 5;
        gbc.ipady = 5;
        gbc.gridwidth = dims[0].length;
        gbc.gridy = dims[1].length;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        p.setBorder(BorderFactory.createEtchedBorder());
        this.add(p, gbc);
        this.pack();
        progressIndicator = p;
    }

	public static boolean promptDraw() { return promptDraw(null, null); }
    public static boolean promptDraw(String prompt, String promptSetting) {
		if (prompt == null) prompt = "Do you want to draw structures?";
        if (promptSetting == null) promptSetting = RNAstructure.options.drawStructures;
        switch(promptSetting) {
			case UserOptions.DRAW_ALWAYS:
				return true;
			case UserOptions.DRAW_NEVER:
				return false;
			default: // PROMPT
				return Dialogs.showConfirm(prompt);
		}
	}
    protected void showDrawProgress() { showDrawProgress(null);  }
    protected void showDrawProgress(String message) {
        if (message == null) message = "Calculation Complete.  Now Drawing...";
        showProgress(false, false, message);
    }

    public void close() {
        try {
            this.setClosed(true);
        } catch (PropertyVetoException ex) {
            // OK fine. Don't close. Whatever.
        }
    }
}
