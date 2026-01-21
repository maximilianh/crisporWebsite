 package ur_rna.GUITester.GuiTools;

 import abbot.finder.AWTHierarchy;
 import ur_rna.GUITester.GuiTools.Matchers.DescriptionMatcher;
 import ur_rna.GUITester.RuntimeTools;
 import ur_rna.RNAstructureUI.AppMainFrame;
 import ur_rna.RNAstructureUI.RNAstructure;
 import ur_rna.Utilities.AppLog;

 import javax.accessibility.AccessibleContext;
 import javax.swing.*;
 import javax.swing.text.JTextComponent;
 import java.awt.*;

 public class GuiAppManager {
	private final AppLog _log;
	private GuiFinder _finder;
	private GuiRobot _robot;
	private JFrame _rootFrame = null;
	public JFrame getRootFrame() { return _rootFrame; }

	public GuiAppManager(AppLog log) {
		_log = log;
		_finder = new GuiFinder(this);
		_robot = new GuiRobot(this);
	}

	public boolean isShowing() { return isStarted() && _rootFrame.isShowing();}
	public boolean isStarted() { return _rootFrame != null && _rootFrame.isDisplayable();}
	public GuiFinder guiFinder() { return _finder; }
	public GuiRobot guiRobot() { return _robot; }
	public AppLog getLog() { return _log; }


	/**
	 * Close the GUI application including all existing Windows.
	 */
	public void close() {
		try {
			_robot.waitForIdle();
			if (_rootFrame != null)
				_rootFrame.dispose();
			_rootFrame = null;
			for (Window w : java.awt.Window.getWindows()) {
				if (!"Abbot Robot Verification".equals(w.getName()))
					w.dispose();
			}
		} catch (Exception e) {
			_log.warn("Exception when closing application.", e);
		}
	}
	
//	/**
//	 * Start the application. If another instance is already showing the behavior depends on {@code preserveExisting}.
//	 * @param preserveExisting If true, an existing visible application will be re-used instead of the default behavior, which is to close any existing instance and start a new one. This can improve performance, but it may also result in unwanted side-effects if the state of the application has been altered by a previous test.
//	 */
//	public void start(boolean preserveExisting) {
//		if (preserveExisting && isShowing())
//			return;
//	}

	/**
	 * Start the GUI application. If another instance is already showing, close it and start a new one.
	 */
	public void start(String args) {
		start(args.split("\\s+"));
	}

	public void start(String[] args) {
		close();
		launchSubjectApp(args);
	}

	private void launchSubjectApp(String[] args) {
		try {
			RNAstructure.main(args);
			_rootFrame = AppMainFrame.getFrame();

//			if (_log.isTraceEnabled())
//				printHierarchyTrace();
		} catch (Exception e) {
			throw new RuntimeException("Failed to create or detect new RNAstructure application window.", e);
		}
	}

	public void printHierarchyTrace() {
		if (_log.isTraceEnabled())
			_log.trace(getHierarchy(null).toString());
	}

	public StringBuilder getHierarchy(java.util.Collection roots) {
		if (roots == null) roots = AWTHierarchy.getDefault().getRoots();
		StringBuilder sb = new StringBuilder();
		String prefix = "\n  ";
		sb.append("Application Hierarchy :").append(roots.size()).append(" roots.");
		for ( Object root : roots ) {
			Component rc = (Component) root;
			sb.append("\n  +Root: ");
			addComponentInfo(rc, sb);
			dumpHierarchy(sb, ((Container)rc), prefix);
//			if (rc instanceof Window)
//				dumpHierarchy(sb, ((Window)rc).getOwnedWindows(), prefix);
//			if (rc instanceof Container)
//				dumpHierarchy(sb, ((Container)rc).getComponents(), prefix);
//			dumpHierarchy(sb, ((Container)rc).getComponents(), prefix);
		}
		return sb;
	}

	public void dumpHierarchy(StringBuilder sb, Container parent, String prefix) {
//		while (prefixes.size() <= level)
//			prefixes.add(prefixes.get(prefixes.size()-1) + "|  ");
//		String prefix = prefixes.get(level++);
		String childPrefix = prefix + "|  ";
		String lastChildPrefix = childPrefix; //prefix.replace("|", " ") + "|   ";

		int item = 0;
		if (parent instanceof Window) {
			Window[] windows = ((Window)parent).getOwnedWindows();
			int lastWindow = windows.length - 1;
			for (int i = 0; i <= lastWindow; i++) {
				Window c = windows[i];
				sb.append(prefix).append("|--")
						.append("[").append(item++).append("] ");
				addComponentInfo(c, sb);
				dumpHierarchy(sb, c, (i == lastWindow) ? childPrefix : lastChildPrefix);
			}
		}

		Component[] list = parent.getComponents();
		int lastChild = list.length -1;
		for (int i = 0; i <= lastChild; i++) {
			Component c = list[i];
			sb.append(prefix).append("|--")
					.append("[").append(item++).append("] ");
			addComponentInfo(c, sb);
			if (c instanceof Container)
				dumpHierarchy(sb, (Container)c, (i == lastChild)?childPrefix:lastChildPrefix);
		}
	}

	private void appendText(StringBuilder sb, String propertyName, String text) {
		if (text != null && !text.isEmpty())
			sb.append(" ").append(propertyName).append(": \"").append(text).append("\"");
	}
	public String getComponentInfo(Component c) {
		StringBuilder sb = new StringBuilder();
		addComponentInfo(c, sb);
		return sb.toString();
	}
	 public void addComponentInfo(Component c, StringBuilder sb) {
		 sb.append(c.getClass().getSimpleName());
		 appendText(sb, "Name", c.getName());
		 if (c instanceof JComponent) appendText(sb, "ToolTip", ((JComponent) c).getToolTipText());
		 if (c instanceof JTextComponent) appendText(sb, "Text", ((JTextComponent) c).getText());
	 	 if (c instanceof TextComponent) appendText(sb, "Text", ((TextComponent) c).getText());
		 if (c instanceof AbstractButton) appendText(sb, "Caption", ((AbstractButton) c).getText());
		 if (c instanceof Button) appendText(sb, "Caption",(((Button) c).getLabel()));
		 if (c instanceof JLabel) appendText(sb, "LabelText", ((JLabel) c).getText());
		 if (c instanceof Label) appendText(sb, "LabelText", ((Label) c).getText());
		 if (c instanceof Frame) appendText(sb, "Title", ((Frame) c).getTitle());
		 if (c instanceof Dialog) appendText(sb, "Title", ((Dialog) c).getTitle());
		 if (c instanceof JMenuItem) appendText(sb, "MenuPath", DescriptionMatcher.getMenuPath((JMenuItem) c));
		 if (c instanceof JList) appendText(sb, "Selected", ""+((JList) c).getSelectedValue());
		 if (c instanceof JComboBox) appendText(sb, "Selected", "" + ((JComboBox) c).getSelectedItem());
		 appendText(sb, "Visible", Boolean.toString(c.isVisible()));
		 AccessibleContext ac = c.getAccessibleContext();
		 if (ac != null) {
			appendText(sb, "AccessibleName", ac.getAccessibleName());
			appendText(sb, "AccessibleDesc", ac.getAccessibleDescription());
//			if (ac.getAccessibleText() != null)
//				appendText(sb, "AccessibleText", ac.getAccessibleText().get .toString());
		}
	}
	 /**
	  * Peform a named action on a GuiRef.
	  * @param action
	  * @param args
	  * @throws RuntimeTools.ActionException
	  * @throws GuiItemNotFoundException
      */
	 public void performGuiAction(final String action, final Object[] args)
			 throws RuntimeTools.ActionException, GuiItemNotFoundException {
		 if (!isStarted()) throw new RuntimeTools.ActionException("GUI components cannot be accessed because no window exists.");
		 _robot.doAction(action, args);
	 }
	 /**
	  * @see GuiRobot#finalShutdown()
	  */
	 public void shutdown() {
		 GuiRobot.finalShutdown();
	 }
 }
