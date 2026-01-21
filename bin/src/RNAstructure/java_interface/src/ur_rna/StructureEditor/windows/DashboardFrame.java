package ur_rna.StructureEditor.windows;

import ur_rna.StructureEditor.AppActions;
import ur_rna.StructureEditor.Program;
import ur_rna.StructureEditor.Settings;
import ur_rna.StructureEditor.services.RecentFileList;
import ur_rna.Utilities.swing.ActionHelper;
import ur_rna.Utilities.swing.JFlatButton;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.plaf.basic.BasicInternalFrameUI;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.HierarchyBoundsAdapter;
import java.awt.event.HierarchyEvent;
import java.awt.event.KeyEvent;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * This screen is shown on startup and whenever there are no files loaded.
 * It presents a list of recent and example files, which the user can click to open.
 */
public class DashboardFrame extends JInternalFrame {
    public static final String EXAMPLE_INDICATOR = "example://";
    private JButton btnWebsite;
    private JButton btnHelp;
    private JList<Object> lstRecent;
    private JPanel pnlMain;
    private JList<Object> lstExamples;
    private JButton btnCleanupRecent;
    private JButton btnClearRecent;
    private JCheckBox chkDisableDashboard;
    private JButton btnOpenFile;
    private JButton btnFeedback;
    private boolean updatingSettings;

    //private Map<Object, UiAction> actionMap = new HashMap<>();

    public DashboardFrame() {
        super("Home", false);
        $$$setupUI$$$();
        //actionMap.put(btnWebsite, AppActions.OPEN_WEBSITE);
        //actionMap.put(btnHelp, AppActions.SHOW_LOCAL_HELP);

        btnWebsite.addActionListener(AppActions.OPEN_WEBSITE);
        btnHelp.addActionListener(AppActions.SHOW_LOCAL_HELP);
        btnFeedback.addActionListener(AppActions.SEND_FEEDBACK);

        btnCleanupRecent.setAction(AppActions.CLEANUP_RECENT_FILES);
        btnClearRecent.setAction(AppActions.CLEAR_RECENT_FILES);
        btnOpenFile.setAction(AppActions.SHOW_OPEN_FILE);

        chkDisableDashboard.addActionListener(this::disableOptionClick);
        Program.getInstance().settings().addChangeHandler(this::settingsChanged);

        setPreferredSize(new Dimension(600, 400));
        setBackground(Color.WHITE);
        BasicInternalFrameUI ui = ((BasicInternalFrameUI) this.getUI());
        ui.setNorthPane(null);
        this.setBorder(BorderFactory.createEtchedBorder());
        this.addHierarchyBoundsListener(ancestorResizeListener);
        this.add(pnlMain, BorderLayout.NORTH);

        Cursor handCursor = Cursor.getPredefinedCursor(Cursor.HAND_CURSOR);
        lstExamples.setCursor(handCursor);
        lstRecent.setCursor(handCursor);

        lstRecent.addListSelectionListener(this::listItemSelected);
        lstExamples.addListSelectionListener(this::listItemSelected);

        loadExamples();
        settingsChanged();

        ActionHelper.addKeyAction(this, KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), "closeOnESC",
                e -> Program.getInstance().getMainFrame().hideDashboard());

//        MouseListener ma = new MouseAdapter() {
//            @Override
//            public void mouseClicked(final MouseEvent e) { actionMap.get(e.getSource()).invoke(this); }
//        };
//        for (JComponent c : new JComponent[]{btnHelp, btnWebsite}) {
//            c.addMouseListener(ma);
//            c.setCursor(handCursor);
//        }
        //btnWebsite.setCursor(handCursor);
        //btnHelp.setCursor(handCursor);
    }
    private void settingsChanged() {
        Settings settings = Program.getInstance().settings();
        updatingSettings = true;
        chkDisableDashboard.setSelected(!settings.ShowDashboard);
        updatingSettings = false;
    }
    private void disableOptionClick(final ActionEvent event) {
        if (updatingSettings) return;
        boolean value = chkDisableDashboard.isSelected();
        Program.getInstance().settings().ShowDashboard = !value;
        Program.getInstance().settings().notifyChanged();

        if (value)
            JOptionPane.showMessageDialog(this, "You can re-enable the dashboard screen later from the Settings menu.");
    }
    //    private void disableOptionChanged(final ChangeEvent event) {
//        if (chkDisableDashboard.cha)
//
//    }
    private void listItemSelected(final ListSelectionEvent e) {
        JList list = (JList) e.getSource();
        if (list.getSelectedIndex() == -1) {
            //No selection
        } else {
            //Selection, enable the fire button.
            ListItem li = (ListItem) list.getSelectedValue();
            list.clearSelection();
            AppActions.openFile(li.file.path, li.file.type);
        }
    }
    private HierarchyBoundsAdapter ancestorResizeListener = new HierarchyBoundsAdapter() {
        @Override
        public void ancestorResized(final HierarchyEvent e) {
            //System.out.println("Resized: changed:" + e.getChangedParent() + " parnet: " + getParent());
            updateLayout();
            super.ancestorResized(e);
        }
    };

//    private ComponentListener parentChangeListener = new ComponentAdapter() {
//        @Override
//        public void componentResized(final ComponentEvent e) {
//
//        }
//    };
//    private void hierarchyEvent(final HierarchyEvent e) {
//        if ((e.getChangeFlags() & HierarchyEvent.PARENT_CHANGED) != 0) {
//            if (owner != null)
//                owner.removeComponentListener();
//            System.out.println("Changed Parent: " + e.getChangedParent());
//        }
//    }
//    private Component owner;

    //    @Override
//    public void doLayout() {
//        Container p = this.getParent();
//        if (p != null)
//            updateLayout(p);
//        super.doLayout();
//    }
    public void updateLayout() {
        Container p = this.getParent();
        if (p == null) return;
        Dimension desktopSize = p.getSize(), min = this.getPreferredSize();
        int w = Math.round(desktopSize.width * 0.9f);
        int h = Math.round(desktopSize.height * 0.9f);
        w = Math.min(desktopSize.width, Math.max(min.width, w));
        h = Math.min(desktopSize.height, Math.max(min.height, h));
        reshape((desktopSize.width - w) / 2, (desktopSize.height - h) / 2, w, h);
    }

    private void createUIComponents() {
        btnHelp = new JFlatButton();
        btnWebsite = new JFlatButton();
        btnFeedback = new JFlatButton();
        // TODO: place custom component creation code here
    }

    public static String extractExample(String name) {
        name = name.substring(EXAMPLE_INDICATOR.length());
        File outFile = Paths.get(System.getProperty("java.io.tmpdir"), Program.APPNAME, "examples", name).toFile();
        if (!outFile.exists() || outFile.length() == 0)
            try {
                Files.createDirectories(outFile.toPath().getParent());
                try (InputStream in = Program.getResourceStream("examples/" + name)) {
                    if (in == null) {
                        System.err.println("Resource " + "examples/" + name + " was not found.");
                        return null;
                    }
                    try (FileOutputStream out = new FileOutputStream(outFile)) {
                        byte[] buffer = new byte[1024];
                        int bytesRead;
                        while ((bytesRead = in.read(buffer)) != -1)
                            out.write(buffer, 0, bytesRead);
                        out.flush();
                    }
                }
                if (!outFile.setReadOnly()) System.out.println("Failed to set example file as readonly.");
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        return outFile.toString();
    }
    /**
     * Method generated by IntelliJ IDEA GUI Designer
     * DO NOT edit this method OR call it in your code!
     *
     * @noinspection ALL
     */
    private void $$$setupUI$$$() {
        createUIComponents();
        pnlMain = new JPanel();
        pnlMain.setLayout(new GridBagLayout());
        pnlMain.setBackground(new Color(-1));
        pnlMain.setForeground(new Color(-16777216));
        final JPanel panel1 = new JPanel();
        panel1.setLayout(new GridBagLayout());
        panel1.setBackground(new Color(-14667710));
        panel1.setForeground(new Color(-16777216));
        panel1.setInheritsPopupMenu(false);
        GridBagConstraints gbc;
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 5;
        gbc.weightx = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        pnlMain.add(panel1, gbc);
        panel1.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), null));
        final JLabel label1 = new JLabel();
        label1.setFont(new Font("Segoe UI", Font.BOLD, 14));
        label1.setForeground(new Color(-1));
        label1.setText("RNAstucture is developed by the David H. Mathews Lab at the University of Rochester.");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 4;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.insets = new Insets(10, 0, 5, 0);
        panel1.add(label1, gbc);
        final JLabel label2 = new JLabel();
        label2.setFont(new Font("Segoe UI", Font.BOLD, 12));
        label2.setForeground(new Color(-1));
        label2.setText("The RNA Drawing Editor was written by Richard M. Watson");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.gridwidth = 4;
        gbc.anchor = GridBagConstraints.NORTH;
        gbc.insets = new Insets(2, 0, 5, 0);
        panel1.add(label2, gbc);
        btnHelp.setActionCommand("help");
        btnHelp.setAlignmentY(0.0f);
        btnHelp.setFont(new Font("Segoe UI", Font.BOLD, 14));
        btnHelp.setForeground(new Color(-256));
        btnHelp.setIcon(new ImageIcon(getClass().getResource("/ur_rna/StructureEditor/resources/images/open-link_yellow.png")));
        btnHelp.setInheritsPopupMenu(true);
        btnHelp.setText("Help/Documentation");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 2;
        gbc.weightx = 0.3;
        gbc.insets = new Insets(5, 0, 5, 0);
        panel1.add(btnHelp, gbc);
        btnWebsite.setActionCommand("website");
        btnWebsite.setAlignmentY(0.0f);
        btnWebsite.setFont(new Font("Segoe UI", Font.BOLD, 14));
        btnWebsite.setForeground(new Color(-256));
        btnWebsite.setIcon(new ImageIcon(getClass().getResource("/ur_rna/StructureEditor/resources/images/open-link_yellow.png")));
        btnWebsite.setText("Lab Website");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 2;
        gbc.weightx = 0.3;
        gbc.insets = new Insets(5, 0, 5, 0);
        panel1.add(btnWebsite, gbc);
        btnFeedback.setActionCommand("feedback");
        btnFeedback.setAlignmentY(0.0f);
        btnFeedback.setFont(new Font("Segoe UI", Font.BOLD, 14));
        btnFeedback.setForeground(new Color(-256));
        btnFeedback.setIcon(new ImageIcon(getClass().getResource("/ur_rna/StructureEditor/resources/images/open-link_yellow.png")));
        btnFeedback.setText("Send Feedback");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 2;
        gbc.weightx = 0.3;
        gbc.insets = new Insets(5, 0, 5, 0);
        panel1.add(btnFeedback, gbc);
        final JLabel label3 = new JLabel();
        label3.setFont(new Font("Segoe UI", Font.BOLD, 14));
        label3.setForeground(new Color(-16777216));
        label3.setText("Recent Files:");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 0.1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(10, 5, 0, 0);
        pnlMain.add(label3, gbc);
        lstRecent = new JList();
        lstRecent.setForeground(new Color(-16776961));
        lstRecent.setMinimumSize(new Dimension(100, 10));
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.gridwidth = 5;
        gbc.weighty = 0.9;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.VERTICAL;
        gbc.insets = new Insets(0, 10, 0, 0);
        pnlMain.add(lstRecent, gbc);
        final JLabel label4 = new JLabel();
        label4.setFont(new Font("Segoe UI", Font.BOLD, 14));
        label4.setForeground(new Color(-16777216));
        label4.setText("Example Files:");
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.gridwidth = 5;
        gbc.weighty = 0.1;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(10, 5, 0, 0);
        pnlMain.add(label4, gbc);
        lstExamples = new JList();
        lstExamples.setForeground(new Color(-16776961));
        gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 4;
        gbc.gridwidth = 5;
        gbc.weighty = 0.9;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.VERTICAL;
        gbc.insets = new Insets(0, 10, 0, 0);
        pnlMain.add(lstExamples, gbc);
        btnCleanupRecent = new JButton();
        btnCleanupRecent.setBackground(new Color(-1));
        btnCleanupRecent.setHideActionText(true);
        btnCleanupRecent.setIcon(new ImageIcon(getClass().getResource("/ur_rna/StructureEditor/resources/images/clear-recent.png")));
        btnCleanupRecent.setText("");
        btnCleanupRecent.setToolTipText("Cleanup missing files");
        gbc = new GridBagConstraints();
        gbc.gridx = 2;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlMain.add(btnCleanupRecent, gbc);
        chkDisableDashboard = new JCheckBox();
        chkDisableDashboard.setBackground(new Color(-1));
        chkDisableDashboard.setForeground(new Color(-16777216));
        chkDisableDashboard.setMargin(new Insets(2, 14, 2, 14));
        chkDisableDashboard.setText("Disable this Screen");
        chkDisableDashboard.setToolTipText("This dashboard screen is usually shown when no files are open, unless you check this box which will disable it.");
        gbc = new GridBagConstraints();
        gbc.gridx = 4;
        gbc.gridy = 1;
        gbc.anchor = GridBagConstraints.WEST;
        pnlMain.add(chkDisableDashboard, gbc);
        btnClearRecent = new JButton();
        btnClearRecent.setBackground(new Color(-1));
        btnClearRecent.setHideActionText(true);
        btnClearRecent.setIcon(new ImageIcon(getClass().getResource("/ur_rna/StructureEditor/resources/images/delete.png")));
        btnClearRecent.setText("");
        btnClearRecent.setToolTipText("Clear ALL recent files");
        gbc = new GridBagConstraints();
        gbc.gridx = 3;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlMain.add(btnClearRecent, gbc);
        btnOpenFile = new JButton();
        btnOpenFile.setBackground(new Color(-1));
        btnOpenFile.setHideActionText(true);
        btnOpenFile.setIcon(new ImageIcon(getClass().getResource("/ur_rna/StructureEditor/resources/images/file-open.png")));
        btnOpenFile.setText("");
        btnOpenFile.setToolTipText("Cleanup missing files");
        gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        pnlMain.add(btnOpenFile, gbc);
    }
    /** @noinspection ALL */
    public JComponent $$$getRootComponent$$$() { return pnlMain; }

    private static class ListItem {
        public final String desc;
        public final RecentFileList.RecentFile file;
        public ListItem(final RecentFileList.RecentFile file, final String desc) {
            this.desc = desc;
            this.file = file;
        }
        @Override
        public String toString() { return desc; }
    }
    public void setRecent(final RecentFileList items) {
        int i = 0;
        Object[] listItems = new ListItem[items.size()];
        for (RecentFileList.RecentFile rf : items) {
            File f = new File(rf.path);
            ListItem item = new ListItem(rf, String.format("%s    in    %s",
                    f.getName(),
                    f.getParent()
            ));
            listItems[i++] = item;
        }
        lstRecent.setListData(listItems);
    }
    public void loadExamples() {
        try {
            try (InputStream stream = Program.getResourceStream("examples/examples.dir")) {
                if (stream == null) return;
                java.util.List<ListItem> files = new ArrayList<>();

                try (Scanner scanner = new Scanner(stream)) {
                    while (scanner.hasNextLine()) {
                        String line = scanner.nextLine().trim();
                        if (line.isEmpty()) continue;
                        RecentFileList.RecentFile r = new RecentFileList.RecentFile(EXAMPLE_INDICATOR + line);
                        files.add(new ListItem(r, line));
                    }
                }
                lstExamples.setListData(files.toArray(new ListItem[files.size()]));
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    //    @Override
//    public void dispose() {
//        super.dispose();
//    }
}
