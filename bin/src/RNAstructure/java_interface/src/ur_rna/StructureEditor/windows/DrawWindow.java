package ur_rna.StructureEditor.windows;

import ur_rna.RNAstructure.RnaBackendException;
import ur_rna.StructureEditor.*;
import ur_rna.StructureEditor.menus.FileMenu;
import ur_rna.StructureEditor.models.*;
import ur_rna.StructureEditor.services.RnaDrawController;
import ur_rna.StructureEditor.services.SceneColorizer;
import ur_rna.StructureEditor.services.SceneUpdateEvent;
import ur_rna.StructureEditor.services.drawing.NucLayout;
import ur_rna.StructureEditor.services.drawing.View2D;
import ur_rna.StructureEditor.services.drawing.export.Graphics2DRecorder;
import ur_rna.StructureEditor.services.drawing.export.PageSize;
import ur_rna.StructureEditor.services.drawing.export.SvgGraphicsExporter;
import ur_rna.StructureEditor.services.fileIO.RnaFileIO;
import ur_rna.StructureEditor.ui.DrawPanel;
import ur_rna.Utilities.Colors;
import ur_rna.Utilities.PathTools;
import ur_rna.Utilities.StopWatch;
import ur_rna.Utilities.swing.*;
import ur_rna.Utilities.swing.FileDialog;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.Timer;
import javax.swing.border.BevelBorder;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
import java.awt.*;
import java.awt.event.AdjustmentListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.prefs.Preferences;

import static ur_rna.StructureEditor.Program.getIcon;
import static ur_rna.Utilities.Strings.*;
import static ur_rna.Utilities.swing.AcceleratorKey.resetTabTraversalKeys;

/**
 * Basic child window for drawing.
 */
public class DrawWindow extends ChildFrame {
    private JScrollPane scrollPane;
    private DrawPanel panel;
    public final RnaDrawController controller = new RnaDrawController();
    private List<MergeMenu> menus = new ArrayList<>();
    private MergeMenu structureMenu = new MergeMenu("Choose &Structure"), formatBondType;
    private ButtonGroupList structureList = new ButtonGroupList();
    private ToolButtonList toolButtons = new ToolButtonList();
    private RnaSceneGroup scenes;
    private RnaScene currentScene;
    private SceneInfo[] sceneInfo;
    private NucLayout layout;
    private JButton structureMenuToolbarButton;
    //private AbstractButton btnModeSingle, btnModeBranch;
    private AbstractButton btnSelModeNuc, btnSelModeHelix, btnSelModeBranch;
    private DrawSettings drawSettings;
    private Timer asyncUpdateTimer;
    private JPanel pnlTaskProgress;
    private JLabel lblTaskStatus;
    private JProgressBar prgTaskProgress;
    private JLabel lblTaskProgress;
    private RnaFileIO.AsyncTask _asyncTask;

//    public void setBehaviorMode(final RnaDrawController.BehaviorMode behaviorMode) {
//        controller.setBehaviorMode(behaviorMode);
//        // We SHOULD be able to affect the selected state of the buttons just by changing the selected state of the actions.
//        // Unfortunately, the authors of JButton forgot to override the shouldUpdateSelectedStateFromAction function, so this won't work.
//        Action mode = behaviorMode == RnaDrawController.BehaviorMode.BRANCH ? actionModeBranch : actionModeSingle;
//        if (!UiAction.isSelected(mode))
//            UiAction.setSelected(mode, true);
//    }

    public void setSelMode(final RnaDrawController.SelectionType selMode) {
        controller.setSelType(selMode);
        Action a;
        switch (selMode) {
            case Branch: a = actionSelectBranch; break;
            case HelixOrLoop: a = actionSelectHelix; break;
            case Individual: a = actionSelectNuc; break;
            default: a = null; break;
        }
        // We SHOULD be able to affect the selected state of the buttons just by changing the selected state of the actions.
        // Unfortunately, the authors of JButton forgot to override the shouldUpdateSelectedStateFromAction function, so this won't work.
        if (a==null) {
            UiAction.setSelected(a, false);
        }
        if (!UiAction.isSelected(a))
            UiAction.setSelected(a, true);
    }

    public RnaSceneGroup getScenes() {
        return scenes;
    }
    public RnaScene getCurrentScene() {
        return currentScene;
    }

    private static class SceneInfo {
        Point2D offset;
        float scale;
        public Point scrollPos;
    }

    public static class ToolButtonList extends ArrayList<Component> {
        private int mergePos = IMergeItem.MergePosDefault;
        public void setMergePos(int pos) { mergePos = pos; }
        public <T extends JComponent> T add(T component) { super.add(component); return component; }
        public MergeButton add(Action a) {            return add(new MergeButton(a, mergePos));        }
        public MergeToggleButton addToggle(Action a) {            return add(new MergeToggleButton(a, mergePos));        }
        public JSeparator addSeparator() {            return add(new MergeButton.Separator(mergePos));        }
    }

    public void loadFile(final RnaSceneGroup file) {
        this.scenes = file;
        sceneInfo = new SceneInfo[this.scenes.size()];
        buildStructureMenu();
        loadScene(scenes.get(0));
        //controller.addUndo(SceneUpdateInfo.View);
    }

    RnaDrawController.SelectionType selType = RnaDrawController.SelectionType.Individual;
    // File
    UiAction actionSave = new UiAction("&Save Drawing", e->save(FileType.NsdDrawing, false, false), getIcon("save-drawing")).setKeyStroke("*S");
    UiAction actionSaveAs = new UiAction("&Save As..", e->save(FileType.NsdDrawing, true, false), getIcon("save-as")).setKeyStroke("*+S");
    UiAction actionExport = new UiAction("&Export (Image, FASTA, CT, etc)", e->save(FileType.AnySaveable, true, true), getIcon("export-as")).setKeyStroke("*+E");
    UiAction actionClose = new UiAction("Close &Window", this::close).setKeyStroke("*W");
    // Edit
    UiAction actionEditUndo = new UiAction("&Undo", controller::undo, getIcon("edit-undo")).setKeyStroke("*Z");
    UiAction actionEditRedo = new UiAction("&Redo", controller::redo, getIcon("edit-redo")).setKeyStroke("*Y");

    UiAction actionSelectAll = new UiAction("Select &All",controller::selectAll).setKeyStroke("*A");
    UiAction actionSelectExpand = new UiAction("&Expand Selection",controller::expandSelection, getIcon("edit-sel-expand")).setKeyStroke(" ");
    UiAction actionSelectNuc = new UiAction("Select &Individual Nucleotides",e->setSelMode(RnaDrawController.SelectionType.Individual), getIcon("edit-sel-single")).setKeyStroke("N");
    UiAction actionSelectHelix = new UiAction("Select &Helices and Loops",e->setSelMode(RnaDrawController.SelectionType.HelixOrLoop), getIcon("edit-sel-helix")).setKeyStroke("H");
    UiAction actionSelectBranch = new UiAction("Select &Branches",e->setSelMode(RnaDrawController.SelectionType.Branch), getIcon("edit-sel-branch")).setKeyStroke("B");
//    UiAction actionModeBranch = new UiAction("&Branch Slide Mode",e-> setBehaviorMode(RnaDrawController.BehaviorMode.BRANCH), getIcon("mode-branch")).setKeyStroke("B");
//    UiAction actionModeSingle = new UiAction("&Nucleotide Mode",e-> setBehaviorMode(RnaDrawController.BehaviorMode.NUC), getIcon("mode-single")).setKeyStroke("N");

    // View
    UiAction actionZoomAuto = new UiAction("&Auto-Fit", this::autoScale, getIcon("zoom-fit")).setKeyStroke("*V");
    UiAction actionZoomIn = new UiAction("&Zoom In", e->panel.zoom(1), getIcon("zoom-in")).setKeyStroke("*EQUALS").setKeyStroke("*UP");
    UiAction actionZoomOut = new UiAction("Zoom &Out", e->panel.zoom(-1), getIcon("zoom-out")).setKeyStroke("*MINUS").setKeyStroke("*DOWN");
    UiAction actionZoomReset = new UiAction("&Reset Zoom", e-> { panel.setScale(1); panel.resetOffset(); }, getIcon("zoom-reset")).setKeyStroke("*0");
    UiAction actionChooseStructure = new UiAction("&Choose Structure", this::showStructureMenu, getIcon("structure-choose"));
    UiAction actionNextStructure = new UiAction("&Next Structure", this::showNextStructure, getIcon("structure-next")).setKeyStroke("*RIGHT");
    UiAction actionPrevStructure = new UiAction("&Previous Structure", this::showPrevStructure, getIcon("structure-prev")).setKeyStroke("*LEFT");
    UiAction actionShowEnergy = new UiAction("&Calculate Energy", this::showCalcEnergy).setKeyStroke("*E");

    // Structure
    //UiAction actionIdentifyPseudoKnots2 = new UiAction("&Identify Pseudo-Knot Bonds (c++)", this::findPseudoKnotsBackend); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionIdentifyPseudoKnots = new UiAction("&Identify Pseudo-Knot Bonds", this::findPseudoKnots); //, getIcon("zoom-fit")).setKeyStroke("*V")

    UiAction actionSetBondTypePseduo = new UiAction("&Pseudo-Knot Bond", ()->setSelectedBondType(BondType.Pseudo), getIcon("edit-bond-pseudo")); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionSetBondTypeNormal = new UiAction("&Normal Basepair", ()->setSelectedBondType(BondType.Default), getIcon("edit-bond-normal")); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionSetBondTypeForced = new UiAction("&Forced pair (Folding Contstraint)", ()->setSelectedBondType(BondType.Forced), getIcon("edit-bond-forced")); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionSetBondTypeBanned = new UiAction("Forbi&dden pair (Folding Contstraint)", ()->setSelectedBondType(BondType.Prohibited), getIcon("edit-bond-banned")); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionSetBondTypeSpecial = new UiAction("Special pair (alternate color)", ()->setSelectedBondType(BondType.Special), getIcon("edit-bond-special")); //, getIcon("zoom-fit")).setKeyStroke("*V")

    UiAction actionBreakBonds = new UiAction("Remove Basepairs (on selected bases)", controller::breakSelectedBonds, getIcon("edit-rem-bonds")).setKeyStroke("DELETE"); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionSplitStrand = new UiAction("&Split Strand (after each selected base)", controller::splitStrandsAtSelection, getIcon("edit-split-strand")); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionJoinStrands = new UiAction("&Join Selected Strands", controller::joinSelectionStrands, getIcon("edit-join-strand")); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionRemoveNucs = new UiAction("Delete Selected Bases", controller::deleteSelectedBases, getIcon("edit-rem-nucs")); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionInsertNucs = new UiAction("&Insert Bases", this::insertBases, getIcon("edit-insert-nucs")); //, getIcon("zoom-fit")).setKeyStroke("*V")

    UiAction actionEditNucs = new UiAction("Edit Selected Bases (change base type)", this::editSelectedBases, getIcon("edit-base")); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionFold = new UiAction("&Fold Sequence", this::foldSequence, getIcon("edit-fold")).setKeyStroke("*+F");

    //UiAction actionZoomIn = new UiAction("&Zoom In", e->panel.zoom(1), getIcon("zoom-in")).setKeyStroke("*EQUALS").setKeyStroke("*UP");

    // Format
    UiAction actionFormatColor = new UiAction("&Color-Annotate Bases", this::showColorizeDialog, getIcon("format-color")).setKeyStroke("*F"); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionFormatRemoveColor = new UiAction("Remove Color A&nnotations from Selected Bases", e->controller.removeColors(false), getIcon("format-remove-color")); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionFormatRemoveAllColor = new UiAction("Remove &ALL Color Annotations", e->controller.removeColors(true), getIcon("format-remove-color")); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionColorHelices = new UiAction("Color &Helices", this::findHelices); //, getIcon("zoom-fit")).setKeyStroke("*V")
    UiAction actionRedraw = new UiAction("&Redraw", this::redrawScene, getIcon("edit-redraw")).setKeyStroke("*R");
    UiAction actionRedrawCirc = new UiAction("Redraw C&ircular", this::redrawCircular, getIcon("edit-redraw-circ")).setKeyStroke("*+R");
    UiAction actionRedrawLinear = new UiAction("Redraw &Linear", this::redrawLinear, getIcon("edit-redraw-linear")).setKeyStroke("*+L");
    UiAction actionFlipScene = new UiAction("&Flip Scene", e->controller.flip(true, false), getIcon("layout-flip-scene")).setDesc("Flipping a scene is equivalent to converting between Clockwise and Counter-Clockwise layouts.");
    UiAction actionLocalRedraw = new UiAction("&Repair Selected Helices and Closed-Loops", e->controller.redrawSelection()).setKeyStroke("R").setDesc(
            "Every helix in which at least one base is selected will be \"repaired\" by lining up its bases.\n"
           +"Every closed loop in which at least one base is selected will be arranged into an optimal-sized loop with branches spaced according to intervening nucleotides.");
    UiAction actionFormatRotateNumber = new UiAction("Rotate Number Label", controller::rotateNumber).setKeyStroke("*K");

    UiAction testAction  = new UiAction("Test Operation", this::test).setKeyStroke("*F12");


    public DrawWindow() {
        super("Untitled Drawing", true);
        setSize(600,500);
        setLocation(0,0);

        panel = new DrawPanel();
        drawSettings = Program.getInstance().settings().getDrawSettings();
        layout = new NucLayout(drawSettings);

        controller.setSettings(drawSettings);
        controller.setCanvas(panel);
        panel.setRenderer(controller);
        panel.setPreferredSize(new Dimension(200, 200));

//        addKeyBinding(KeyEvent.VK_F5, ActionHelper.COMMAND_MODIFIER_MASK, "refresh");
//        addAction("refresh", e->{ System.out.println("repaint"); repaint(); updateFileUI();} );

        //addKeyBinding('Q', 0, "info");
        //addAction("info", this::showInfo);
        controller.programSettingsChanged(Program.getInstance().settings());
        createUI();
        createMenus();
        createToolButtons();

//        panel.addMouseListener(new MouseAdapter() {
//            /**
//             * {@inheritDoc}
//             *
//             * @param e
//             */
//            @Override
//            public void mouseClicked(final MouseEvent e) {
//                super.mouseClicked(e);
//                if (e.getButton() == MouseEvent.BUTTON3) {
//                    switch (selType) {
//                        case Individual:
//                            setSelType(RnaDrawController.SelectionType.HelixOrLoop);
//                            break;
//                        case HelixOrLoop:
//                            setSelType(RnaDrawController.SelectionType.Branch);
//                            break;
//                        case Branch:
//                            setSelType(RnaDrawController.SelectionType.Individual);
//                            break;
//                    }
//                }
//            }
//        });

        controller.sceneUpdated.add(this::onSceneUpdate);
        controller.historyEvent.add(this::onHistoryUpdate);
        panel.viewChanged.add(this::updateFileUI);

        AdjustmentListener updateView = e->panel.updateView(); // listen for scroll events.
        scrollPane.getVerticalScrollBar().addAdjustmentListener(updateView);
        scrollPane.getHorizontalScrollBar().addAdjustmentListener(updateView);

        panel.setFocusable(true);
        //setBehaviorMode(RnaDrawController.BRANCH_BEHAVIOR_MODE);
        setSelMode(RnaDrawController.SelectionType.Individual);
        loadSettings();
        updateFileUI();
        resetTabTraversalKeys(this);
        resetTabTraversalKeys(panel);
        resetTabTraversalKeys(scrollPane);
    }

    @Override
    public void programSettingsUpdated() {
        controller.programSettingsChanged(Program.getInstance().settings());
    }

    private void createUI() {
        pnlTaskProgress = new JPanel();
        pnlTaskProgress.setLayout(new GridBagLayout());
        lblTaskStatus = new JLabel("Task Status"); lblTaskStatus.setFont(lblTaskStatus.getFont().deriveFont(14f));
        prgTaskProgress = new JProgressBar(JProgressBar.HORIZONTAL);
        lblTaskProgress = new JLabel("0%"); lblTaskProgress.setFont(lblTaskStatus.getFont());
        lblTaskProgress.setHorizontalTextPosition(SwingConstants.CENTER);

        pnlTaskProgress.setVisible(false);

        GridBagConstraints gbc = new GridBagConstraints();

        pnlTaskProgress.add(lblTaskStatus, gbc);
        gbc.insets = new Insets(4,3,4,3);
        gbc.weightx = 1;
        gbc.gridx = 1;
        gbc.fill = GridBagConstraints.BOTH;
        pnlTaskProgress.add(prgTaskProgress, gbc);
        gbc.weightx = 0;
        gbc.gridx = 2;
        gbc.weightx = 0.25;
        gbc.fill = GridBagConstraints.BOTH;
        pnlTaskProgress.add(lblTaskProgress, gbc);
        pnlTaskProgress.setBorder(new BevelBorder(BevelBorder.LOWERED));

        scrollPane = new JScrollPane(panel);
        scrollPane.getVerticalScrollBar().setUnitIncrement(16);
        scrollPane.getHorizontalScrollBar().setUnitIncrement(16);

        add(pnlTaskProgress, BorderLayout.NORTH);
        add(scrollPane, BorderLayout.CENTER);
    }

    private void loadSettings() {
        //Settings s = Program.getInstance().settings();
    }

    private void onHistoryUpdate(final HistoryUpdateEvent event) {
        updateFileUI();
    }

    private void onSceneUpdate(final SceneUpdateEvent event) {
        if (event.scene == currentScene)
            switch (event.type) {
                case Selection:
                    updateSelectionUI();
                    break;
            }
    }
    private void updateSelectionUI() {
        Nuc[] sel = controller.getSelected();
        boolean any = sel.length != 0;
        boolean hasBonds = false;
        if (any)
            for(Nuc n : sel)
                if (n.isPaired()) {
                    hasBonds = true;
                    break;
                }

        actionBreakBonds.setEnabled(hasBonds);
        formatBondType.setEnabled(hasBonds);
        actionJoinStrands.setEnabled(sel.length>1);
        actionSplitStrand.setEnabled(any);
        actionRemoveNucs.setEnabled(any);
        actionInsertNucs.setEnabled(any);
        actionFormatRotateNumber.setEnabled(any);
    }

    private void createToolButtons() {
//        Collections.addAll(menus, jm, new DrawingFileMenu(this), new DrawingEditMenu(this), new DrawingViewMenu(this), new DrawingFormatMenu(this));
        toolButtons.setMergePos(MainFrame.FileMenuSection);
        toolButtons.addSeparator();
        toolButtons.add(actionSave);

        toolButtons.setMergePos(MainFrame.EditMenuSection);
        toolButtons.addSeparator();
        toolButtons.add(actionEditUndo);
        toolButtons.add(actionEditRedo);
        toolButtons.addSeparator();
        toolButtons.add(actionRedraw);
        toolButtons.add(actionRedrawCirc);
        toolButtons.add(actionRedrawLinear);
        toolButtons.addSeparator();
        toolButtons.add(actionFold);
        toolButtons.addSeparator();
//        ButtonGroup grpMode = new ButtonGroup();
//        grpMode.add(btnModeSingle = toolButtons.addToggle(actionModeSingle));
//        grpMode.add(btnModeBranch = toolButtons.addToggle(actionModeBranch));
        ButtonGroup grpSelMode = new ButtonGroup();
        grpSelMode.add(btnSelModeNuc = toolButtons.addToggle(actionSelectNuc));
        grpSelMode.add(btnSelModeHelix = toolButtons.addToggle(actionSelectHelix));
        grpSelMode.add(btnSelModeBranch = toolButtons.addToggle(actionSelectBranch));

        toolButtons.addSeparator();
        toolButtons.add(actionSelectExpand);

//        toolButtons.add(actionSelectHelix);
//        toolButtons.add(actionSelectLoop);
//        toolButtons.add(actionSelectBranch);

        toolButtons.setMergePos(MainFrame.ViewMenuSection);
        toolButtons.addSeparator();
        toolButtons.add(actionZoomAuto);
        toolButtons.add(actionZoomIn);
        toolButtons.add(actionZoomOut);
        toolButtons.add(actionZoomReset);
        toolButtons.addSeparator();
        toolButtons.add(actionPrevStructure);
        structureMenuToolbarButton = toolButtons.add(actionChooseStructure);
        toolButtons.add(actionNextStructure);
    }

    private void createMenus() {
        // File Menu
        MergeMenu m = new MergeMenu("File");
        m.setSubItemMergePos(FileMenu.SAVE_FILE_POS);
        m.addSeparator();
        m.add(actionSave);
        m.add(actionSaveAs);
        m.add(actionExport);
        m.setSubItemMergePos(FileMenu.CLOSE_FILE_POS);
        m.addSeparator();
        m.add(actionClose);
        menus.add(m);

        // Edit Menu
        m = new MergeMenu("Edit");
        m.add(actionEditUndo);
        m.add(actionEditRedo);
        m.addSeparator();
        m.addCheck(actionSelectNuc);
        m.addCheck(actionSelectHelix);
        m.addCheck(actionSelectBranch);
        m.addSeparator();
//        m.addCheck(actionModeSingle);
//        m.addCheck(actionModeBranch);
        m.add(actionSelectExpand);
        m.add(actionSelectAll);
        menus.add(m);

        // View menu
        m = new MergeMenu("&View", MainFrame.ViewMenuSection);
        //m.setSubItemMergePos(FileMenu.CLOSE_FILE_POS);
        m.add(actionZoomAuto);
        m.add(actionZoomIn);
        m.add(actionZoomOut);
        m.add(actionZoomReset);
        m.addSeparator();
        m.add(structureMenu);
        m.add(actionNextStructure);
        m.add(actionPrevStructure);
        m.addSeparator();
        m.add(actionShowEnergy);
        menus.add(m);

        // Structure menu
        m = new MergeMenu("&Structure", MainFrame.FormatMenuSection);
        formatBondType = new MergeMenu("Set Basepair Type");
        formatBondType.add(actionSetBondTypeNormal);
        formatBondType.add(actionSetBondTypePseduo);
        formatBondType.add(actionSetBondTypeForced);
        formatBondType.add(actionSetBondTypeBanned);
        formatBondType.add(actionSetBondTypeSpecial);
        m.add(formatBondType);
        m.add(actionIdentifyPseudoKnots);
         m.add(actionBreakBonds);
        m.addSeparator();
        m.add(actionSplitStrand);
        m.add(actionJoinStrands);
        m.addSeparator();
        m.add(actionEditNucs);
        m.add(actionInsertNucs);
        m.add(actionRemoveNucs);
        m.addSeparator();
        m.add(actionFold);
        menus.add(m);

        m = new MergeMenu("&Format", MainFrame.FormatMenuSection);
        m.add(actionFormatColor);
        m.add(actionFormatRemoveColor);
        m.add(actionFormatRemoveAllColor);
        m.add(actionColorHelices);
        m.addSeparator();
        m.add(actionRedraw);
        m.add(actionRedrawCirc);
        m.add(actionRedrawLinear);
        m.addSeparator();
        m.add(actionFlipScene);
        m.add(actionLocalRedraw);
        m.addSeparator();
        m.add(actionFormatRotateNumber);
        if (asBool(System.getProperty("debug")))
            m.add(testAction);
        menus.add(m);

        structureMenu.getPopupMenu().addPopupMenuListener(new PopupMenuListener() {
            @Override
            public void popupMenuWillBecomeVisible(final PopupMenuEvent e) {

            }
            @Override
            public void popupMenuWillBecomeInvisible(final PopupMenuEvent e) {
                ((JPopupMenu)e.getSource()).setInvoker(structureMenu);
            }
            @Override
            public void popupMenuCanceled(final PopupMenuEvent e) {

            }
        });
    }


    private void save(FileType type, final boolean forceSaveAs, final boolean export) {
        String path = scenes.getFilePath();
        FileType currentType = scenes.getFileType();
        boolean getPath = forceSaveAs || export || isEmpty(path) || !type.includes(currentType);
        FileType defaultExportType = FileType.SVG;
        Preferences mru = Program.getInstance().mruPaths();
        if (getPath) {
            if (isEmpty(path)) {
                path = scenes.getTitle()==null?"":scenes.getTitle().trim();
                path = path.replaceAll("[\\\\/?:\"']|\\.\\.", "_"); // replaces '\' '/' '?' ':' '"' "'" and ".."  with "_"
                if (path.isEmpty()) path = "untitled";
                path += "." + type.firstExtension();
            }

            if (!type.includes(currentType)) {
                if (!PathTools.getExt(path).isEmpty())
                    path = PathTools.changeExtension(path, null);
            }

            // Choose a default for export
            if (export) {
                String lastExportPath = mru.get("export-drawing", null);
                String lastExportType = mru.get("export-drawing-type", null);

                if (!isEmpty(lastExportType)&&FileType.fromName(lastExportType)!=null)
                    defaultExportType = FileType.fromName(lastExportType);

                String dir = null, ext = defaultExportType.firstExtension();
                if (!isEmpty(lastExportPath)) {
                    PathTools.FileName ex = PathTools.parse(lastExportPath);
                    dir = ex.dir();
                    ext = ex.ext(false);
                }

                if ((dir == null || !PathTools.isDir(dir)) && !isEmpty(scenes.getFilePath()))
                    dir = PathTools.getDir(scenes.getFilePath());

                if (dir == null || !PathTools.isDir(dir))
                    dir = PathTools.getDocumentsDir().toString();

                path = dir + File.separatorChar + PathTools.getBaseName(path) + "." + ext;
            }
            String context = export ? "export" : "drawing";
            String title = (export ? "Export":"Save") + " Drawing";
            path = FileDialog.getSaveName(type.getFilterString(true) + "|*", path, title, context, this, defaultExportType.firstExtension());
            if (path == null)
                return;
            type = FileType.fromFilter(FileDialog.getLastChosenFilter());
            if (type == null) type = FileType.NsdDrawing;
        }

        try {
            if (FileType.VectorImage.includes(type))
                saveSceneAsVectorImage(path, type);
            else if (FileType.Image.includes(type))
                saveSceneAsRastorImage(path, type);
            else
                AppActions.writeFile(scenes, path, type, export, false);

            if (export) {
                mru.put("export-drawing", path);
                mru.put("export-drawing-type", type.name());
            }
            if (!export)
                scenes.setSource(path, type);
            if (FileType.AnyOpenable.includes(type))
                Program.getInstance().getMainFrame().addRecentFile(path, type);
        } catch (Exception ex) {
            Dialogs.showWarning(AppActions.formatFileError("save", type, path, ex));
        }
        updateFileUI();
    }
    private Rectangle renderFullSceneToImage(Graphics2D imageGraphics, boolean opaque) {
        Rectangle bounds = getSceneBounds(false);
        AffineTransform tr = new AffineTransform();
        tr.translate(-bounds.x, -bounds.y); // moves image to draw in top-left corner
        bounds.translate(-bounds.x, -bounds.y); // reset to (0, 0)
        View2D view = new View2D(tr, bounds);
        if (opaque) {
            imageGraphics.setColor(panel.getBackground());
            imageGraphics.fill(bounds);
        }
        controller.render(imageGraphics, view, false);
        return bounds;
    }
    private void saveSceneAsRastorImage(final String path, final FileType type) throws IOException {
        // Previously the code below was used to save the image in the viewport.
        // But feedback from the group suggested that the rastor file should be saved at full resolution (with the full image)
        //Rectangle bounds = scrollPane.getViewport().getViewRect();
        //AffineTransform trToScreen = panel.getView().trToScreen;
//        if (bounds.x != 0 || bounds.y != 0) {
//            trToScreen.translate(-bounds.x, -bounds.y);
//        }
        //bounds.setLocation(0, 0);
        //View2D view = new View2D(trToScreen, bounds);

        String format; boolean opaqueBackground = false;
        switch (type) {
            case PNG: format = "PNG"; break;
            case JPEG: format = "JPG"; opaqueBackground = true; break;
            case GIF: format = "GIF"; break;
            default: throw new IllegalArgumentException("Unsupported export type: " + type);
        }

        Rectangle bounds = getSceneBounds(false);
        BufferedImage img = ImageUtil.createCompatibleImage(bounds.width, bounds.height, !opaqueBackground);
        Graphics2D g = img.createGraphics();
        // g.translate(-bounds.x, -bounds.y);
        //controller.render(g, view, false);
        renderFullSceneToImage(g, opaqueBackground);
        ImageIO.write(img, format, new File(path));
    }
    private void saveSceneAsVectorImage(final String path, final FileType type) throws IOException {
        Graphics2DRecorder g = new Graphics2DRecorder();
        Rectangle bounds = renderFullSceneToImage(g, false);
        PageSize pg = new PageSize(bounds.getWidth(), bounds.getHeight());
        if (type == FileType.SVG)
            new SvgGraphicsExporter(pg, SvgGraphicsExporter.TrueScreenDPI ).write(g.getGraphicsOps(), path);
        else
            throw new IllegalArgumentException("Unsupported export type: " + type);
    }
    private Rectangle getSceneBounds(boolean includeInteractive) {
        Graphics2D g = (Graphics2D) this.panel.getGraphics();
        g.setClip(0,0,0,0);
        Rectangle rc = controller.calcBounds(g, new View2D(), includeInteractive);
        g.dispose();
        return rc;
    }

    public void redrawScene() {
        try {
            layout.redrawRadial(currentScene);
            controller.layoutUpdated(SceneUpdateInfo.Layout.subType("Full Redraw"));
            autoScale();
        } catch (RnaBackendException ex) {
            Dialogs.showWarning("Redraw operation failed: \n" + ex.getMessage());
        }
    }

    private void redrawCircular() {
        layout.redrawCircular(currentScene);
        controller.layoutUpdated(SceneUpdateInfo.Layout.subType("Redraw Circular"));
        autoScale();
    }

    private void redrawLinear() {
        layout.redrawLinear(currentScene);
        controller.layoutUpdated(SceneUpdateInfo.Layout.subType("Redraw Linear"));
        autoScale();
    }

    private boolean isTaskRunning() {
        return _asyncTask!=null&&!_asyncTask.isDone();
    }
    private void monitorAsyncTask(RnaFileIO.AsyncTask task) {
        if (isTaskRunning()) throw new UnsupportedOperationException("Another task is already running.");
        _asyncTask = task;
        if (asyncUpdateTimer==null) {
            asyncUpdateTimer = new Timer(100, t->asyncTimer_elapsed());
        }
        asyncUpdateTimer.start();
    }
    private void asyncTimer_elapsed() {
        if (_asyncTask.isDone()) {
            asyncUpdateTimer.stop();
            pnlTaskProgress.setVisible(false);
        } else {
            pnlTaskProgress.setVisible(true);
            lblTaskStatus.setText(_asyncTask.getStatus());
            int progress = _asyncTask.getProgress();
            if (progress < 0) {
                prgTaskProgress.setIndeterminate(true);
                lblTaskProgress.setText("...");
            } else {
                if (progress>100) progress = 100;
                prgTaskProgress.setIndeterminate(false);
                prgTaskProgress.setValue(progress);
                lblTaskProgress.setText(progress+"%");
            }
        }

    }

    private void insertBases() {
//        JPanel pnl = new JPanel();
//        pnl.setLayout(new GridBagLayout());
//        GridBagConstraints gbc = new GridBagConstraints();
        if (controller.getFocused()==null) {
            Dialogs.showInfo("Please select the base at the position where the new bases should be inserted.");
            return;
        }
        InsertBasesDialog dlg = new InsertBasesDialog();
        int result = dlg.showDialog(this, "Insert Bases");
        if (JOptionPane.OK_OPTION == result)
            controller.insertBases(dlg.getSequence(), null, dlg.isInsertAfter());
    }

    private void foldSequence() {
        if (isTaskRunning()) {
            Dialogs.showInfo("Another task is already running. Please wait for it to complete.");
            return;
        }
        // TODO: ask the user if they want to REPLACE or ADD-TO the existing list of structures.
        RnaFileIO.BackendCalc<RnaSceneGroup> task = RnaFileIO.foldSeq(RnaFileIO.getRNASequence(currentScene), 1);
        task.whenDone(t-> {
            if (t.hadError())
                Dialogs.showWarning(t.getError().getMessage());
            else {
                RnaSceneGroup results = t.getResult();
                // TODO: Possibly handle multiple Fold results. (Currently just the first is used.)
                currentScene.clearBonds();
                currentScene.copyBonds(results.get(0));
                try {
                    layout.redrawRadial(currentScene);
                } catch (RnaBackendException ex) {
                    ex.printStackTrace();
                    layout.redrawCircular(currentScene);
                }
                // scenes.addAll(getCurrentSceneIndex()+1, results);
                buildStructureMenu(); // list of scenes have been updated, so reload the menu.
                controller.structureUpdated(SceneUpdateInfo.Bonds.subType("Re-Fold Sequence"));
                updateFileUI();
                autoScale();
            }
        });
        task.start();
        monitorAsyncTask(task);
    }

    private void setSelType(RnaDrawController.SelectionType selType) {
        this.selType = selType;
        controller.setSelType(selType);
        updateFileUI();
    }

    private int getCurrentSceneIndex() {
        return scenes == null || currentScene == null ? -1 : scenes.indexOf(currentScene);
    }

    private void updateFileUI() {
        // ******** Update TITLE ********
        StringBuilder sb = new StringBuilder();
        int currentIndex = getCurrentSceneIndex();
        if (scenes == null)
            sb.append("Untitled Drawing");
        else {
            sb.append("Drawing - ").append(strOrDefault(RnaFileIO.stripEnergyTitle(scenes.getTitle()), "Untitled"));
            if (currentIndex != -1 && currentIndex < scenes.structureCount())
                sb.append(" [").append(strOrDefault(currentScene.title, "Structure " + (currentIndex+1) + " of " + scenes.structureCount())).append("] ");
        }
//        sb.append(" Select: ").append(Strings.toFriendlyName(selType.name()))
//                .append(" Scale: ").append(panel.getScale())
//                .append(" Offset: ").append(panel.getOffset());

        setTitle(sb.toString());

        // ******** Update Next/Prev structure buttons and structure menu ********
        if (currentIndex == -1 || currentIndex >= structureList.getButtonCount())
            structureList.clearSelection();
        else
            structureList.select(currentIndex);
        actionChooseStructure.setEnabled(scenes != null && scenes.structureCount() > 1);
        actionNextStructure.setEnabled(scenes != null && currentIndex < scenes.structureCount() - 1);
        actionPrevStructure.setEnabled(scenes != null && currentIndex > 0);

        actionEditUndo.setEnabled(currentScene != null && controller.canUndo());
        actionEditRedo.setEnabled(currentScene != null && controller.canRedo());
    }

    private String strOrDefault(String value, String defaultIfEmpty) {
        if (isWhiteSpace(value))
            return defaultIfEmpty;
        return value;
    }

    public void loadScene(final RnaScene scene) {
        handleView(currentScene, true);
        currentScene = scene;
        // structureIndex = scene.getStructureIndex();
        controller.setScene(scene);
        handleView(currentScene, false);
        updateFileUI();
        updateSelectionUI();
    }
    private void handleView(final RnaScene scene, boolean save) {
        if (scene == null) return;
        int pos = scenes.indexOf(scene);
        if (pos == -1) return;
        SceneInfo si = sceneInfo[pos];
        if (save) {
            if (si == null)
                sceneInfo[pos] = si = new SceneInfo();
            si.scale = panel.getScale();
            si.offset = panel.getOffset();
            si.scrollPos = scrollPane.getViewport().getViewPosition();
        } else {
            if (si == null) {
                panel.setView(1, new Point2D.Float());
                scrollPane.getViewport().setViewPosition(new Point());
                SwingUtilities.invokeLater(() -> {
                    if (currentScene == scene) autoScale();
                });
            } else {
                panel.setView(si.scale, si.offset);
                final Point viewPos = si.scrollPos;
                SwingUtilities.invokeLater(() -> {
                    if (currentScene == scene)
                        scrollPane.getViewport().setViewPosition(viewPos);
                });
            }
        }
    }

    public void autoScale() {
        panel.autoScale(scrollPane.getViewport().getSize());
        updateFileUI();
    }

    @Override
    public Collection<? extends JMenu> getMenus() { return menus; }

    @Override
    public Collection<? extends Component> getToolbarButtons() { return toolButtons; }

    private void showStructureMenu() {
        JPopupMenu m = structureMenu.getPopupMenu();
        System.out.println("before popup: " + structureMenu.getPopupMenu());
        structureMenu.getPopupMenu().show(structureMenuToolbarButton, 0, structureMenuToolbarButton.getHeight());
        System.out.println("after popup same: " + (structureMenu.getPopupMenu() == m));
        System.out.println("after popup items: " + (structureMenu.getPopupMenu().getComponentCount()));
    }
    private void showNextStructure() { loadStructure(getCurrentSceneIndex() + 1); }
    private void showPrevStructure() { loadStructure(getCurrentSceneIndex() - 1); }
    private void loadStructure(final int structureIndex) {
        if (structureIndex == getCurrentSceneIndex())
            return;
        if (structureIndex > -1 && structureIndex < scenes.structureCount())
            loadScene(scenes.get(structureIndex));
        else
            updateFileUI();
    }
    private void buildStructureMenu() {
        structureMenu.removeAll();
        structureList.removeAll();
        for(int i = 0; i < scenes.structureCount(); i++) {
            JMenuItem mi = new JRadioButtonMenuItem((i+1) + ". " + strOrDefault(scenes.get(i).title, "Structure " + i));
            mi.setActionCommand(Integer.toString(i+1));
            mi.addActionListener(e->loadStructure(Integer.parseInt(e.getActionCommand())-1));
            structureMenu.add(mi);
            structureList.add(mi);
        }
    }


    private void findPseudoKnotsBackend() {
        StopWatch sw = new StopWatch(true);
        try {
            RnaFileIO.identifyPseudoKnotsBackend(currentScene);
            controller.structureUpdated(SceneUpdateInfo.BondType.subType("Identified PseudoKnots"));
        } catch (Exception ex) {
            Dialogs.showWarning("Identifying PseudoKnot Bonds failed: \n" + ex.getMessage());
        }
        sw.println("findPseudoKnots");
    }
    private void findHelices() {
        Color[] colors = new Color[] {
                Colors.Red,
                Colors.Orange,
                Colors.Yellow,
                Colors.DarkGreen,
                Colors.LimeGreen,
                Colors.DodgerBlue,
                Colors.DarkBlue,
                Colors.Purple,
                Colors.Indigo,
                Colors.SaddleBrown,
                Colors.Gray,
        };
        int num = 0;
        for (Motif.Helix h : currentScene.getHelices()) {
            for (Nuc n : h.getBases()) {
                n.style().fillColor = colors[num % colors.length];
            }
            num++;
        }
        controller.styleUpdated(SceneUpdateInfo.FormatBases);
    }
    private void findPseudoKnots() {
        StopWatch sw = new StopWatch(true);

        try {
            Set<Bond> bonds = Motif.findNonCrossingBonds(currentScene);
            for(Bond b : currentScene.getBonds())
                b.type = bonds.contains(b) ?BondType.Default : BondType.Pseudo;
            controller.structureUpdated(SceneUpdateInfo.BondType.subType("Identified PseudoKnots"));
        } catch (Exception ex) {
            Dialogs.showWarning("Identifying PseudoKnot Bonds failed: \n" + ex.getMessage());
        }
        sw.println("findPseudoKnots2");
    }
    public void setSelectedBondType (final BondType type) {
        controller.setSelectedBondType(type);
    }

    private void showColorizeDialog() {
        Preferences prefs = Program.getInstance().prefs().user().node("ColorizeDialog");
        new ColorizeDialog(this, prefs, false).showDialog();
    }
    public void colorizeScenes(SceneColorizer c, boolean applyToAllScenes, boolean selectedBasesOnly) {
        if (applyToAllScenes) {
            RnaScene current = currentScene;
            for (RnaScene scene : scenes) {
                if (scene == current) continue;
                loadScene(scene);
                controller.colorize(c, !selectedBasesOnly);
            }
            loadScene(current);
        }
        controller.colorize(c, !selectedBasesOnly);
    }
    public void showCalcEnergy() {
        try {
            double energy = RnaFileIO.calculateEnergy(currentScene);
            Map<String, Integer> bpCount = new HashMap<>();
            int bp = 0;
            for (Nuc n : currentScene.allNucs()) {
                if (n.isPaired() && n.indexInScene() < n.getPaired().indexInScene()) {
                    bp++;
                    String key = n.symbol + "-" + n.getPaired().symbol;
                    if (n.symbol.compareTo(n.getPaired().symbol) > 0)
                        key = n.getPaired().symbol + "-" + n.symbol;
                    Integer value = bpCount.get(key);
                    if (value == null) value = 0;
                    bpCount.put(key, value + 1);
                }
            }
            StringBuilder sb = new StringBuilder();
            sb.append("Sequence Length: ").append(currentScene.getNucCount()).append('\n')
                    .append("Basepairs: ").append("\n");
            for(String key : bpCount.keySet())
                    sb.append("\t").append(key).append(":\t").append(bpCount.get(key)).append("\n");
            sb.append("\t").append("Total:\t").append(bp).append("\n");
            sb.append("\nGibbs Free Energy: ").append(String.format("%01.2f", energy)).append(" kcal/mol\n");
            Dialogs.showInfo(sb.toString(), "Structure Energy");
        } catch (Exception ex) {
            Dialogs.showWarning("Unable to calculate energy due to an error: " + ex.toString(), "Energy Calculation Failed");
        }
    }
    public void editSelectedBases() {
        Nuc[] sel = controller.getSelected();
        Arrays.sort(sel);
        if (sel.length == 0)
            Dialogs.showInfo("No bases were selected.");
        if (sel.length < 4) {
            for (int i = 0; i < sel.length; i++) {
                String s = Dialogs.input("Enter the replacement for " + sel[i].toString("$s$N"), "Edit Base", sel[i].symbol);
                if (s == null) return;
                s = s.trim();
                if (s.isEmpty()) {
                    Dialogs.showWarning("You cannot enter an empty symbol. Use the \"Delete Selected Bases\" tool to remove bases.");
                    i--;
                    continue;
                }
                if (s.length()!=1) {
                    Dialogs.showWarning("The base symbol can only be a single character.");
                    i--;
                } else {
                    String prev = sel[i].toString("$s$N");
                    sel[i].symbol = s;
                    controller.structureUpdated(SceneUpdateInfo.SequenceText.subType("Changed Base "+prev+" to "+s));
                }
            }
        } else if (sel.length <= 20) {
            StringBuilder sb = new StringBuilder();
            StringBuilder sbAns = new StringBuilder();
            for (int i = 0; i < sel.length; i++) {
                sb.append(sel[i].toString("$s$N  "));
                sbAns.append(sel[i].symbol).append(' ');
            }
            String s = sbAns.toString(); char[] found = new char[sel.length + 1];
            while(true) {
                s = Dialogs.input("Enter the replacements, in order, for the following bases. You must enter exactly " + sel.length + " bases.\n(Spaces are allowed and will be ignored.)\n\n" + sb.toString(), "Edit Base", s);
                if (s == null) return;
                int pos = 0;

                for (int i = 0; i < s.length(); i++) {
                    char c = s.charAt(i);
                    if (!Character.isWhitespace(c))
                        found[pos++] = c;
                    if (pos > sel.length)
                        break;
                }
                if (pos != sel.length)
                    Dialogs.showWarning("You have entered too "+(pos<sel.length?"few":"many")+" base characters.", "Edit Bases");
                else
                    break;
            }
            for (int i = 0; i < sel.length; i++)
                sel[i].symbol = Character.toString(found[i]);
            controller.structureUpdated(SceneUpdateInfo.SequenceText);
        } else {
            sel = currentScene.allNucs().toArray(new Nuc[currentScene.getNucCount()]);
            StringBuilder sb = new StringBuilder(sel.length + sel.length / 10 + 1);
            char[] found = new char[sel.length + 1];
            for (int i = 0; i < sel.length; i++) {
                sb.append(sel[i].symbol);
                if (i != 0 && i % 10 == 0)
                    sb.append(' ');
            }
            String s = sb.toString();
            while(true) {
                s = Dialogs.input("You have selected too many bases to easily edit individually.\nPlease edit the entire sequence below, making sure to preserve the total number of bases.\nA space has been inserted every 10 bases for convenience.", "Edit Base", s);
                if (s == null) return;
                int pos = 0;
                for (int i = 0; i < s.length(); i++) {
                    char c = s.charAt(i);
                    if (!Character.isWhitespace(c))
                        found[pos++] = c;
                    if (pos > sel.length)
                        break;
                }
                if (pos != sel.length)
                    Dialogs.showWarning("You have entered too "+(pos<sel.length?"few":"many")+" base characters. Required: " + sel.length, "Edit Bases");
                else
                    break;
            }
            for (int i = 0; i < sel.length; i++)
                sel[i].symbol = Character.toString(found[i]);
            controller.structureUpdated(SceneUpdateInfo.SequenceText);
        }
    }

    private void test() {
        for(Bond b : currentScene.getBonds()) {
            if (b.isFirst())
                b.n5.style().fillColor = Color.BLUE;
            if (b.isLast())
                b.n3.style().fillColor = Color.ORANGE;
            if (b.isFirst() || b.isLast()) {
                b.right(null).style().textColor = Color.RED;
                b.left(null).style().textColor = Color.GREEN;
            }
        }
        controller.controlsUpdated();

    }
}
