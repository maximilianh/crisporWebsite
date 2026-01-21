package ur_rna.StructureEditor;

import ur_rna.RNAstructure.RnaBackendException;
import ur_rna.StructureEditor.models.RnaSceneGroup;
import ur_rna.StructureEditor.services.drawing.NucLayout;
import ur_rna.StructureEditor.services.fileIO.DrawingFileIO;
import ur_rna.StructureEditor.services.fileIO.RnaFileIO;
import ur_rna.StructureEditor.windows.AboutWindow;
import ur_rna.StructureEditor.windows.DashboardFrame;
import ur_rna.StructureEditor.windows.DrawWindow;
import ur_rna.StructureEditor.windows.NewFileWindow;
import ur_rna.Utilities.*;
import ur_rna.Utilities.swing.Dialogs;
import ur_rna.Utilities.swing.FileDialog;
import ur_rna.Utilities.swing.UiAction;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.net.URLEncoder;

import static ur_rna.Utilities.ObjTools.asBool;

public abstract class AppActions {
    private static Program _app;
    private static MainFrame _window;
    private static AppLog log = AppLog.getDefault();
    protected static Program app() { return _app == null ? _app = Program.getInstance() : _app; }
    protected static MainFrame window() { return _window == null ? _window = app().getMainFrame() : _window; }

    private static final String WEB_URL = "http://rna.urmc.rochester.edu/";
    private static final String HELP_URL = WEB_URL + "GUI/html/StructureEditor.html";
    private static final String FEEDBACK_URL = "http://rna2.urmc.rochester.edu/RNAstructureWeb/feedback.php";
    private static final String FEEDBACK_SECRET = "18qY19wP21fR21uCuL1v";
    private static final String LOCAL_HELP_PAGE = "StructureEditor.html";
    private static final String[] LOCAL_HELP_URLS = new String[] { "manual/html", "manual/GUI/html" };


    public static final UiAction
            SHOW_ABOUT_WINDOW = new UiAction("About", AppActions::showAboutWindow, Program.getIcon("help-about")),
            OPEN_WEBSITE = new UiAction("Mathews Lab &Website", e->browse(WEB_URL)),
            SHOW_ONLINE_HELP = new UiAction("&Online Help", e->browse(HELP_URL)),
            SHOW_LOCAL_HELP = new UiAction("Program &Help", AppActions::browseLocalHelp, Program.getIcon("help-docs")),
            SEND_FEEDBACK = new UiAction("Send &Feedback or Bug Report", AppActions::sendFeedbackFromBrowser, Program.getIcon("help-feedback")),
            NEW_FILE = new UiAction("&New File (enter sequence...)", AppActions::createNewFile, Program.getIcon("file-new")),
            SHOW_OPEN_FILE = new UiAction("&Open File", AppActions::showOpenFile, Program.getIcon("file-open")),
            EXIT_PROGRAM = new UiAction("&Exit", AppActions::exitApp),
            //OPEN_DEMO = new UiAction("&Demo", e-> createEditorWindow(createDemoScene(0))),
            SHOW_DASHBOARD = new UiAction("&Home", e-> window().showDashboard(true)),
            CLEANUP_RECENT_FILES = new UiAction("&Remove missing files", e-> window().cleanRecent(false), Program.getIcon("clear-recent")),
            CLEAR_RECENT_FILES = new UiAction("&Remove ALL recent files", e-> window().cleanRecent(true), Program.getIcon("delete"))
    ;

    static {
        SHOW_OPEN_FILE.setDesc("Browse to open a Drawing, Structure, or Sequence file.").setKeyStroke("*o");
        NEW_FILE.setDesc("Create a new drawing by entering a sequence and optional base pairing information.").setKeyStroke("*N");
        //SHOW_OPEN_SEQ.setDesc("Open a SEQ or FASTA file and fold it into a probable structure.").setKeyStroke("*+o");
        EXIT_PROGRAM.setDesc("Close all open windows and exit the program.").setKeyStroke("*+w");
        SHOW_ONLINE_HELP.setKeyStroke("+F1");
        SHOW_LOCAL_HELP.setKeyStroke("F1");
    }

//    protected static UiAction createAction(String name, Consumer<MainFrame> c) { return new UiAction(e -> c.accept(window()), name, null); }
//    protected AppActions(String name, Runnable r) { super(r, name, null); }
//    protected AppActions(String name, Runnable r, Icon i) { super(r, name, i); }
//    protected AppActions(String name, ActionListener l, Icon i) { super(l, name, i); }

//    public static List<UiAction> getHandlers() {
//        if (handlers == null) {
//            List<UiAction> list = new ArrayList<>();
//            for (Field f : AppActions.class.getFields()) {
//                if (Modifier.isStatic(f.getModifiers()) && UiAction.class.isAssignableFrom(f.getType()))
//                    try {
//                        list.add((UiAction) f.get(null));
//                    } catch (IllegalAccessException ex) {
//                        ex.printStackTrace();
//                    }
//            }
//            handlers = Collections.unmodifiableList(list);
//        }
//        return handlers;
//    }

    public static void exitApp() {
        window().close();
    }

    public static void showAboutWindow() {
        window().showCentered(new AboutWindow());
    }

    private static void browseLocalHelp() {
        for(String dir : LOCAL_HELP_URLS) {
            File page = PathTools.findLocalPath(dir + '/' + LOCAL_HELP_PAGE);
            if (page != null) {
                browse(page.toURI().toString());
                return;
            }
        }
        browse(HELP_URL);
    }
    private static String urlenc(String value) {
        try {
            return URLEncoder.encode(value, "UTF-8");
        } catch (UnsupportedEncodingException ex) {
            return value;
        }
    }

    private static void sendFeedbackFromBrowser() {
        String url = FEEDBACK_URL;
        try {
            // Attempt to get a new one-time-password token
            String token;
            try (InputStream netStream = new java.net.URL(FEEDBACK_URL + "?action=get-token").openStream()) {
                token = Strings.readAll(netStream);
            } catch (IOException ex) {
                ex.printStackTrace();
                token = null;
            }

            StringBuilder sb = new StringBuilder();
            sb.append(FEEDBACK_URL)
                    .append("?program=").append(urlenc(Program.TITLE))
                    .append("&version=").append(urlenc(Program.getVersion()));
            if (token != null)
                sb.append("&passtoken=").append(urlenc(token)).append("&passresult=").append(Strings.md5(FEEDBACK_SECRET+token));

            url = sb.toString();
            Desktop.getDesktop().browse(new URI(url));
        } catch( Exception e ) {
            Dialogs.showWarning( "error showing feedback webpage in browser.\nError: " + e.toString() + "\nURL: " + url);
        }
    }

    public static void browse(String uri) {
        try {
            Desktop.getDesktop().browse(new URI(uri));
        } catch( Exception e ) {
            Dialogs.showWarning( "error showing webpage in browser: " + uri);
        }
    }

    public static void openFile(final String file, final String type) {
        try {
            openFile(file, FileType.fromAny(type, null, null, file, FileType.AnyOpenable));
        } catch (IllegalArgumentException ex) {
            ex.printStackTrace();
            Dialogs.showWarning("Unknown file type: " + type);
        }
    }

    public static DrawWindow createEditorWindow(RnaSceneGroup sceneGroup) {
        DrawWindow f = new DrawWindow();
        f.loadFile(sceneGroup);
        window().show(f);
        f.autoScale();
        f.controller.clearUndo();
        return f;
    }

//    public static void showOpenSequenceFile() {
//        String file = FileDialog.getOpenName(FileType.Sequence.getFilterString(), null, "Open Sequence", "sequence", window());
//        if (file != null) openSequenceFile(file);
//    }

    public static void createNewFile() {
        window().showCentered(new NewFileWindow());
    }
    public static void showOpenFile() {
        String path = FileDialog.getOpenName(FileType.AnyOpenable.getFilterString() + "|*", null, "Open Drawing", "drawing", window());
        if (path == null) return;
        openFile(path, FileType.AnyOpenable);
    }

    public static void openFile(String path, FileType type) {
        try {
            if (type == null || type == FileType.AnyOpenable)
                type = FileType.fromExt(PathTools.getExt(path));
            if (type == null) type = showFileTypePrompt();
            if (type == null) return;

            boolean isExample = path.startsWith(DashboardFrame.EXAMPLE_INDICATOR);
            if (isExample) {
                path = DashboardFrame.extractExample(path);
                if (path == null) {
                    Dialogs.showWarning("The example file was not found or could not be read.");
                    return;
                }
            }
            RnaSceneGroup file = readSceneFile(path, type);
            if (file != null) {
                if (isExample) // remove the path from example files, and do not add them to the list of recent files.
                    file.setSource(new File(path).getName(), type);
                else
                    window().addRecentFile(path, type);
                createEditorWindow(file);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            Dialogs.showWarning(formatFileError("open", type, path, ex));
        }
    }

    public static FileType showFileTypePrompt() {
        Object[] selectionValues = {FileType.NsdDrawing, FileType.CteDrawing, FileType.CT, FileType.DotBracket, FileType.Sequence};
        return (FileType)JOptionPane.showInputDialog(Program.getInstance().getMainFrame(), "Please indicate the file type:", "Select File Type", JOptionPane.QUESTION_MESSAGE, null, selectionValues, FileType.NsdDrawing);
    }

    public static RnaSceneGroup readSceneFile(final String path, final FileType type)
            throws RnaBackendException, SyntaxErrorException, IOException {
        //DEBUG_TIMING: System.out.println("Reading file @"+System.currentTimeMillis()); long time = System.currentTimeMillis();
        System.out.printf("Reading file %s (%s)\n", path, type);
        RnaSceneGroup group;
        if (type == null)
            throw new IllegalArgumentException("File type is required.");
        switch (type) {
            case Drawing:
            case NsdDrawing:
                group = DrawingFileIO.readNsdDrawingFile(path);
                break;
            case CteDrawing:
                group = DrawingFileIO.readCteDrawingFile(path);
                break;
            case Structure:
            case CT:
                //DEBUG_TIMING: StopWatch w = new StopWatch(true);
                group = DrawingFileIO.readCteDrawingFile(path);
                //DEBUG_TIMING: w.println("Read CTE").restart();
                createLayout().redrawRadial(group);
                //DEBUG_TIMING: w.println("Redraw").restart();
                break;
            case DotBracket:
            case FoldingSav:
            case PartitionSav:
            case Sequence:
            case Seq:
            case Fasta:
            case SeqText:
                group = RnaFileIO.readFileAsScene(path, type);
                createLayout().redrawRadial(group);
                break;
            default:
                throw new RuntimeException(String.format("%s cannot open %s as drawings.", Program.TITLE, type.getDescription()));
        }
        group.setSource(path, type);
        return group;
    }

    private static NucLayout createLayout() {
        return new NucLayout(app().settings().getDrawSettings());
    }


    public static boolean writeFile(RnaSceneGroup scenes, final String path, final FileType type, boolean asCopy, boolean append)
            throws IOException, RnaBackendException, FormatterException {

        FileType originalType = scenes.getFileType();
        boolean canUpdatePath = originalType == FileType.AnyOpenable;
        if (type == null)
            throw new IllegalArgumentException("File type is required.");

        // canUpdatePath is a boolean that indicates whether the path and filetype of the
        // Group should be updated by this save operation.
        // i.e. When working in Excel, .xls file and then save it as a .xlsm file, the path and type change,
        //      but if you save it as a CSV file, the path and type remain as the xls file, because CSV
        //      is in a lower tier of data storage. Similarly a Drawing is in a higher tier than a structure, which is
        //      still higher than a sequence. So when a drawing is saved as a sequence, the original drawing filename
        //      and type should not be modified. Conversely, when a sequence or structure are modified and saved as a
        //      drawing, it is fine to change the path and type to the drawing file.
        //
        //      If the user has indicated they are saving a COPY, then the path and type are never updated.
        //      In contrast, if the user has created a new file, then any save operation will update the path and type.

        switch (type) {
            case NsdDrawing:
                DrawingFileIO.writeNsdDrawingFile(scenes, path, append);
                canUpdatePath |= originalType.isSubTypeOf(FileType.NsdDrawing, FileType.Structure, FileType.Sequence);
                break;
            case CteDrawing:
                DrawingFileIO.writeCteDrawingFile(scenes, path, append);
                canUpdatePath |= originalType.isSubTypeOf(FileType.CteDrawing, FileType.Structure, FileType.Sequence);
                break;
            case DotBracket:
            case CT:
                new RnaFileIO().writeStructureFile(scenes, path, type, append);
                canUpdatePath |= originalType.isSubTypeOf(FileType.Structure, FileType.Sequence);
                break;
            case Fasta:
            case Seq:
            case SeqText:
                new RnaFileIO().writeSequenceFile(scenes.get(0), path, type, append);
                canUpdatePath |= FileType.Sequence.includes(originalType);
                break;
            default:
                throw new IllegalArgumentException("Programming Error: writeFile cannot handle files of type: " + type);
        }
//            if (err != 0)
//                throw new RnaBackendException(fmt("Failed to save %s File \"%s\": %s\nFull Path: \"%s\"\nInternal Error Number: %s", type.name(), outFile.getName(), rna.GetFullErrorMessage(), getCanonicalPath(outFile), err), err, path);
        if (canUpdatePath && !asCopy)
            scenes.setSource(path, type);
        return true;
    }

    public static String formatFileError(String operation, FileType type, String path, Exception ex) {
        File file = new File(path);
        String msg = String.format("Failed to %s %s \"%s\":\n\n%s\n\nFull Path: %s",
                operation,
                type==null?"File":type.getDescription(),
                file.getName(),
                ExceptionHelper.getMessage(ex).trim().replace("Error opening file: ","").replace("Error saving file: ",""),
//                ex.getMessage()==null ?
//                        ex.getClass().getSimpleName() :
//                        ex.toString().trim().replace("Error opening file: ","").replace("Error saving file: ",""),
                PathTools.getCanonicalPath(file));

        if (asBool(System.getProperty("debug")))
            msg += "\n\n\n\nError Stack Trace (for debugging):\n" + ExceptionHelper.getStackTrace(ex);

        // fmt("Failed to load %s File \"%s\": %s\nFull Path: \"%s\"\nInternal Error Number: %s", type.name(), tmp.getName(), message, getCanonicalPath(tmp), err)
        return msg;
    }

//    private static RnaSceneGroup createDemoScene(int demoModel) {
//        RnaScene scene  = new RnaScene();
//        Strand s = new Strand();
//        scene.strands.add(s);
//
//        int nucRad = 16;
//
//        switch (demoModel) {
//            case 0:{
//                int count = 10;
//                int circ = (int)(nucRad *  1.4 * (count + 1));
//                int radius = (int)(circ / 2 / Math.PI);
//                int margin = 15;
//                double angle = 2 * Math.PI / (count + 1);
//
//                double x = radius * Math.cos(angle + Math.PI / 2) + nucRad*1.4*0.5;
//                double y = radius * Math.sin(angle + Math.PI / 2);
//
//                int num = 0;
//                for (int i = 0; i < 6; i++)
//                    s.add(new Nuc("" + ++num, margin + (int)x, margin + (int) (y + (6 -i) * nucRad * 1.4)));
//
//                for (int i = 1; i < count+1; i++) {
//                    double theta = angle * i + Math.PI / 2;
//                    s.add(new Nuc("" + ++num, margin + (int) (radius * Math.cos(theta)), margin + (int) (radius * Math.sin(theta))));
//                }
//                for (int i = 0; i < 6; i++)
//                    s.add(new Nuc("" + ++num, margin + (int)(x+nucRad*1.4), margin + (int) (y + (i +1) * nucRad * 1.4)));
//
//                scene.reIndex();
//
//                scene.bonds.add(new Bond(s.get(5), s.get(s.size() - 5 - 1)));
//
////                for (int i = 0; i < 6; i++)
////                    scene.bonds.add(new Bond(s.get(i), s.get(s.size() - i - 1)));
//
//            }
//            break;
//            case 1: {
//                int count = 180;
//                int circ = (int)(nucRad * 1.2 * count);
//                int radius = (int)(circ / 2 / Math.PI);
//                int margin = 5;
//                double angle = 2 * Math.PI / count;
//                for (int i = 0; i < count; i++) {
//                    double theta = angle * i;
//                    s.add(new Nuc("" + i, margin + (int) (radius * Math.cos(theta)), margin + (int) (radius * Math.sin(theta))));
//                }
//                s.add(new Nuc("0,0", 0, 0));
//            }
//            break;
//        }
//        scene.reIndex();
//        return RnaSceneGroup.from(scene);
//    }
}
