package ur_rna.StructureEditor;

import ur_rna.StructureEditor.services.UserPrefs;
import ur_rna.Utilities.*;
import ur_rna.Utilities.annotation.Nullable;
import ur_rna.Utilities.swing.Dialogs;
import ur_rna.Utilities.swing.ImageUtil;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.InputStream;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.prefs.Preferences;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Main class for running the structure editor as a standalone drawing program.
 * author: Richard M. Watson
 *
 */
public class Program implements AppInfo.AppInfoProvider {
    public static final int VERSION_MAJOR = 1;
    public static final int VERSION_MINOR = 0;
    public static final String TITLE = "Structure Editor";
    public static final String APPNAME = "RNAeditor"; // a shorter form for the name of the preferences file etc.
    public static final String RESOURCES = "/ur_rna/StructureEditor/resources/";
    public static final Version VERSION = new Version(VERSION_MAJOR, VERSION_MINOR);
    public static final String C_YEAR = "2016 - 2017"; // copyright year(s)
    public static String DrawingFileExtension = "nsd"; // the default extension for Structure Drawing files.
    public static final AppLog log = AppLog.getDefault();

    private static final Program instance = new Program();
    public static Program getInstance() { return instance; }
    public static String getVersion() {
        return VERSION_MAJOR + "." + VERSION_MINOR;
    }

    private UserPrefs userPrefs = new UserPrefs(APPNAME);
    private Settings settings = new Settings();
    private MainFrame mainFrame;

    public UserPrefs prefs() {
        return userPrefs;
    }
    public Settings settings() { return settings; }

    /**
     * The main program entry-point.
     *
     * @param args   The command line arguments
     */
    public static void main( String[] args ) {
        instance.settings.load(instance.userPrefs.user().node("settings"));
        StartupArgs startup = parseCommandArgs(args); //checks for some command-line arguments and sets System properties.
        log.readSystemProperties(); //update values based on command-line arguments.
        OSInfo.applyNativeLookAndFeel();
        OSInfo.useNativeMenus();
        MainFrame frame = instance.mainFrame = new MainFrame();
        frame.setVisible(true);
        if (!loadNativeLib()) {
            System.exit(2);
            return;
        }
        frame.addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(final WindowEvent e) {
                instance.settings.save(instance.userPrefs.user().node("settings"));
                super.windowClosing(e);
            }
        });
        instance.doInitialStartup();
        frame.showDashboard(false);
        for(String file : startup.files)
            AppActions.openFile(file, FileType.AnyOpenable);
    }
    public Preferences mruPaths() {
        return prefs().user().node("MRUPaths");
    }

    private static class StartupArgs {
        public List<String> files = new ArrayList<>();
    }
    private static StartupArgs parseCommandArgs(String[] args) {
        StartupArgs info = new StartupArgs();
        try {
            boolean endOfOptions = false;
            for (int i = 0; i < args.length; i++) {
                String arg = args[i];
                if (arg.startsWith("-")&&!endOfOptions)
                    switch (arg.toLowerCase()) {
                        case "-d": System.setProperty("debug", "true"); break;
                        case "-v": System.setProperty("verbose", args[++i]); break;
                        case "-sfc": System.setProperty("simple-file-chooser", "true"); break; // helps with automated testing.
                        case "--": endOfOptions = true; break;
//                    case "-jlaf": options.useNativeLAF = false; break;
//                    case "-jmnu": options.useNativeMenus = false; break;
//                    case "-prefs": preferencesNode = args[++i]; break;
//                    case "-opt": {// override a user option (stored in preferences)
//                        int pos = args[++i].indexOf('=');
//                        if (pos == -1)
//                            throw new InvalidPreferencesFormatException("The -opt flag requires an argument in the form \"name=value\".");
//                        UserOptions.override(args[i].substring(0, pos), args[i].substring(pos + 1));
//                    }
                        default: log.error("Unknown command-line option: \"%s\".", args[i]); break;
                    }
                else // plain argument -- not a -flag
                    info.files.add(arg);
            }
        } catch (Throwable e) {
            log.error("error parsing command line.", e);
        }
        return info;
    }

    public MainFrame getMainFrame() {
        return mainFrame;
    }

    //Implementations of AppInfoProvider
    public Version getAppVersion() {
        return VERSION;
    }
    public String getAppTitle() {
        return TITLE;
    }
    public String getAppName() {
        return APPNAME;
    }
    public AppInfo.AppInfoProvider getAppInfo() { return instance; }

    public static boolean loadNativeLib() {
        try {
            final String libName = "RNAstructure_GUI";
            JniLoader loader = new JniLoader("RNAstructure",
                    // RNAstructureBackendCalculator::getVersion, REQUIRED_DLL_VERSION,
                    libName + JniLoader.getJvmBits(), libName);
            loader.load();
            AppLog.getDefault().info("Loaded Native DLL: " + loader.getLoadedLib());
            return true;
        } catch (JniLoader.NativeLibLoadError ex) {
            Dialogs.showWarning(ex.getMessage());
            return false;
        }
    }

    private void doInitialStartup() {
//        final String INIT_RUN_KEY = "initial_run_" + getVersion();
//        String val = userPrefs.get(INIT_RUN_KEY);
//        if (!asBool(val)) {
//            userPrefs.put(INIT_RUN_KEY, "1");
//            AppActions.SHOW_ABOUT_WINDOW.invoke(this);
//        }
    }
    public static BufferedImage getImage(String name) {
        if (name.indexOf('.') == -1)
            name += ".png";
        try {
            URL res = getResource("images/"+ name);
            if (res != null) return ImageUtil.toBufferedImage(ImageIO.read(res), true);
            log.error("Error loading image resource: \"" + name + "\". Resource not found.");
        } catch (Exception ex) {
            log.error("Error loading image resource: \"" + name + "\"", ex);
        }
        return null;
    }

    @Nullable
    public static ImageIcon getIcon(String name) {
        Image image = getImage(name);
        return image == null ? null : new ImageIcon(image, name);
    }
    @Nullable
    public static URL getResource(String relativePath) { return Program.class.getResource(RESOURCES+relativePath); }
    @Nullable
    public static File getResourceFile(String relativePath) {
        URL res = getResource(relativePath);
        return res==null?null:new File(res.toString());
    }
    public static InputStream getResourceStream(String relativePath) { return Program.class.getResourceAsStream(RESOURCES+relativePath); }
    public static List<String> getImageResources() {
        List<String> files = new ArrayList<>();
        try {
            File res = getResourceFile("images");
            if (res == null) return Collections.emptyList();
            File[] list = res.listFiles();
            if (list == null) return Collections.emptyList();
            for(File f : list)
                if (f.getName().endsWith(".png"))
                    files.add(f.getName());
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        return files;
    }
    public static String getResourceText(final String messageName) { return getResourceText(messageName, null); }
    public static String getResourceText(final String messageName, String sectionName) {
        String path = null;
        String[] exts = new String[]{".txt", ".htm", ".html", ""};
        for (String ext : exts)
            if (getResource("messages/" + messageName + ext) != null) {
                path = "messages/" + messageName + ext;
                break;
            }
        if (path != null)
            try (InputStream stream = getResourceStream(path)) {
                //try (InputStream stream = new FileInputStream("E:\\home\\jobs\\RNA\\RNAstructure\\alignment-editor\\src\\ur_rna\\StructureEditor\\resources\\messages\\new-drawing-help.txt")) {
                String content = Strings.readAll(stream, StandardCharsets.UTF_8);
                if (sectionName != null && !sectionName.isEmpty()) {
                    Pattern p = Pattern.compile("^<section +" + sectionName + " *>(?:\\r?\\n)?(.*?)(?:\\r?\\n)?^</section>", Pattern.DOTALL | Pattern.MULTILINE);
                    Matcher m = p.matcher(content);
                    content = m.find() ? replaceResUrls(m.group(1)) : "<Failed to find section " + sectionName + ">";
                }
                return content;
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        return "<Failed to open message " + messageName + ">";
    }

    static Pattern pattern = Pattern.compile("\\bres://([^\\s\"']+)");
    public static String replaceResUrls(String input) {
        StringBuffer output = new StringBuffer();
        Matcher matcher = pattern.matcher(input);
        while (matcher.find()) {
            URL res = getResource(matcher.group(1));
            matcher.appendReplacement(output, res == null ? matcher.group() : res.toString());
        }
        matcher.appendTail(output);
        return output.toString();
    }

//    public void invoke(String textCommand) {
//        int pos = textCommand.indexOf(' ');
//        if (pos == -1)
//            invoke(AppActions.valueOf(textCommand));
//        else {
//            try {
//                String command = textCommand.substring(0, pos);
//                NamedArg[] args = parseArgs(textCommand.substring(pos + 1));
//                invoke(AppActions.valueOf(command), args);
//            } catch (SyntaxErrorException ex) {
//                JOptionPane.showMessageDialog(mainFrame, "Error invoking application command: " + ex.toString(), "Application Error", JOptionPane.ERROR_MESSAGE);
//            }
//        }
//    }
//    private NamedArg[] parseArgs(final String text) throws SyntaxErrorException {
//        CommandLineParser p = new CommandLineParser(true, CommandLineParser.ParamAffinity.OPTIONAL, true);
//        CommandLineParser.ParseResults r = p.parse(text);
//        NamedArg[] args = new NamedArg[r.arguments.size()];
//        for (int i = 0; i < args.length; i++) {
//            CommandLineParser.Argument arg = r.arguments.get(i);
//            if (arg instanceof CommandLineParser.FlagArg) {
//                CommandLineParser.FlagArg fa = (CommandLineParser.FlagArg)arg;
//                args[i] = new NamedArg(fa.flag.name, fa.valueText);
//            } else
//                args[i] = new NamedArg(null, arg.fullText);
//        }
//        return args;
//    }
//
//    public void invoke(AppActions command, NamedArg... namedArgs) {
//        Map<String,Object> args = NamedArg.toMap(namedArgs, false);
//        switch (command) {
//            case SHOW_ABOUT_WINDOW:
//                mainFrame.showCentered(new AboutWindow());
//                break;
//            case SHOW_OPEN_SEQ:
//                showOpenSequenceFile();
//                break;
//            case SHOW_OPEN_STRUCTURE:
//                showOpenFile();
//                break;
//            case EXIT_PROGRAM:
//                mainFrame.close();
//                break;
//            case OPEN_RECENT_FILE:
//                openRecentFile((String)args.get("file"), (String)args.get("type"));
//                break;
//            case OPEN_DEMO:
//                createEditorWindow(createDemoScene(Convert.toInt(args.get("type"))));
//                break;
//        }
//    }
//
//
//    /**
//     * Invoked when an action occurs.
//     */
//    @Override
//    public void actionPerformed(final ActionEvent e) {
//        invoke(e.getActionCommand());
//    }
}
