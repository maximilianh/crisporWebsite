/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package ur_rna.RNAstructureUI;

import ur_rna.RNAstructure.RNAstructureVersion;
import ur_rna.RNAstructureUI.ui.Dialogs;
import ur_rna.RNAstructureUI.ui.StandardFileChooser;
import ur_rna.RNAstructureUI.utilities.*;
import ur_rna.Utilities.*;
import ur_rna.Utilities.annotation.ApplicationInfo;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.prefs.BackingStoreException;
import java.util.prefs.InvalidPreferencesFormatException;
import java.util.prefs.Preferences;

import static ur_rna.Utilities.Strings.asBool;
import static ur_rna.Utilities.Strings.isEmpty;

/**
 * A class that starts the RNAstructure interface.
 *
 * @author Jessica S. Reuter
 * @author Richard M. Watson
 */
@ApplicationInfo(name = "RNAstructure", title = "RNAstructure GUI", version = RNAstructureVersion.VERSION)
public class RNAstructure
	implements Serializable {
	private static final long serialVersionUID = 20120802;
	public static final AppLog log = AppLog.getDefault();
	public static MRUFileStorage MRUFiles = new MRUFileStorage();
	public static UserOptions options = new UserOptions();
	public static String preferencesNode = "default";

	// ### Important: The following lines should be edited prior to each RNAstructure Release ###
	public static final String VERSION = RNAstructureVersion.VERSION;
	public static final String REQUIRED_DLL_VERSION = RNAstructureVersion.VERSION;
	public static final String RELEASE_DATE = RNAstructureVersion.RELEASE_DATE;

	/**
	 * The main method.
	 *
	 * @param args   The command line arguments
	 */
	public static void main( String[] args ) {
		try {
			parseCommandArgs(args); //checks for some command-line arguments and sets System properties.
			log.readSystemProperties(); //update values based on command-line arguments.

			// If the OS is a type of Mac, certain things need to be set in
			// order to make the native Java more "Mac-like."
			if (options.useNativeLAF)
				OSInfo.applyNativeLookAndFeel();

			if (options.useNativeMenus)
				OSInfo.useNativeMenus();

			log.debug("OS: %s Arch: %s Bits: %s%n", System.getProperty("os.name"), System.getProperty("os.arch"), System.getProperty("sun.arch.data.model"));
			log.debug("user.home: " + System.getProperty("user.home"));
			log.debug("user.dir: " + System.getProperty("user.dir"));
			log.debug("CurDir: " + new File(".").getCanonicalPath());
			log.debug("java.library.path: " + System.getProperty("java.library.path"));
			log.debug("java.class.path: " + System.getProperty("java.class.path"));
			if (log.isTraceEnabled())
				System.getProperties().list(log.getTrStream());
			log.debug("Headless: " + java.awt.GraphicsEnvironment.isHeadless());

			try {
				final String libName = "RNAstructure_GUI";
				JniLoader loader = new JniLoader("RNAstructure",
						RNAstructureBackendCalculator::getVersion, REQUIRED_DLL_VERSION,
						libName + JniLoader.getJvmBits(), libName);
				loader.load();
				// This line was commented out avoid an illegal reflective access operation that was prohibited in Java 10. 
				// The function loader.getLoadedLib should be modified to fix this in the long term.
				// log.debug("Loaded Native DLL: " + loader.getLoadedLib());
			} catch (JniLoader.NativeLibLoadError ex) {
				showWarning(ex.getMessage());
                System.exit(2);
                return;
			}

			options.load(getPrefs().node("options"));
			Resources.verify();
			verifyDataPath();
			StandardFileChooser.loadRecentPaths(getPrefs().node(StandardFileChooser.PREFS));
			setDefaultFileOpenDir();
			MRUFiles.loadFromStorage();

			// The splash screen window is loaded by a native Java library
			// before the JVM loads (based on the corresponding entry
			// in the Manifest file). It is closed automatically as soon
			// as the first window is displayed by Swing/AWT.
			// So there is no need to manipulate the splash screen directly.
			// If running directly, without compiling to JAR you can specify
			// the splash image by passing the flag  -splash:<path_to_image.gif>
			// on the java command line.

			// Get a reference to the splash screen if it has been shown.
			// getSplashScreen will return null if the splash image is not specified in the
			//   manifest or passed in on the command-line.
			// This would happen for example, if the application is launched
			//   in an IDE directly, as opposed to via a JAR file.
			java.awt.SplashScreen splash = java.awt.SplashScreen.getSplashScreen();
			//Note: Splash screen fails on mac app-bundle: https://bugs.openjdk.java.net/browse/JDK-8090606
			if (splash != null && !silentMode())
				Thread.sleep(2000); //give the splash screen 2 seconds to display before launching the main window.

			// Start the GUI.
			new AppMainFrame();
			//JFrame frm = new AppMainFrame();
//			frm.getRootPane().addPropertyChangeListener(new PropertyChangeListener() {
//				@Override
//				public void propertyChange(final PropertyChangeEvent evt) {
//					log.debug("PropertyChangeEvent: " + evt.toString());
//				}
//			});
		} catch( Exception e ) {
			log.error("error loading RNAstructure.", e);
			try {
				java.io.FileWriter w = new java.io.FileWriter("RNAstructure_error.log");
				w.write(DateHelper.getFormattedDate(DateHelper.ISO8601) + "\n");
				e.printStackTrace(new java.io.PrintWriter(w));
				w.flush();
				w.close();
			} catch (Throwable ex) {
				showWarning("Could not write error log:\n" + AppLog.getErrorInfo(ex));
			}
			System.exit( 10 );
		}
	}

	public static Preferences getPrefs() {
		return Preferences.userNodeForPackage(RNAstructure.class).node(preferencesNode);
	}
	public static String getPref(String name) { return getPref(name, null); }
	public static String getPref(String name, String defaultValue) {
//		try {
			Preferences p = getPrefs();
			return p.get(name, defaultValue);
			//if (p.nodeExists(name))
//		} catch (BackingStoreException ex) {
//			// return default below.
//			ex.printStackTrace();
//		}
//		return defaultValue;
	}
	public static void setPref(String name, String value) {
		Preferences p = getPrefs();
		p.put(name, value);
		try {
			p.flush();
		} catch (BackingStoreException ex) {
			ex.printStackTrace();
		}
	}

	private static void setDefaultFileOpenDir() {
		File examples = getExamplesDir();
		if (examples != null)
			StandardFileChooser.setDefaultDir(examples);
	}

	public static File getExamplesDir() {
		File examples = StandardFileChooser.getRecentDir("examples", true);
		String[] folderNames = new String[] { "RNAstructureExamples", "RNAstructure-Examples", "RNAstructure-examples",  "RNAstructure/Examples", "RNAstructure/examples", "examples" };
		if (examples == null)
			examples = PathTools.getFirstExisting(OSInfo.getDocumentsDir(), true, false, folderNames);

		if (examples == null)
			examples = PathTools.getFirstExisting(new File(System.getProperty("user.home")), true, false, folderNames);

		if (examples == null)
			examples = PathTools.findLocalPath("examples");
		if (examples != null && examples.exists())
			return examples;
		return null;
	}

	public static boolean silentMode() {
		return asBool(System.getProperty("silentLaunch"), false);
	}

	private static void parseCommandArgs(String[] args) {
		try {
			for(int i = 0; i < args.length; i++) {
				switch (args[i].toLowerCase()) {
					case "-d":  // debug.
						System.setProperty("log-debug", "true");
						System.setProperty("show-gui-info-menu","true");
						break;
					case "-s": System.setProperty("silentLaunch","true"); break;
					case "-v": System.setProperty("log-verbosity", args[++i]);
					case "-sfc": AppMainFrame.setUseSimpleFileChooser(true); break;
					case "-jlaf": options.useNativeLAF = false; break;
					case "-jmnu": options.useNativeMenus = false; break;
					case "-prefs": preferencesNode = args[++i]; break;
					case "-opt": {// override a user option (stored in preferences)
						int pos = args[++i].indexOf('=');
						if (pos == -1)
							throw new InvalidPreferencesFormatException("The -opt flag requires an argument in the form \"name=value\".");
						UserOptions.override(args[i].substring(0, pos), args[i].substring(pos + 1));
					}
					break;
					default:
						log.error("Unknown command-line option: \"%s\" (position: %s).",  args[i], i);
				}
			}
		} catch (Throwable e) {
			log.error("error parsing command line.", e);
		}
	}

	static String findDataPath() {
		String paths[] = { ".", "..", "./resources", "../Resources" };
		String fullPaths[] = { "{path}/", "{path}/data_tables", "{app}/{path}/", "{app}/{path}/data_tables" };
		String glob = "{*.specification,autodetect}.dat"; // find any file of the form: *.specification.dat  or  autodetect.dat

		String appPath = PathTools.getAppPath(RNAstructure.class);
        String realpath;
		for (String path : paths) {
			for (String full : fullPaths) {
                String fullPath = full;
                if (appPath!=null) fullPath=fullPath.replace("{app}", appPath);
				fullPath=fullPath.replace("{path}", path);
                File dir = new File(fullPath);
                //System.out.println("Trying: " + fullPath + " ==> " + dir.getAbsolutePath());
                if (dir.exists())
                    try {
                        realpath = dir.getCanonicalPath();
                        if (PathTools.anyFilesMatch(realpath, glob))
                            return realpath;
                    } catch (Exception ex) {
                        //Do nothing. continue the loop.
                    }
            }
		}
		return null;
	}

	static String getFullPath(String path) {
		File f = new File(path);
		if (isEmpty(path))
			return "<EMPTY>";

		if (f.exists())
			try {
				return f.getCanonicalPath();
			} catch(IOException ex) {
				return "ERROR: " + ex.getMessage();
			}
		else
			return "<File not found>";
	}
	static boolean fileExists(String path) { return new File(path).exists(); }

	static boolean verifyDataPath() {
		final String APP_PATH = "{APP_PATH}";

		String path = System.getenv("DATAPATH");
		if (path == null) path = "";
		if (path.contains(APP_PATH)) {
			String appPath = PathTools.getAppPath(RNAstructure.class);
			if (appPath == null) appPath = ".";
			log.debug("Using " + APP_PATH + " => " + appPath);
			log.debug("Before replacement: DATAPATH=" + path);
			path = path.replace(APP_PATH, appPath);
			log.debug("After replacement: DATAPATH=" + path);
			setDataPath(path);
		}
		log.debug(String.format("DATAPATH='%s' (%s)", path, getFullPath(path)));

		if (fileExists(path))
			return true;

		// DATAPATH was either NULL or not found. In either case, search for one in standard directories.
		path = findDataPath();
		if (path == null)
			showWarning("The thermodynamic data tables could not be found.\nPlease set the DATAPATH environment variable to the directory where they are stored.\n\t(e.g. <PATH-TO-RNAstructure>/data_tables).");
		else {
			log.debug("The thermodynamic data tables were located at \"" + path + "\". The DATAPATH environment variable will be set to this location.");
			setDataPath(path);
		}
		return path != null;
	}

	static String arrJoin(Object[] arr, String sep) {
		StringBuilder sbStr = new StringBuilder();
		for (Object o : arr) {
			sbStr.append(o.toString());
			sbStr.append(sep);
		}
		if (sbStr.length() > 0)
			sbStr.setLength(sbStr.length()-sep.length());
		return sbStr.toString();
	}

	static void setDataPath(String path) {
		RNAstructureBackendCalculator.setEnvVar("DATAPATH=" + path);
		System.setProperty("DATAPATH", path);
        log.debug("DATAPATH Verification (setDataPath): " + RNAstructureBackendCalculator.getEnvVar("DATAPATH"));
	}


	static void showWarning(String message) {
		log.error("Warning: " + message);
		if (!silentMode())
			Dialogs.showWarning(message, true);
	}
	static void showInfo(String message) {
		log.info("Information: " + message);
		if (!silentMode())
			Dialogs.showMessage(message);
	}
	public static void ExtractExamples() {
		try {
			String targetPath;
			File src = PathTools.findLocalPath("examples");
			if (src == null) {
				Dialogs.showError("The example files source folder could not be found."); //, but you can download them online and extract them manually
				return;
			}
			File docs = OSInfo.getDocumentsDir();
			if (new File(docs, "RNAstructure").exists())
				targetPath = "/RNAstructure/Examples";
			else
				targetPath = "RNAstructure-Examples";
			targetPath = docs.getPath() + targetPath.replace('/', File.separatorChar);

			String path = StandardFileChooser.getFileName(true, null, "Extract Examples Files", targetPath, "*", null);
			if (path == null) return;
			File[] list = src.listFiles();
			if (list == null) return;

			File dir = new File(path);
			if (!dir.mkdirs()) throw new IOException("Failed to create requested directory: " + path);
			Path dst = dir.toPath();
			for (File sub : list)
				try {
					Files.copy(sub.toPath(), dst.resolve(sub.getName()) , StandardCopyOption.REPLACE_EXISTING);
				} catch (IOException e) {
					e.printStackTrace();
				}
			StandardFileChooser.setRecentDir("examples", dir);
			StandardFileChooser.setDefaultDir(dir);
		} catch (Exception ex) {
			Dialogs.showError("Failed to extract example files: " + ex.toString());
		}
	}
	private static String examplesDir;
	public static boolean isExampleFile(final File f) {
		try {
			if (examplesDir == null) {
				File dir = getExamplesDir();
				if (dir == null)
					examplesDir = "<>"; //cannot exist.
				else
					examplesDir = dir.getCanonicalPath();
			}
			return examplesDir.equals(f.getCanonicalFile().getParent());
		} catch (IOException ex) {
			return false;
		}
	}

	/**
	 * Compares two version strings.
	 *
	 * Use this instead of String.compareTo() for a non-lexicographical
	 * comparison that works for version strings. e.g. "1.10".compareTo("1.6").
	 *
	 * Note: It does not work if "1.10" is supposed to be equal to "1.10.0".
	 *
	 * @param str1 a string of ordinal numbers separated by decimal points.
	 * @param str2 a string of ordinal numbers separated by decimal points.
	 * @return The result is a negative integer if str1 is _numerically_ less than str2.
	 *         The result is a positive integer if str1 is _numerically_ greater than str2.
	 *         The result is zero if the strings are _numerically_ equal.
	 */
	public static int versionCompare(String str1, String str2) {
		String[] vals1 = str1.split("\\.");
		String[] vals2 = str2.split("\\.");
		int i = 0;
		// set index to first non-equal ordinal or length of shortest version string
		while (i < vals1.length && i < vals2.length && vals1[i].equals(vals2[i])) {
			i++;
		}
		// compare first non-equal ordinal number
		if (i < vals1.length && i < vals2.length) {
			int diff = Integer.valueOf(vals1[i]).compareTo(Integer.valueOf(vals2[i]));
			return Integer.signum(diff);
		}
		// the strings are equal or one string is a substring of the other
		// e.g. "1.2.3" = "1.2.3" or "1.2.3" < "1.2.3.4"
		return Integer.signum(vals1.length - vals2.length);
	}

	/**
	 * This is called from within AppMainFrame during startup when
	 * the program is being run for the first time. If the user
	 * deletes it, it will not be created again.
	 */
	public static void createDocsFolder() {
		File dir = new File(OSInfo.getDocumentsDir(), "RNAstructure");
		if (!dir.exists())
			dir.mkdirs();
	}
}
