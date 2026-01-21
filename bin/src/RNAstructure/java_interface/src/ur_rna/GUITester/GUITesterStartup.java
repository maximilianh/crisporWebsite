/**
 * 
 */
 package ur_rna.GUITester;

import ur_rna.GUITester.GuiTools.GuiAppManager;
import ur_rna.Utilities.*;
import ur_rna.Utilities.annotation.*;


import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * The main entry-point for the test application. 
 */
@ApplicationInfo(name = "GUITester", title="GUITester", version = "1.0.1")
public class GUITesterStartup {
	//package-private
	static final AppOptions options = new AppOptions();
	static final AppLog log = new AppLog();

	//package-private
	static class AppOptions {
		public Verbosity logVerbosity = Verbosity.Info;
		public final ArrayList<Pair<String,ScriptSource>> scriptList = new ArrayList<>();
		public boolean parseOnly;
	}

	/**
	 * @param args Command line arguments.
	 */
	public static void main(String[] args) {
		System.setProperty("abbot.robot.default_delay", "10000"); // set default abbot Robot delay (time to wait for GUI components to be ready) to 10 seconds.

		//ClassLoader cl = ClassLoader.getSystemClassLoader();
		//URL[] urls = ((URLClassLoader)cl).getURLs();
		//for(URL url: urls){
		//	System.out.println("CP: " + url.getFile());
		//}

		// Parse the command line to set fields in the AppOptions object.
		// This includes adding script arguments to the AppOptions scriptList.
		if (!parseArgs(args)) return;
		log.setVerbosity(options.logVerbosity); // sets default verbosity (as defined in AppOptions)
		log.readSystemProperties();

		log.debug("CP:" + System.getProperty("java.class.path"));

		if (options.scriptList.size() == 0) {
			log.error("No test scripts were specified.");
			return;
		}

		int errorCode = 0;
		GuiAppManager app = new GuiAppManager(log);
		try {
			log.debug("UI Testing started.");
			for (Pair<String, ScriptSource> s : options.scriptList) {
				String fileName = s.getKey();
				ScriptSource srcType = s.getValue();
				ArrayList<Script> scripts = new ArrayList<>(1);

				//String scriptText = null;
				if (srcType == ScriptSource.INLINE) {
					// The script has been specified in-line. fileName is actually the full text.
					scripts.add(new Script(null, new StringReader(fileName)));
				} else {
					try {
						String baseName = PathTools.getBaseName(fileName);
						scripts.addAll(splitScriptSuite(new FileReader(fileName), fileName, baseName));
					} catch (FileNotFoundException e) {
						errorCode = 1;
						log.error("Script file not found. Path: " + fileName, e);
						continue;
					}
				}
				for (Script script : scripts) {
					try {
						log.debug("Running UI test script: " + script.toString());
						new ScriptRunner(script, app, log).run();
					} catch (Exception ex) {
						log.error("Exception in Script file: " + script.toString(), ex);
						errorCode = 1;
					} finally {
						script.input.close();
					}
				}
			}
		} catch (Throwable ex) {
			log.error("Unhandled error during test. Test aborted.", ex);
			errorCode = 10;
		} finally {
			app.close();
			app.shutdown();
			log.debug("UI Testing ended.");
		}
		if (errorCode != 0)
			System.exit(errorCode);
	}
	private static boolean parseArgs(String[] commandArgs) {
		AppOptions opts = GUITesterStartup.options;
		CommandLineParser p = new CommandLineParser(true); //ignore the case of flags
		try {
			p.parseFlagDefs("help,h  debug,d  trace,t  log,l:R  suite,s:R  inline,i:R noexec");
			CommandLineParser.ParseResults r = p.parse(commandArgs);

			if (r.hasFlag("help")) {
				ShowAppHelp();
				return false;
			}

			// If the -L flag was specified, set the log-level to the flag's parameter.
			// (The value will be null if it was not found, which will result in the default value, LOG_INFO being used.)
			if (r.hasFlag("log"))
				opts.logVerbosity = Verbosity.fromImportance(Convert.toInt(r.getFlagValue("log"), Verbosity.Info.importance));

			if (r.hasFlag("debug")) opts.logVerbosity = Verbosity.Debug;
			if (r.hasFlag("trace")) opts.logVerbosity = Verbosity.Trace;
			opts.parseOnly = r.hasFlag("noexec");

			ArrayList<Pair<String,ScriptSource>> scripts = opts.scriptList;
			for (CommandLineParser.Argument a : r.arguments) {
				String source = a.fullText;
				ScriptSource srcType;
				if (a instanceof CommandLineParser.UnnamedArg)
					srcType = ScriptSource.FILE;
				else if (a instanceof CommandLineParser.FlagArg) {
					CommandLineParser.FlagArg f = (CommandLineParser.FlagArg)a;
					source = f.valueText;
					switch (f.flag.name) {
						case "suite": srcType= ScriptSource.SUITE_FILE; break;
						case "inline": srcType= ScriptSource.INLINE; break;
						default:
							continue;
					}
				} else
					continue; // There are only two types of Arguments, so we'd never get here.

				scripts.add(new Pair<String, ScriptSource>(source, srcType));
			}

			return true;
		} catch (SyntaxErrorException e) {
			log.error("Command-line error: " + e.getMessage());
			log.info(p.getUsageMessage());
			return false;
		}
	}

	private enum ScriptSource {
		FILE,
		SUITE_FILE,
		INLINE
	}

	private static void ShowAppHelp() {
		System.out.println(GUITesterStartup.class.getPackage().getName() + "  [-L <log-level>] [-h|-help] ( <file-path> | @\"<script-text>\" ...)*");
	}

	private static List<Script> splitScriptSuite(Reader input, String file, String defaultTitle) throws IOException {
		LineNumberReader r = new LineNumberReader(input);
		StringBuilder sb = new StringBuilder();
		ArrayList<Script> scripts = new ArrayList<>();
		long startOffset = r.getLineNumber();
		long endOffset;
		boolean scriptStarted = false;
		String s, title = null;
		Pattern p = Pattern.compile("<--([^\\n]*)-->");
		boolean atEOF = false;
		while (!atEOF) {
			// this gets the position BEFORE the line is read, because if the line matches
			// the pattern, the NEW line number indicates the start of the NEXT script, not the current one.
			endOffset = r.getLineNumber();

			s = r.readLine();
			if (s == null) {
				s = "<--FINAL-->";
				atEOF = true;
			}

			Matcher m = p.matcher(s);
			if (m.matches()) {
				// This could either be on the first line (in which case there is NO previous script
				//   or there could be a script above this, possibly with no header-line.
				if (scriptStarted) {
					if (title == null) title = defaultTitle + " part #" + (scripts.size() + 1);
					String header = "//Script: " + title + "\n";
					Script scr = new Script(file, new StringReader(header + sb.toString()), startOffset, endOffset);
					scr.title = title;
					scripts.add(scr);
					sb.setLength(0); // discard previous text
					startOffset = endOffset + 1;
				}
				//get title for upcoming script
				title = m.group(1);
			} else {
				if (scriptStarted)
					sb.append('\n');
				sb.append(s);
				scriptStarted = true;
			}
		}
		return scripts;
	}
//		int errors = 0;
//		AppManager app = new AppManager();
//		for (Script ts : scripts) {
//			try {
//				AppLog.info("Running Test: %s.", ts.sourceInfo);
//				ts.parseTree = ParseScript(ts);
//				ExecuteScript(ts, app);
//				AppLog.info("Test Successful: %s.", ts.sourceInfo);
//				catch (ScriptRuntimeException ex) {
//					System.out.flush();
//					AppLog.error(String.format("Runtime error in Script %s.", ts.sourceInfo), ex);
//					errors++;
//				}
//				System.out.flush();
//				System.err.flush();
//
//				if (errors != 0) break; //just during initial testing
//			}
//			return scripts;
//		}
//	private BufferedReader buffered(final Reader input) {
//		if (input instanceof BufferedReader)
//			return (BufferedReader)_input;
//		return new BufferedReader(_input);
//	}
}
