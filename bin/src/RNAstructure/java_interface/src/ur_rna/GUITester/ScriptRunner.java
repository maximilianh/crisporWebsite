package ur_rna.GUITester;

import ur_rna.GUITester.GuiTools.GuiAppManager;
import ur_rna.GUITester.ScriptParser.*;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.Verbosity;

import java.io.PrintStream;

public class ScriptRunner implements Runnable {
	private Script _src;
	private ScriptRuntimeException _runError;
	private ParseException _parseError;
	private GuiAppManager _subjectApp;
	private AppLog _log;
	private static final boolean TRACE_PARSE_TREE = false;

	public ScriptRunner(Script src, GuiAppManager subject, AppLog log) { _src = src; _subjectApp = subject; _log = log; }

	public Script getSource() { return  _src;}
	public ScriptRuntimeException getRunTimeError() { return _runError; }
	public ParseException getParseError() { return _parseError; }
	public Throwable getError() { return _runError == null ? _parseError : _runError; }
	public boolean hadError() { return _runError != null || _parseError != null; }

	@Override
	public void run() {
		try {
			ScriptNode ast = parse(false);
			exec(ast);
		} catch (ParseException ex) {
			_parseError = ex;
			_log.error(String.format("Syntax error in script %s.", _src.toString()), ex);
		} catch (ScriptRuntimeException ex) {
			_runError = ex;
			_log.error(String.format("Runtime error in script %s.", _src.toString()), ex);
		}
	}

	public ScriptNode parse(boolean enableTracing) throws ParseException {
		try {
			GuiTestScriptParser p = new GuiTestScriptParser(_src.getInput());
			if (enableTracing)
				p.enable_tracing();
			else
				p.disable_tracing();

			p.disable_tracing();

			return (ASTRoot) p.parse(); //throws ParseException
		} catch (TokenMgrError ex) {
			throw new ParseException("Internal Parser error", ex);
		}
	}

	public void exec(ScriptNode root) throws ScriptRuntimeException {
		if (TRACE_PARSE_TREE) {
			PrintStream tr = _log.getStream(Verbosity.Trace);
			if (_log.isTraceEnabled()) {
				tr.println("Parse Tree for " + _src.toString());
				root.dump("  ", tr);
			}
		}
		ScriptRuntime runtime = new ScriptRuntime(_src, _subjectApp, _log);
		runtime.exec(root);
	}
}


