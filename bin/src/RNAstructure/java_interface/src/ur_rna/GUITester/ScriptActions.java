package ur_rna.GUITester;

import ur_rna.GUITester.GuiTools.GuiAppManager;
import ur_rna.GUITester.GuiTools.GuiItemNotFoundException;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.Verbosity;

import java.io.PrintStream;

import static ur_rna.GUITester.RuntimeTools.*;
import static ur_rna.Utilities.ObjTools.toDisplayString;
import static ur_rna.Utilities.ObjTools.toStr;

/**
 * @author Richard M. Watson
 */
public class ScriptActions {
    final GuiAppManager _gui;
    final ScriptRuntime _rt;
    final AppLog _log;
    public ScriptActions(final ScriptRuntime rt) {
        _rt = rt;
        _gui = rt.getGui();
        _log = rt.getLog();
    }

    public void performAction(final String action, final Object[] args)
        throws GuiItemNotFoundException, ScriptRuntimeException {
        switch (action.toLowerCase()) {
            case "title":
                expectArgs(action, args, 1, 1);
                _rt.getScript().title = toStr(args[0]);
                break;
            case "cout":
            case "echo":
            case "message":
                System.out.println(joinArgs(args, false)); break;
            case "cerr":
            case "errmsg":
                System.err.println(joinArgs(args, false)); break;
            case "debug": _log.debug(joinArgs(args, true)); break;
            case "trace": _log.trace(joinArgs(args, true)); break;
            case "log": _log.info(joinArgs(args, true)); break;
            case "warn": _log.warn(joinArgs(args, true)); break;
            case "logerr": _log.error(joinArgs(args, true)); break;

            case "sleep":
                expectArgs(action, args, 1, 1);
                expectArgType(action, args, 1, "millis", false, Number.class);
                try {
                    Thread.sleep(((Number) args[0]).longValue());
                } catch (InterruptedException ignored) {}
                break;
            case "startgui":
                expectArgs(action, args, 0, 1);
                if (args.length == 0)
                    _gui.start("");
                else {
                    expectArgType(action, args, 1, "args", false, String.class);
                    _gui.start((String)args[0]);
                }
                break;
            case "closegui":
                _gui.close();
                break;
            case "waitidle":
                _gui.guiRobot().waitForIdle();
                break;
            case "listgui":
                System.out.println(_gui.getHierarchy(null).toString());
                break;
            default:
                _gui.performGuiAction(action, args);
                break;
        }
    }
    private String joinArgs(final Object[] args, boolean quoted) {
        StringBuilder sb = new StringBuilder();
        for (Object o : args) {
            if (sb.length() != 0)
                sb.append("\t");
            if (quoted)
                sb.append(toDisplayString(o));
            else
                sb.append(o);
        }
        return sb.toString();
    }
}
