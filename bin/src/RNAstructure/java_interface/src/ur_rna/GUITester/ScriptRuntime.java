package ur_rna.GUITester;

import abbot.finder.Matcher;
import ur_rna.GUITester.GuiTools.GuiAppManager;
import ur_rna.GUITester.GuiTools.GuiItemNotFoundException;
import ur_rna.GUITester.GuiTools.GuiItemRef;
import ur_rna.GUITester.GuiTools.GuiRelative;
import ur_rna.GUITester.GuiTools.Matchers.*;
import ur_rna.GUITester.ScriptParser.*;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.Verbosity;

import java.io.IOException;
import java.util.*;

import static ur_rna.GUITester.RuntimeTools.*;
import static ur_rna.Utilities.ObjTools.toDisplayString;
import static ur_rna.Utilities.ObjTools.toStr;

class ScriptRuntime {
    private HashMap<String, Object> _vars = new HashMap<String, Object>();
    private Script _script;
    private AppLog _log;
    private GuiAppManager _gui;

    private final boolean DEBUG_GUI_MATCHER = false;

    private boolean _doParsingOnly;
    private ScriptActions _actions;
    private boolean _disableTrace;
    private int delayBetweenSteps = 0;

    public void trace(String s, Object ... args) {
        if (!_disableTrace) _log.trace(s, args);
    }
    public void trace(String s) { if (!_disableTrace) _log.trace(s); }
    public boolean shouldTrace() {
        return _log.isEnabled(Verbosity.Trace) && ! _disableTrace;
    }
    public Script getScript() { return _script; }
    public AppLog getLog() { return _log; }
    public GuiAppManager getGui() { return _gui; }

    public ScriptRuntime(Script script, GuiAppManager guiApp, AppLog log) {
        _script = script;
        _gui = guiApp;
        _log = log;
        _doParsingOnly = GUITesterStartup.options.parseOnly;
        _actions = new ScriptActions(this);
    }

	public Object tryGetVar(String name)  {
        if (name.startsWith("ENV_"))
            return System.getenv(name.substring(4));
        if (name.startsWith("SV_"))
            return getSpecialValue(name.substring(3));
        return _vars.getOrDefault(name, ObjTools.MISSING);
    }
    private Object getSpecialValue(final String name) {
        switch (name.toLowerCase()) {
            case "timer": return System.currentTimeMillis();
            case "title": return _script.title;
            case "date": {
                Calendar cal = Calendar.getInstance(); // locale-specific
                cal.set(Calendar.HOUR_OF_DAY, 0);
                cal.set(Calendar.MINUTE, 0);
                cal.set(Calendar.SECOND, 0);
                cal.set(Calendar.MILLISECOND, 0);
                return cal.getTime();
            }
            case "now": return new Date();
            case "pwd": case "cwd": try { return new java.io.File( "." ).getCanonicalPath(); } catch (IOException ex) { return ""; }
            default:
                return ObjTools.MISSING;
        }
    }

    public Object getVar(String name) throws ScriptRuntimeException {
        Object value = tryGetVar(name);
        if (value == ObjTools.MISSING) throw new ScriptRuntimeException("Variable is undefined: '" + name + "'.");
        return value;
    }

    public void setVar(String name, Object value) {
        if (name.startsWith("ENV_"))
            throw new UnsupportedOperationException("Setting environment variables is not supported in Java.");

        if (name.startsWith("SV_")) {
            switch(name.substring(3)) {
                case "DELAY":
                    delayBetweenSteps = ((Number) value).intValue();
                    break;
                default:
                    throw new UnsupportedOperationException("Cannot set the value of Special Variable '" + name.substring(3) + "'.");
            }
        }

        _vars.put(name, value);
    }

    public void exec(ScriptNode scriptRoot) throws ScriptRuntimeException {
        boolean exitOK = false;
        try {
            java.util.Objects.requireNonNull(scriptRoot);
            for (Node n : scriptRoot.getChildren()) {
                if (n == null) continue;
                ScriptNode node = (ScriptNode) n;
                try {
                    if (node instanceof ASTAssignment) {
                        exitOK = false;
                        executeScriptAssignment((ASTAssignment) node);
                    } else if (node instanceof ASTAction) {
                        exitOK = false;
                        executeScriptAction((ASTAction) node);
                    }
                    doStepDelay();
                    exitOK = true;
                } catch (ScriptRuntimeException ex) {
                    ex.setNode(node, false);
                    throw ex;
                } catch (Exception ex) {
                    throw new ScriptRuntimeException(ex.getMessage(), ex, node);
                }
            }
        } finally {
            if (!exitOK)
                _log.debug("Aborting Script Execution due to Error");
            _gui.close();
        }
    }

    private void doStepDelay() {
        if (delayBetweenSteps != 0)
            try {
                Thread.sleep(delayBetweenSteps);
            } catch (InterruptedException ignored) {
                // exit sleep
            }
    }

    private void executeScriptAction(ASTAction node) throws ScriptRuntimeException {
        if (shouldTrace())
            _log.trace("Action: %s called with (%s).", node.verb, buildArgListString(node));
        if (_doParsingOnly) return;
        boolean good = false;
        try {
            String action = node.verb;
            Object[] args = evaluateScriptArgList(node, false).values().toArray();
            _actions.performAction(action, args);
            good = true;
        } catch (RuntimeTools.UnsupportedActionException ex) {
            throw new ScriptRuntimeException("Unknown script action: '" + node.verb + "'.", null, node);
        } catch (ActionException ex) {
            _log.trace("ActionException" + ex.getMessage());
            ex.setMessage("Error performing script action '" + node.verb + "': " + ex.getMessage());
            ex.setNode(node, false);
            throw ex;
        } catch (GuiItemNotFoundException ex) {
            _log.trace("Gui Item Not Found:" + ex.getMessage());
            throw new ScriptRuntimeException(ex.getMessage(), ex, node);
        } finally {
            if (good)
                _log.trace("Action Completed Successfully.");
            else {
                String message = "HAD ERROR!!";
                _log.trace(message);
                _log.error(message);
            }
        }
    }
//    private String[] toStringArray(final Object[] args) {
//        String[] sargs = new String[args.length];
//        for(int i = 0; i < args.length; i++)
//            sargs[i] = toStr(args[i], null);
//        return sargs;
//    }

//    private void requireActionArgs(String action, int actualLen, int numArgs)  throws ScriptRuntimeException { requireActionArgs(action, actualLen, numArgs, numArgs);}
//    private void requireActionArgs(String action, int actualLen, int minArgs, int maxArgs) throws ScriptRuntimeException {
//        if (actualLen < minArgs) throw new ScriptRuntimeException(String.format("The action \"%s\" requires at least %d arguments, but only %d was given.", action, minArgs, actualLen));
//        if (maxArgs != -1 && actualLen > maxArgs) throw new ScriptRuntimeException(String.format("The action \"%s\" can receive at most %d argument(s), but it was passed %d argument(s).", action, maxArgs, actualLen));
//    }

    private void executeScriptAssignment(ASTAssignment node) throws ScriptRuntimeException {
        if (shouldTrace())
            _log.trace("Assignment: " +
                    node.subject.toString() + " assigned " +
                    toDisplayString(safeEvaluateScriptNode(node.newValue)));
        if (_doParsingOnly) return;

        if (node.subject instanceof ASTVarRef) {
            Object value = node.newValue == null ? null : evaluateScriptNode(node.newValue);
            if (value instanceof GuiItemRef) {
                // if we assign the GuiItemRef to a variable, we should
                // attempt to evaluate it on the spot, so that if conditions change later
                // and the GuiItemRef no longer matches the original criteria,
                // the variable should STILL contain the reference to the
                // Component that was valid at the time of assignment.
                GuiItemRef ref = (GuiItemRef)value;
                try {
                    _gui.guiFinder().find(ref, false); // this stores the current value, but does not overwrite any existing one.
                } catch (GuiItemNotFoundException ex) {
                    //no need to handle this.
                }
            }
            setVar(((ASTVarRef) node.subject).varName, value);
        } else
            throw new ScriptRuntimeException("Cannot assign to entity of type " + node.getNodeName());
    }

    private Object safeEvaluateScriptNode(ScriptNode n)  {
        if (n == null) return null;
        try {
            return evaluateScriptNode(n);
        } catch (ScriptRuntimeException ex) {
            return ex;
        }
    }

    private Object evaluateScriptNode(ScriptNode n) throws ScriptRuntimeException {
        if (n == null)
            throw new IllegalArgumentException("Cannot evaluate a null ScriptItem");

        if (n instanceof ASTVarRef)
            return getVar(((ASTVarRef) n).varName);
        if (n instanceof ASTLiteral)
            return ((ASTLiteral)n).getValue();
        if (n instanceof ASTGuiReference)
            return createGuiReference((ASTGuiReference) n);

        throw new ScriptRuntimeException("Unable to evaluate node type '" + n.getClass().getName() + "'.", null, n);
    }

//    private Component evaluateGuiReference(GuiItemRef ref, boolean supressError) throws ScriptRuntimeException {
//        try {
//            return ref.find(!supressError);
//        } catch (GuiItemNotFoundException ex) {
//            throw new ScriptRuntimeException(ex.getMessage(), ex);
//        }
//    }

    private GuiItemRef createGuiReference(ASTGuiReference n) throws ScriptRuntimeException {
        GuiItemRef ref;
        GuiRelative rel = convertGuiSearchRelationship(n.relationshipKind);

        if (n.guiItem instanceof ASTVarRef) {
            Object value = getVar(((ASTVarRef) n.guiItem).varName);
            if (value instanceof GuiItemRef)
                ref = (GuiItemRef)value;
            else
                throw new ScriptRuntimeException("The variable '"+((ASTVarRef) n.guiItem).varName+"' does not contain a valid GUI Component.", null, n.guiItem);
        } else if (n.guiItem instanceof ASTGuiSearchCriteria) {
            ASTGuiSearchCriteria c = (ASTGuiSearchCriteria) n.guiItem;

            if (DEBUG_GUI_MATCHER && shouldTrace()) {
                String traceStr = String.format("Creating GUI item reference. Search Criteria: {{ %s }}", buildArgListString(c));
                //if (rel != null) traceStr += "( " + rel.name() + " of " + relative + " )";
                _log.trace(traceStr + ".");
            }

            try {
                ref = new GuiItemRef(buildMatcher(c));
            } catch (ScriptRuntimeException ex) {
                ex.setNode(n.guiItem, false);
                throw ex;
//            } catch (Throwable ex) {
//                throw new ScriptRuntimeException("error constructing Gui item reference: {" + buildArgListString(c) + "}", ex, n.guiItem);
            }
        } else
            throw new ScriptRuntimeException("Invalid gui item type: " + (n.guiItem == null ? "null" : n.guiItem.getClass().getName()) + ".",null , n.guiItem);

        if (rel != null)
            ref.setRelative(createGuiReference(n.relative), rel);

        if (DEBUG_GUI_MATCHER) trace("Created " + ref.toString(false) + ".");
        return ref;
    }

    private String buildArgListString(ScriptNode node) {
        StringBuilder sb = new StringBuilder();
        Map<String,Object> args = null;

        // If tracing is turned on, turn it off temporarily.
        _disableTrace = true;
        try {
            // evaluateScriptArgList will NOT throw an error when useSafeEval is true.
            args = evaluateScriptArgList(node, true);
        } catch (ScriptRuntimeException ex) {
            ex.printStackTrace();
            return String.format("{ ERROR: %s }", ex.getMessage().replace("\n", "; ").replace("\r", ""));
        } finally {
            _disableTrace = false;
        }

        int pos = 0;
        for (Map.Entry<String,Object> arg : args.entrySet()) {
            if (pos != 0)
                sb.append(", ");
            if (!arg.getKey().startsWith("#"))
                sb.append(arg.getKey()).append(": ");
            sb.append(toDisplayString(arg.getValue()));
            pos++;
        }
        return sb.toString();
    }

    private static GuiRelative convertGuiSearchRelationship(int rel) {
        switch (rel) {
            case 0: return null;
            case GuiTestScriptParserConstants.CHILD: return GuiRelative.Child;
            case GuiTestScriptParserConstants.PARENT: return GuiRelative.Parent;
            case GuiTestScriptParserConstants.ANCESTOR: return GuiRelative.Ancestor;
            case GuiTestScriptParserConstants.DESCENDANT: return GuiRelative.Descendant;
            case GuiTestScriptParserConstants.SIBLING: return GuiRelative.Sibling;
            //case GuiTestScriptParserConstants.PREV_SIBLING: return Relationship.PrevSibling;
            case GuiTestScriptParserConstants.REL_LABEL: return GuiRelative.LabelTarget;
            default: throw new IllegalArgumentException("Invalid GUI component relationship: " + rel);
        }
    }

    private ClassMatcher getClassMatcher(String name) throws ScriptRuntimeException {
        switch (name.toLowerCase()) {
            case "menu": return SwingClassMatchers.menu;
            case "button": return SwingClassMatchers.button;
            case "text":
            case "input": return SwingClassMatchers.text;
            case "textfield":
            case "field": return SwingClassMatchers.field;
            case "textarea": return SwingClassMatchers.textArea;

            case "tick":
            case "check": return SwingClassMatchers.check;

            case "radio":
            case "option": return SwingClassMatchers.radio;

            case "list": return SwingClassMatchers.list;

            case "dropdown":
            case "droplist":
            case "combo": return SwingClassMatchers.combo;

            case "label": return SwingClassMatchers.label;
            case "panel": return SwingClassMatchers.panel;

            case "window": return SwingClassMatchers.window;
            case "dialog": return SwingClassMatchers.dialog;

            case "spin":
            case "spinner":
            case "updown":return SwingClassMatchers.spinner;

            default: {
                ClassMatcher cm = ClassMatcher.ofType(name);
                if (cm == null)
                    throw new ScriptRuntimeException("Unknown GUI component type: " + name);
                return cm;
            }
        }
    }

    private abbot.finder.Matcher buildMatcher(ASTGuiSearchCriteria searchCriteria) throws ScriptRuntimeException {
        Map<String, Object> criteria = evaluateScriptArgList(searchCriteria, false);
        ComposableMatcher m = ComposableBase.TRUE;

        for (Map.Entry<String, Object> c : criteria.entrySet()) {
            String key = c.getKey();
            if (key.startsWith("#")) key = null;
            m = appendGuiSearchParam(m, key, c.getValue());
        }

        return m;
    }

    private Map<String,Object> evaluateScriptArgList(ScriptNode listNode, boolean useSafeEval) throws ScriptRuntimeException {
        LinkedHashMap<String, Object> list = new LinkedHashMap<String, Object>();
        if (listNode instanceof  ASTGuiSearchCriteria) {
            String type = ((ASTGuiSearchCriteria)listNode).itemType;
            if (type != null)
                list.put("@Type", type);
        }
        int pos = 0;
        for (Node arg : listNode.getChildren()) {
            String name= "#"+(++pos);
            ScriptNode argVal;
            if (arg instanceof ASTNamedArgument) {
                ASTNamedArgument nArg = (ASTNamedArgument) arg;
                argVal = nArg.val;
                if (nArg.name != null)
                    name = nArg.name;
            } else
                argVal = (ScriptNode)arg;

            Object value;
            if (useSafeEval)
                value = safeEvaluateScriptNode(argVal);
            else
                value = evaluateScriptNode(argVal);

            list.put(name, value);
        }
        return list;
    }

    private ComposableMatcher appendGuiSearchParam(ComposableMatcher subject, String argName, Object argVal) throws ScriptRuntimeException {
        argName = argName == null ? "" : argName;
        String sval = argVal == null ? "" : argVal.toString();
        switch (argName.toLowerCase()) {
            case "":
            case "*":
            case "desc":
                return subject.and(new DescriptionMatcher(sval));
            case "type":
            case "@type":
                return subject.and(getClassMatcher(toStr(argVal)));
            case "name":
                return subject.and(new DescriptionMatcher(sval, DescriptionMatcher.SEARCH_NAME));
            case "text":
            case "value":
                return subject.and(new DescriptionMatcher(sval, DescriptionMatcher.SEARCH_INPUT));
            case "caption":
            case "title":
                return subject.and(new DescriptionMatcher(sval, DescriptionMatcher.SEARCH_CAPTION));
            case "comment":
            case "note":
            case "tip":
                return subject.and(new DescriptionMatcher(sval, DescriptionMatcher.SEARCH_COMMENT));
            case "selection":
            case "selected":
            case "sel":
                return subject.and(new DescriptionMatcher(sval, DescriptionMatcher.SEARCH_SELECTION));
            case "vis":
            case "visible": {
                ClassMatcher cm = FindClassMatcher(subject);
                switch(toStr(argVal).toLowerCase()) {
                    case "true":
                    case "1":
                        if (cm == null)
                            subject = subject.and(new VisibilityMatcher(true));
                        else
                            cm.allowHidden = false;
                        break;
                    case "false":
                    case "0":
                        if (cm != null)
                            cm.allowHidden = true;
                        subject = subject.and(new VisibilityMatcher(false));
                        break;
                    case "any":
                        if (cm != null)
                            cm.allowHidden = true;
                        break;
                    default:
                        throw new ScriptRuntimeException("Unknown Visibility value: " + argVal);
                }
                return subject;
            }
            default:
                throw new ScriptRuntimeException("Unknown GUI Search parameter: " + argName);
        }
    }
    private ClassMatcher FindClassMatcher(final Matcher subject) {
        if (subject instanceof ClassMatcher)
            return (ClassMatcher)subject;
        if (subject instanceof ComposableMatcher) {
            for (Matcher m : ((ComposableMatcher) subject).getChildMatchers()) {
                ClassMatcher cm = FindClassMatcher(m);
                if (cm != null) return cm;
            }
        }
        return null;
    }
}
