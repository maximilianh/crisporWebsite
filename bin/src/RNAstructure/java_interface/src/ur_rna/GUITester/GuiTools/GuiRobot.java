package ur_rna.GUITester.GuiTools;

import abbot.finder.Matcher;
import abbot.util.AWT;
import ur_rna.GUITester.GuiTools.Matchers.DescriptionMatcher;
import ur_rna.GUITester.GuiTools.Matchers.SwingClassMatchers;
import ur_rna.GUITester.RuntimeTools;
import ur_rna.Utilities.*;

import javax.swing.*;
import javax.swing.text.JTextComponent;
import java.awt.*;
import java.awt.event.KeyEvent;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.lang.reflect.Field;
import java.text.ParseException;
import java.util.Objects;
import java.util.concurrent.ThreadPoolExecutor;

import static ur_rna.GUITester.RuntimeTools.expectArgType;
import static ur_rna.GUITester.RuntimeTools.expectArgs;
import static ur_rna.Utilities.ObjTools.toStr;
import static ur_rna.Utilities.Strings.escapeStringLiteral;

public class GuiRobot extends abbot.tester.Robot {

    /**
     * Total number of milliseconds to wait for a GUI component, when
     * simply testing for its existence. (i.e. for "PROHIBIT" etc)
     */
    public static final int DEFAULT_TEST_WAIT_TIME = 250;
    /**
     * Total number of milliseconds to wait for a GUI component to be created, before returning an error.
     * Can be overridden by setting {@link #GuiWaitTime}.
     */
    public static final int DEFAULT_WAIT_TIME = 2000;
    /**
     * Number of milliseconds to sleep between attempts to poll for a component.
     */
    public static final int DEFAULT_POLL_TIME = 50;

    private GuiAppManager _owner;
    private AppLog _log;

    /* The default amount of time to wait for a required GUI item.
     */
    public int GuiWaitTime = DEFAULT_WAIT_TIME;

    public GuiRobot(final GuiAppManager owner) {
        _owner = owner;
        _log = owner.getLog();
        abbot.tester.Robot.setAutoDelay(1);
    }

    public void doAction(String action, Object[] args) throws RuntimeTools.ActionException, GuiItemNotFoundException {
        Component cmp;
        int i=0, arglen = args.length;
        Class guiClass = GuiItemRef.class;

        switch (action.toLowerCase()) {
            case "menu": // MENU "MainMenu", "SubMenu", ...
                expectArgs(action, args, 1, -1);
                GuiItemRef parentMenu = new GuiItemRef(_owner.getRootFrame().getJMenuBar()); // Also works, but is less direct: SwingClassMatchers.ofType(DynamicMenuBar.class.getName())
                for (i=0; i < arglen; i++) {
                    expectArgType(action, args, i + 1, "menu-item", false, String.class);
                    Matcher m = (i == 0 ? SwingClassMatchers.menu : SwingClassMatchers.menuItem)
                            .withText(toStr(args[i]), DescriptionMatcher.SEARCH_CAPTION);
                    GuiItemRef menu = new GuiItemRef(m);
                    menu.setRelative(parentMenu, GuiRelative.Child);
                    if (_log.isDebugEnabled()) _log.debug("Click Menu \"" + args[i] + "\"");
                    waitForIdle();
                    click(evalGuiRef(menu));
                    parentMenu = menu;
                }
                break;
            case "click":
                expectArgs(action, args, 1, 1);
                _log.debug("Click  \"" + args[0] + "\"");
                _log.trace("Click  \"" + evalGuiRef(args[0], false) + "\"");
                try {
                    click(evalGuiRef(args[0]));
                } catch (Throwable t) {
                    _log.trace("Click  Error!!");
                    throw t;
                }
                break;
            case "settext":
                expectArgs(action, args, 2, 2);
                expectArgType(action, args, 1, "guiRef", false, guiClass);
                _log.debug("Set Text of \"" + args[0] + "\"");
                setText(evalGuiRef(args[0]), toStr(args[1]));
                break;
            case "setvalue":
            case "setval":
                expectArgs(action, args, 2, 2);
                expectArgType(action, args, 1, "guiRef", false, guiClass);
                _log.debug("Set Value of \"" + args[0] + "\"");
                setValue(evalGuiRef(args[0]), args[1]);
                break;
            case "focus":
                expectArgs(action, args, 1, 1);
                expectArgType(action, args, 1, "guiRef", false, guiClass);
                _log.debug("Focus \"" + args[0] + "\"");
                focus(evalGuiRef(args[0]));
                break;
            case "cleartext":
                expectArgs(action, args, 0, 1);
                if (arglen == 1) {
                    _log.debug("Clear Text of \"" + args[0] + "\"");
                    expectArgType(action, args, 1, "guiRef", false, guiClass);
                    setText(evalGuiRef(args[0]), "");
                } else {
                    // Clear the text of the focused field.
                    _log.debug("Clear Text of Focused Control");
                    cmp = AWT.getFocusOwner();
                    if (cmp == null)
                        throw new RuntimeTools.ActionArgumentException("In ClearText: No component currently has the focus.");
                    setText(cmp, "");
                }
                break;
            case "typekeys":
                expectArgs(action, args, 1, -1);
                if (args[0] instanceof GuiItemRef) {
                    _log.debug("Focus " + args[0]);
                    expectArgs(action, args, 2, -1);
                    focus(evalGuiRef(args[0]));
                    i=1;
                }
                for (; i < arglen; i++) {
                    _log.debug("Type Key " + args[i]);
                    sendKeyString(toStr(args[i]));
                }
                break;
            case "typetext":
                expectArgs(action, args, 1, -1);
                if (args[0] instanceof GuiItemRef) {
                    _log.debug("Focus " + args[0]);
                    expectArgs(action, args, 2, -1);
                    focus(evalGuiRef(args[0]));
                    i=1;
                }
                for (; i < arglen; i++) {
                    _log.debug("Type Text " + args[i]);
                    keyString(toStr(args[i]));
                }
                break;
            case "typeenter": keyClick('\n'); _log.debug("Type ENTER"); break;
            case "typeesc": keyClick(KeyEvent.VK_ESCAPE); _log.debug("Type ESC "); break;
            case "wait": // Wait for a GuiRef search to succeed. I.e. waits for a control to come into existance.
                int timeout = GuiRobot.DEFAULT_WAIT_TIME; // 250 ms
                expectArgs(action, args, 1, 2);
                if (arglen == 2)
                    timeout = Integer.decode(args[1].toString());
                expectArgType(action, args, 1, "guiRef", false, guiClass);
                _log.debug("Wait for " + args[0] + " (time " + timeout + ")");
                evalGuiRef(args[0], true, false, timeout); // forces required evaluation. prints trace if failed.
                break;
            case "require":  // throw an error if the GuiRef is not found
                expectArgs(action, args, 1, -1);
                for (; i < arglen; i++) {
                    expectArgType(action, args, i + 1, "guiRef", false, guiClass);
                    _log.debug("Require: " + args[i]);
                    evalGuiRef(args[i]);
                }
                break;
            case "prohibit": // the opposite of "require" -- throw an error if found.
                expectArgs(action, args, 1, -1);
                for (; i < arglen; i++) {
                    expectArgType(action, args, i + 1, "guiRef", false, guiClass);
                    GuiItemRef ref = (GuiItemRef) args[i];
                    _log.trace("Prohibit: " + ref.toString());
                    ref.clearFound(true);
                    if (evalGuiRef(ref, false, false, DEFAULT_TEST_WAIT_TIME) != null) {
                        _log.trace("Prohibited item found!");
                        _owner.printHierarchyTrace();
                        throw new GuiItemNotFoundException("The GUI item was prohibited, but was found: " + ref.toString());
                    }
                }
                break;
            case "capture":
                expectArgs(action, args, 1, 2);
                expectArgType(action, args, 1, "filename", false, String.class);
                String file = toStr(args[0]);
                if (PathTools.getExt(file).length()==0)
                    file += ".png";

                if (arglen == 2) {
                    if (args[1] instanceof GuiItemRef)
                        cmp = AWT.getWindow(evalGuiRef(args[1]));
                    else if ("screen".equals(args[1]))
                        cmp = null;
                    else
                        throw new RuntimeTools.ActionArgumentException("Unknown screen-capture instruction:  " + toString(args[1]));
                } else
                    cmp = AWT.getFocusedWindow();
                try {
                    BufferedImage img = cmp == null ? ScreenCapture.captureDesktopImage() : ScreenCapture.captureImage(cmp);
                    ScreenCapture.writeImage(img, file);
                } catch (IOException ex) {
                    throw new RuntimeTools.ActionException("Failed to save screen-capture file: " + ex.getMessage(), ex);
                } catch (AWTException ex) {
                    throw new RuntimeTools.ActionException("Failed capture GUI component image: " + ex.getMessage(), ex);
                }
                break;
            case "forget":  // Clear the cached GuiRef
                // GuiRefs are cached in the script variable. Forget causes the cached value to be cleared,
                // so that a new search will be performed using the origincal criteria the next time the
                // GuiRef is used.
                expectArgs(action, args, 1, -1);
                for (; i < arglen; i++) {
                    expectArgType(action, args, i + 1, "guiRef", false, guiClass);
                    ((GuiItemRef) args[i]).clearFound(true);
                }
                break;
            default:
                throw new RuntimeTools.UnsupportedActionException(action);
        }
        try {
            waitForIdle(); //wait until the action has completed.
        } catch (UnsupportedOperationException ex) {
            // Seems that abbot has a bug that throws an exception:
            //    Caused by: java.lang.UnsupportedOperationException
            //    at java.lang.Thread.stop(Thread.java:872)
            //    at abbot.tester.Robot.waitForIdle(Robot.java:790)
            //    at abbot.tester.Robot.waitForIdle(Robot.java:908)
            ex.printStackTrace(_log.getTrStream());
        }
    }

    public void setText(Component cmp, String text)     throws RuntimeTools.ActionArgumentException {
        if (cmp instanceof JTextComponent)
            ((JTextComponent) cmp).setText(text);
        else if (cmp instanceof TextComponent)
            ((TextComponent) cmp).setText(text);
        else if (cmp instanceof JSpinner) {
            JComponent editor = ((JSpinner) cmp).getEditor();
            if (editor instanceof JSpinner.DefaultEditor) {
                JFormattedTextField field = ((JSpinner.DefaultEditor) editor).getTextField();
                field.setText(text);
                try { field.commitEdit(); } catch (ParseException ex) {
                    throw new RuntimeTools.ActionArgumentException(String.format("Parse Error when setting the text of a %s component.", cmp.getClass().getSimpleName()), ex);
                }
            }
        } else if (cmp instanceof JComboBox)
            ((JComboBox)cmp).setSelectedItem(text);
        else if (cmp instanceof JList)
            ((JList)cmp).setSelectedValue(text, true);

        else if (cmp == null)
            throw new RuntimeTools.ActionArgumentException("In setText: Component cannot be NULL.");
        else
            throw new RuntimeTools.ActionArgumentException(String.format("Setting the text of a %s component is not implemented.", cmp.getClass().getSimpleName()));
    }

    public void setValue(Component cmp, Object value)     throws RuntimeTools.ActionArgumentException {
        if (cmp instanceof JFormattedTextField)
            ((JFormattedTextField) cmp).setValue(value);
        if (cmp instanceof JTextComponent)
            ((JTextComponent) cmp).setText(value==null?"":value.toString());
        else if (cmp instanceof TextComponent)
            ((TextComponent) cmp).setText(value==null?"":value.toString());
        else if (cmp instanceof JSpinner) {
            SpinnerModel model = ((JSpinner)cmp).getModel();
            Object current = model.getValue();
            if (current!=null)
                value = Convert.toType(value, current.getClass());
            model.setValue(value);
        } else if (cmp instanceof JComboBox)
            ((JComboBox)cmp).setSelectedItem(value);
        else if (cmp instanceof JList)
            ((JList)cmp).setSelectedValue(value, true);

        else if (cmp == null)
            throw new RuntimeTools.ActionArgumentException("In setText: Component cannot be NULL.");
        else
            throw new RuntimeTools.ActionArgumentException(String.format("Setting the value of a %s component is not implemented.", cmp.getClass().getSimpleName()));
    }

    /**
     * These constants correspond to numbers used in abbot.Robot to set and test for modifier keys.
     * Unfortunately abbot does not define and enum or constants, but simply uses the plan, undocumented numbers.
     */
    @SuppressWarnings("MagicNumber")
    public enum ModifierKey {
        CTRL(2),
        SHIFT(1),
        ALT(8),
        META(4),
        ALTG(32);
        public final int code;
        ModifierKey(int code) {this.code = code;}
    }

    public void sendKeyString(String s) {
        final String logpfx = "  *";
        int len = s.length();
        char[] chars = s.toCharArray();
        boolean inSpecial = false;
        int modifiers = 0;
        StringBuilder specialPart = new StringBuilder();
        for (int i = 0; i < len; i++) {
            char c = chars[i];
            // "{CTRL}{HOME}{CTRL}{SHIFT}{END}"
            if (inSpecial) {
                if (c == '}') {
                    inSpecial = false;
                    String skey = specialPart.toString();
                    specialPart.setLength(0); //clear it
                    boolean isModifier = true;
                    switch (skey.toLowerCase()) {
                        case "ctrl": case "ctl": modifiers |= ModifierKey.CTRL.code; break;
                        case "shift": case "shft": modifiers |= ModifierKey.SHIFT.code; break;
                        case "alt": modifiers |= ModifierKey.ALT.code; break;
                        case "meta": modifiers |= ModifierKey.META.code; break;
                        case "altg": case "altgraph":modifiers |= ModifierKey.ALTG.code; break;
                        default: isModifier = false;
                    }
                    if (isModifier)
                        _log.trace(logpfx+"Added modifier: " + skey);

                    if (!isModifier) {
                        int key = getSpecialKey(skey);
                        if (modifiers != 0) _log.trace(logpfx+"Setting Modifers: " + modifiers);
                        setModifiers(modifiers, true);
                        if (key == 0) {
                            _log.trace(logpfx+"Sending KeyString: " + escapeStringLiteral(skey));
                            keyString(skey);
                        } else {
                            _log.trace(logpfx+"Sending KeyClick %s (%d)", KeyEvent.getKeyText(key),  key);
                            keyClick(key);
                        }
                        setModifiers(modifiers, false);
                        modifiers = 0;
                    }
                } else
                    specialPart.append(c);
            } else switch (c) {
                case '+': modifiers |= ModifierKey.SHIFT.code; _log.trace(logpfx+"Added modifier: SHIFT"); break;
                case '^': modifiers |= ModifierKey.ALT.code;  _log.trace(logpfx+"Added modifier: ALT"); break;
                case '!': modifiers |= ModifierKey.CTRL.code;  _log.trace(logpfx+"Added modifier: CTRL");  break;
                case '{': inSpecial = true; break;
                default:
                    if (modifiers != 0)_log.trace(logpfx+"Setting Modifers: " + modifiers);
                    setModifiers(modifiers, true);
                    _log.trace(logpfx+"Sending KeyStroke: %s (%d)", escapeStringLiteral(Character.toString(c)),  (int)c);
                    keyStroke(c);
                    setModifiers(modifiers, false);
                    modifiers = 0;
                    break;
            }
        }
    }

    public int getSpecialKey(String keyName) {
        int key;
        switch (keyName.toLowerCase()) {
            case "enter": key = KeyEvent.VK_ENTER; break;
            case "bs": case "backspace": case "bksp": key = KeyEvent.VK_BACK_SPACE; break;
            case "tab": key = KeyEvent.VK_TAB; break;
            case "cancel": key = KeyEvent.VK_CANCEL; break;
            case "clr":case "clear": key = KeyEvent.VK_CLEAR; break;
            case "compose": key = KeyEvent.VK_COMPOSE; break;
            case "pause": key = KeyEvent.VK_PAUSE; break;
            case "caps": key = KeyEvent.VK_CAPS_LOCK; break;
            case "esc":case "escape": key = KeyEvent.VK_ESCAPE; break;
            case "space": key = KeyEvent.VK_SPACE; break;
            case "pgup": key = KeyEvent.VK_PAGE_UP; break;
            case "pgdn": key = KeyEvent.VK_PAGE_DOWN; break;
            case "end": key = KeyEvent.VK_END; break;
            case "home": key = KeyEvent.VK_HOME; break;
            case "left": key = KeyEvent.VK_LEFT; break;
            case "up": key = KeyEvent.VK_UP; break;
            case "right": key = KeyEvent.VK_RIGHT; break;
            case "down": key = KeyEvent.VK_DOWN; break;
            case "begin": key = KeyEvent.VK_BEGIN; break;

            // modifiers
            case "shift": key = KeyEvent.VK_SHIFT; break;
            case "ctrl":case "control": key = KeyEvent.VK_CONTROL; break;
            case "alt": key = KeyEvent.VK_ALT; break;
            case "meta": key = KeyEvent.VK_META; break;
            case "altgraph": key = KeyEvent.VK_ALT_GRAPH; break;

            case "{": key = KeyEvent.VK_OPEN_BRACKET;  break; // '{'
            case "}": key = KeyEvent.VK_CLOSE_BRACKET;  break;//'{';
            case "bkslsh": key = KeyEvent.VK_BACK_SLASH; break; //'{';
            case "+": key = KeyEvent.VK_PLUS;  break;
            case "^": key = KeyEvent.VK_CIRCUMFLEX;  break;
            case "!": key = KeyEvent.VK_CONTROL;  break;

                // punctuation
//            case "comma": key = KeyEvent.VK_COMMA; break;
//            case "period": key = KeyEvent.VK_PERIOD; break;
//            case "slash": key = KeyEvent.VK_SLASH; break;
//            case "semicolon": key = KeyEvent.VK_SEMICOLON; break;
//            case "equals": key = KeyEvent.VK_EQUALS; break;
//            case "openBracket": key = KeyEvent.VK_OPEN_BRACKET; break;
//            case "backSlash": key = KeyEvent.VK_BACK_SLASH; break;
//            case "closeBracket": key = KeyEvent.VK_CLOSE_BRACKET; break;

                // numpad numeric keys handled below
            case "multiply": key = KeyEvent.VK_MULTIPLY; break;
            case "add": key = KeyEvent.VK_ADD; break;
            case "separator": key = KeyEvent.VK_SEPARATOR; break;
            case "subtract": key = KeyEvent.VK_SUBTRACT; break;
            case "decimal": key = KeyEvent.VK_DECIMAL; break;
            case "divide": key = KeyEvent.VK_DIVIDE; break;
            case "del": case "delete": key = KeyEvent.VK_DELETE; break;
            case "numLock": key = KeyEvent.VK_NUM_LOCK; break;
            case "scroll": key = KeyEvent.VK_SCROLL_LOCK; break;

            case "win": case "windows": key = KeyEvent.VK_WINDOWS; break;
            case "ctx":case "context": key = KeyEvent.VK_CONTEXT_MENU; break;

            case "f1": key = KeyEvent.VK_F1; break;
            case "f2": key = KeyEvent.VK_F2; break;
            case "f3": key = KeyEvent.VK_F3; break;
            case "f4": key = KeyEvent.VK_F4; break;
            case "f5": key = KeyEvent.VK_F5; break;
            case "f6": key = KeyEvent.VK_F6; break;
            case "f7": key = KeyEvent.VK_F7; break;
            case "f8": key = KeyEvent.VK_F8; break;
            case "f9": key = KeyEvent.VK_F9; break;
            case "f10": key = KeyEvent.VK_F10; break;
            case "f11": key = KeyEvent.VK_F11; break;
            case "f12": key = KeyEvent.VK_F12; break;
            case "f13": key = KeyEvent.VK_F13; break;
            case "f14": key = KeyEvent.VK_F14; break;
            case "f15": key = KeyEvent.VK_F15; break;
            case "f16": key = KeyEvent.VK_F16; break;
            case "f17": key = KeyEvent.VK_F17; break;
            case "f18": key = KeyEvent.VK_F18; break;
            case "f19": key = KeyEvent.VK_F19; break;
            case "f20": key = KeyEvent.VK_F20; break;
            case "f21": key = KeyEvent.VK_F21; break;
            case "f22": key = KeyEvent.VK_F22; break;
            case "f23": key = KeyEvent.VK_F23; break;
            case "f24": key = KeyEvent.VK_F24; break;

            case "prtscr": case "printscrn":case "printscr": case "printscreen": key = KeyEvent.VK_PRINTSCREEN; break;
            case "ins": case "insert": key = KeyEvent.VK_INSERT; break;
            case "help": key = KeyEvent.VK_HELP; break;
            case "backquote": key = KeyEvent.VK_BACK_QUOTE; break;
            case "quote": key = KeyEvent.VK_QUOTE; break;

//        case "up": key=KeyEvent.VK_KP_UP;break;
//        case "down": key=KeyEvent.VK_KP_DOWN;break;
//        case "left": key=KeyEvent.VK_KP_LEFT;break;
//        case "right": key=KeyEvent.VK_KP_RIGHT;break;

            case "deadgrave": key = KeyEvent.VK_DEAD_GRAVE; break;
            case "deadacute": key = KeyEvent.VK_DEAD_ACUTE; break;
            case "deadcircumflex": key = KeyEvent.VK_DEAD_CIRCUMFLEX; break;
            case "deadtilde": key = KeyEvent.VK_DEAD_TILDE; break;
            case "deadmacron": key = KeyEvent.VK_DEAD_MACRON; break;
            case "deadbreve": key = KeyEvent.VK_DEAD_BREVE; break;
            case "deadabovedot": key = KeyEvent.VK_DEAD_ABOVEDOT; break;
            case "deaddiaeresis": key = KeyEvent.VK_DEAD_DIAERESIS; break;
            case "deadabovering": key = KeyEvent.VK_DEAD_ABOVERING; break;
            case "deaddoubleacute": key = KeyEvent.VK_DEAD_DOUBLEACUTE; break;
            case "deadcaron": key = KeyEvent.VK_DEAD_CARON; break;
            case "deadcedilla": key = KeyEvent.VK_DEAD_CEDILLA; break;
            case "deadogonek": key = KeyEvent.VK_DEAD_OGONEK; break;
            case "deadiota": key = KeyEvent.VK_DEAD_IOTA; break;
            case "deadvoicedsound": key = KeyEvent.VK_DEAD_VOICED_SOUND; break;
            case "deadsemivoicedsound": key = KeyEvent.VK_DEAD_SEMIVOICED_SOUND; break;

            case "ampersand": case "amp": key = KeyEvent.VK_AMPERSAND; break;
            case "asterisk": key = KeyEvent.VK_ASTERISK; break;
            case "dblqt": case "dblquote": case "quoteDbl": key = KeyEvent.VK_QUOTEDBL; break;
            case "lt": key = KeyEvent.VK_LESS; break;case "less": key = KeyEvent.VK_LESS; break;
            case "gt": case "greater": key = KeyEvent.VK_GREATER; break;
            case "braceleft": key = KeyEvent.VK_BRACELEFT; break;
            case "braceright": key = KeyEvent.VK_BRACERIGHT; break;
            case "at": key = KeyEvent.VK_AT; break;
            case "colon": key = KeyEvent.VK_COLON; break;
            case "circumflex": key = KeyEvent.VK_CIRCUMFLEX; break;
            case "dollar": key = KeyEvent.VK_DOLLAR; break;
            case "euro": key = KeyEvent.VK_EURO_SIGN; break;
            case "exclamationmark": key = KeyEvent.VK_EXCLAMATION_MARK; break;
            case "invertedexclamationmark": key = KeyEvent.VK_INVERTED_EXCLAMATION_MARK; break;
            case "leftparenthesis": key = KeyEvent.VK_LEFT_PARENTHESIS; break;
            case "numbersign": key = KeyEvent.VK_NUMBER_SIGN; break;
            case "minus": key = KeyEvent.VK_MINUS; break;
            case "plus": key = KeyEvent.VK_PLUS; break;
            case "rightparenthesis": key = KeyEvent.VK_RIGHT_PARENTHESIS; break;
            case "underscore": key = KeyEvent.VK_UNDERSCORE; break;

            case "final": key = KeyEvent.VK_FINAL; break;
            case "convert": key = KeyEvent.VK_CONVERT; break;
            case "noconvert": key = KeyEvent.VK_NONCONVERT; break;
            case "accept": key = KeyEvent.VK_ACCEPT; break;
            case "modechange": key = KeyEvent.VK_MODECHANGE; break;
            case "kana": key = KeyEvent.VK_KANA; break;
            case "kanji": key = KeyEvent.VK_KANJI; break;
            case "alphanumeric": key = KeyEvent.VK_ALPHANUMERIC; break;
            case "katakana": key = KeyEvent.VK_KATAKANA; break;
            case "hiragana": key = KeyEvent.VK_HIRAGANA; break;
            case "fullWidth": key = KeyEvent.VK_FULL_WIDTH; break;
            case "halfWidth": key = KeyEvent.VK_HALF_WIDTH; break;
            case "romanCharacters": key = KeyEvent.VK_ROMAN_CHARACTERS; break;
            case "allCandidates": key = KeyEvent.VK_ALL_CANDIDATES; break;
            case "previousCandidate": key = KeyEvent.VK_PREVIOUS_CANDIDATE; break;
            case "codeInput": key = KeyEvent.VK_CODE_INPUT; break;
            case "japaneseKatakana": key = KeyEvent.VK_JAPANESE_KATAKANA; break;
            case "japaneseHiragana": key = KeyEvent.VK_JAPANESE_HIRAGANA; break;
            case "japaneseRoman": key = KeyEvent.VK_JAPANESE_ROMAN; break;
            case "kanaLock": key = KeyEvent.VK_KANA_LOCK; break;
            case "inputMethodOnOff": key = KeyEvent.VK_INPUT_METHOD_ON_OFF; break;

            case "again": key = KeyEvent.VK_AGAIN; break;
            case "undo": key = KeyEvent.VK_UNDO; break;
            case "copy": key = KeyEvent.VK_COPY; break;
            case "paste": key = KeyEvent.VK_PASTE; break;
            case "cut": key = KeyEvent.VK_CUT; break;
            case "find": key = KeyEvent.VK_FIND; break;
            case "props": key = KeyEvent.VK_PROPS; break;
            case "stop": key = KeyEvent.VK_STOP; break;
            default: key = 0; break;
        }
        return key;
    }

    public void keyClick(int keyCode) {
        keyPress(keyCode);
        keyRelease(keyCode);
    }

    public Component evalGuiRef(Object gui) throws GuiItemNotFoundException {
        return evalGuiRef(gui, true);
    }
    public Component evalGuiRef(Object gui, boolean required) throws GuiItemNotFoundException {
        return evalGuiRef(gui, required, false, this.GuiWaitTime);
    }
    public Component evalGuiRef(Object gui, boolean required, boolean allowNullGuiRef, int timeout) throws GuiItemNotFoundException {
        if (gui == null && allowNullGuiRef) return null;
        Objects.requireNonNull(gui);
        if (!(gui instanceof GuiItemRef))
            throw new IllegalArgumentException("evalGuiRef requires a GuiItemRef");
        long start = System.currentTimeMillis();
        GuiItemRef guiRef = (GuiItemRef) gui;
        Component found = null;
        try {
            while (System.currentTimeMillis() - start < timeout) {
                try {
                    found = _owner.guiFinder().find(guiRef, false);
                    if (found != null) break;
                    Thread.sleep(DEFAULT_POLL_TIME);
                } catch (InterruptedException ex) {
                    break;
                }
            }
            if (found == null)
                found = _owner.guiFinder().find(guiRef, required); // throw an error if required is true and we still haven't found the component.
            return found;
        } finally {
            if (_log.isTraceEnabled())
                _log.trace(found == null ? "Gui Item NOT found." : "Gui Item FOUND: " + _owner.getComponentInfo(found));
        }
    }
    /**
     * A final shutdown that cleans up some resources the Abbot Robot leaves around.
     * If this is not called, the application will hang for 60 seconds after it closes.
     */
    public static void finalShutdown() {
        try {
            Field poolField = abbot.tester.Robot.class.getDeclaredField("REAL_SYNC_POOL");
            poolField.setAccessible(true);
            ThreadPoolExecutor pool = (ThreadPoolExecutor)poolField.get(null);
            if (pool != null)
                pool.shutdown();
        } catch (ReflectiveOperationException ex) {
            ex.printStackTrace();
        }
    }
}
