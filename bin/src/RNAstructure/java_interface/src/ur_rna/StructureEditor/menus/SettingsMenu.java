package ur_rna.StructureEditor.menus;

import ur_rna.StructureEditor.Program;
import ur_rna.StructureEditor.Settings;
import ur_rna.StructureEditor.services.RnaDrawController;
import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.swing.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;

import static ur_rna.Utilities.Strings.isEmpty;

public class SettingsMenu extends MainMenu {
    //private CheckMergeItem dashboard, drawClockwise, drawHints;
    private Settings settings = Program.getInstance().settings();
    private boolean updatingSettings;

    //public final UiAction actionDrawClockwise = new UiAction("Draw Structures Clock&wise", ()->{}).setDesc("When (re)drawing the structure, it can be drawn (from 5' to 3') either clockwise or counter-clockwise. This setting only has an effect when you click \"Redraw\".");

    public SettingsMenu() {
        super("&Settings");
        addCheck("Enable &Dashboard", "ShowDashboard","Show the Dashboard (aka Start-Up Screen) when no other files are loaded.");
        addCheck("Draw Structures Clock&wise", "DrawClockwise", "<html>When (re)drawing the structure, it can be drawn either clockwise or counter-clockwise.<br>This setting only has an effect when you click <b>\"Redraw\"</b>.");
        addCheck("Show Visual Queues", "ShowHints", "Show temporary visual guides (lines, arcs, etc) while dragging nucleotides during rotations, loop resizing, etc.");
        addCheck("Enable \"Sticky\" Positions", "SnapGuides", "<html>When dragging dragging nucleotides (during rotations, loop resizing, etc.) the mouse will <i>\"stick\"</i> or <i>\"snap\"</i> to special angles or positions.<br>" +
                        "(This can be temporarily disabled by holding down the <b>" + AcceleratorKey.getModifierKeyName(RnaDrawController.DRAG_SPECIAL_MASK) + "</b> Key.");
        addCheck("Resize MultiLoops by Default", "ExpandMotif", "<html>Usually only one loop segment is affected when resizing loops. But if this setting is enabled, the full loop circuit is expanded.<br>Hold down the <b>" + AcceleratorKey.getModifierKeyName(RnaDrawController.DRAG_EXPANDED_MOTIF_MASK) + "</b> Key to temporarily toggle this setting.");

        final MergeMenu tools = new MergeMenu("Disable Layout Tools (aka \"Drag Handles\")");

        tools.addCheck(new UiAction("Disable Free &Rotation (and the Center-of_Rotation cross-hair", this::checkSettingChanged)//, Program.getIcon("rotate-center"))
                .setDesc("Check this setting if the free-rotation feature is distracting or in the way (i.e. the center-of-rotation cross-hair icon or rotation icon).")
                .setCommand("DisableRotation"));
        tools.addCheck(new UiAction("Disable &Loop Resizing", this::checkSettingChanged)//, Program.getIcon("loop-resize-handle"))
                .setCommand("DisableLoopResize"));
        tools.addCheck(new UiAction("Disable &Branch Sliding", this::checkSettingChanged)//, Program.getIcon("slide-branch-handle"))
                .setCommand("DisableBranchSlide"));
        tools.addSeparator();
        tools.add(new UiAction("Enable &All", e->checkAll(tools, false)));
        tools.add(new UiAction("&Disable All", e->checkAll(tools, true)));

        add(tools);
        //add(test=new CheckMergeItem("Test Option", this::checkSettingChanged, "Test")).setToolTipText("This is a test");
        settings.addChangeHandler(this::settingsChanged);
        loadSettings();
    }

    private void checkAll(Container parent, boolean enabled) {
        for(Component c : Menus.getMenuItems(parent))
            if (c instanceof AbstractButton)
                ((AbstractButton) c).setSelected(enabled);
        saveSettings();
    }

    private CheckMergeItem addCheck(String text, String settingName, String tooltip) {
        CheckMergeItem mi = new CheckMergeItem(text, this::checkSettingChanged, settingName);
        mi.setToolTipText(tooltip); add(mi); return mi;
    }

    private void settingsChanged() {
        updatingSettings = true;
        loadSettings();
        updatingSettings = false;
    }
    private void loadSettings() {
        loadSettings(this);
    }
    private void loadSettings(Component c) {
        for (JMenuItem item : Menus.getMenuItems(c)) {
            if (item instanceof JCheckBoxMenuItem && !isEmpty(item.getActionCommand()))
                item.setSelected(getBoolSetting(item.getActionCommand()));
            if (item instanceof JMenu)
                loadSettings(item);
        }
    }
    // This only needs to be called when multiple settings are changed en-mass
    private void saveSettings() {
        settings.suspendChangeNotification();
        try {
            saveSettings(this);
        } finally {
            settings.resumeChangeNotification();
        }
    }
    private void saveSettings(Component c) {
        String name;
        for(JMenuItem mi : Menus.getMenuItems(c)) {
            if (mi instanceof JCheckBoxMenuItem) {
                if (!isEmpty(name=mi.getActionCommand()))
                    settings.set(name, mi.isSelected());
            } else if (mi instanceof JMenu)
                saveSettings(mi); // save sub-menu-items
        }
    }
    private boolean getBoolSetting(final String name) {
        return ObjTools.asBool(settings.get(name));
    }

    private void checkSettingChanged(final ActionEvent e) {
        if (updatingSettings) return;
        String name = e.getActionCommand();
        boolean value = ((AbstractButton)e.getSource()).isSelected();
        settings.set(name, value);
    }
}
