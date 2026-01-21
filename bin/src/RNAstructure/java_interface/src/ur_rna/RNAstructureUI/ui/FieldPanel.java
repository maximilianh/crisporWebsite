/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.RNAstructureUI.ui;

import ur_rna.RNAstructureUI.RNAstructure;
import ur_rna.RNAstructureUI.utilities.FileSelectedEvent;
import ur_rna.RNAstructureUI.utilities.MRUFileStorage;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import static ur_rna.Utilities.Strings.isEmpty;

/**
 * A class that creates a panel which holds aligned input fields.
 *
 * @author Jessica S. Reuter
 * @author Richard Watson
 */
public class FieldPanel
        extends JPanel {
    private static final long serialVersionUID = 20120802;

    /**
     * The specified width of the panel.
     */
    private int specifiedPanelWidth;

    /**
     * The array of text fields.
     */
    protected ArrayList<JTextField> textFields = new ArrayList<>();

    protected FieldPanel() {
        setLayout(new GridBagLayout());
        // Set a border around the panel.
        int gap = 10;
        setBorder(BorderFactory.createEmptyBorder(gap, gap, gap, gap));
    }
    /**
     * Constructor.
     *
     * @param fields The array of text fields in this panel.
     */
    public FieldPanel(JTextField... fields) {
        this();
        for (JTextField f : fields)
            addField(f);
    }
    public void addField(final JTextField field) { addField(field, null, null, null); }
    public void addField(final JTextField field, String labelText) { addField(field, labelText, null, null); }
    public void addField(final JTextField field, String labelText, String actionCommand) { addField(field, labelText, actionCommand, null); }
    public void addField(final JTextField field, String labelText, String actionCommand, ActionListener listener) {
        if (labelText == null) labelText = field.getName();
        if (actionCommand == null) actionCommand = field.getName();
        GridBagConstraints gbc = new GridBagConstraints();
        textFields.add(field);
        addFieldCore(field, labelText, actionCommand, listener, gbc);
    }

    protected void addFieldCore(final JTextField f, final String text, final String cmd, final ActionListener al, final GridBagConstraints gbc) {
        JLabel label = new JLabel(text);
        label.setHorizontalAlignment(JLabel.TRAILING);
        //new BorderBuilder().makeRightBorder(2, label);
        label.setLabelFor(f);

        gbc.gridy = textFields.size();

        gbc.gridx = 0;
        gbc.insets = new Insets(1, 2, 1, 2);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.LINE_START;
        add(label, gbc);
        gbc.gridx = 1;
        gbc.gridwidth = 2;
        gbc.weightx = 0.8;

        add(f, gbc);
    }

    /**
     * Get a  particular input field.
     *
     * @param index The text field to get, one-indexed.
     * @return The input field.
     */
    public JTextField getField(int index) {
        return textFields.get(index - 1);
    }
    public JTextField getField(String name) {
        for (JTextField field : textFields) {
            if (name.equals(field.getName()))
                return field;
        }
        return null;
    }

    /**
     * Build the panel.
     */
    public void makePanel() {
//		// Determine the border that should go around the panel.
//		int border = 10;
//
//		// Determine the width of the widest component associated with the
//		// fields, so all components can be that width and align.
//		int maxWidth = 0;
//		int numFields = textFields.length;
//		for( int i = 1; i <= numFields; i++ ) {
//			int nextWidth = components[i-1].getPreferredSize().width;
//			maxWidth = Math.max( maxWidth, nextWidth );
//		}
//
//		// For each component in the array, place it next to its field.
//		for( int i = 1; i <= numFields; i++ ) {
//
//			// Set the preferred width of the component associated with the
//			// data field.
//			JComponent component = components[i-1];
//			int height = component.getPreferredSize().height;
//			component.setPreferredSize( new Dimension( maxWidth, height ) );
//
//			// Create the panel that holds the new component and text field,
//			// then add it to the input panel.
//			JPanel panel = new JPanel( new BorderLayout() );
//			panel.add( component, BorderLayout.WEST );
//			panel.add( textFields[i-1] );
//			add( panel );
//		}
//
        // Set the full panel size.
//		int preferredHeight =
//			( new JTextField().getPreferredSize().height * numFields ) +
//			( 2 * border );
        Dimension ps = this.getPreferredSize();
        int h = (int) ps.getHeight();
        int w = (int) Math.max(specifiedPanelWidth, ps.getWidth());
        Dimension size = new Dimension(w, h);
        setPreferredSize(size);
        //setMinimumSize( preferredDimension );
        //setMaximumSize( preferredDimension );
        setSize(size);
    }

    /**
     * Set the width of the panel.
     */
    public void setPanelWidth(int width) {
        specifiedPanelWidth = width;
    }

    /**
     * An inner class that handles a panel which contains file fields.
     *
     * @author Jessica S. Reuter
     */
    public static class FilePanel
            extends FieldPanel {
        private static final long serialVersionUID = 20120802;
        private final ActionListener buttonListener;
        private static MRUFileStorage recentFiles = new MRUFileStorage();

        /**
         * Constructor.
         * <br><br>
         * Create a panel of file fields attached to buttons.
         *
         * @param listener The action listener attached to the buttons.
         * @param fields   The file fields.
         */
        public FilePanel(final ActionListener listener, final FileField... fields) {
            super(fields);
            buttonListener = listener;
        }


        public void saveRecent() {
            MRUFileStorage recent = RNAstructure.MRUFiles;
            for (JTextField f : textFields) {
                if (f instanceof FileField) {
                    String file = ((FileField) f).getFile();
                    if (!isEmpty(file))
                        recent.add(file);
                }
            }
            recent.saveToStorage();
        }

        private class ButtonListener implements ActionListener {
            final int index;
            public ButtonListener(int fieldIndex) {
                index = fieldIndex;
            }
            public void actionPerformed(ActionEvent e) {
                // If the button is not the first, check to make sure all previous fields have been filled.
                if (!getField(index).isEnabled()) {
                    for (int i = 1; i <= index - 1; i++) {
                        if (getFile(i).equals("")) {
                            String message = getField(index).getName() + " cannot be set until after all preceding files (e.g. " + getField(i).getName() + ").";
                            Dialogs.showError(message);
                            return;
                        }
                    }
                }

                // The event was fired by the Examples or Recent Files menus. So notify the StandardFileChooser that the request for a file name (which will be triggered by this event) should be
                // ignored, and the selected file should be returned.
                if (e instanceof FileSelectedEvent)
                    StandardFileChooser.preemptNextFileNameRequest(((FileSelectedEvent) e).getFile().toString(), ((FileField)getField(index)).getFilters());

                buttonListener.actionPerformed(e);
            }
        }

        /**
         * Get a file name from a particular input field.
         *
         * @param index The text field whose contents to get, one-indexed.
         * @return The file name.
         */
        public String getFile(int index) { return getField(index).getText().trim(); }

        /**
         * Get whether any fields are empty; an empty field is an error.
         *
         * @return Whether an error happened.
         */
        public boolean isError() {
            for (JTextField element : textFields) {
                if (element.getText().trim().equals("")) {
                    String msg = element.getName() + " must not be empty.";
                    Dialogs.showError(msg);
                    return true;
                }
            }
            return false;
        }

        /**
         * Set a file name in a particular input field.
         *
         * @param index The text field whose contents to set, one-indexed.
         * @param file  The file name to set in the field.
         */
        public void setFile(int index, String file) { getField(index).setText(file); }

        @Override
        protected void addFieldCore(final JTextField f, final String text, final String cmd, final ActionListener al, final GridBagConstraints gbc) {
            super.addFieldCore(f, text, cmd, al, gbc);
            int index = textFields.size();
            gbc.gridx += 2;
            gbc.gridwidth = 1;
            gbc.weightx = 0;
            JButton button = new JButton("...");
            button.setName(f.getName());
            button.setActionCommand(cmd);
            button.setToolTipText("Browse to select file.");
            // Add an action listener to the button.
            ActionListener bl = new ButtonListener(index);
            button.addActionListener(bl);
            add(button, gbc);
            if (f instanceof FileField) {
                FileField ff = (FileField)f;
                if (ff.isInputFile() && ff.getFilters() != null) {
                    RecentFileButton b = new RecentFileButton("Recent", ff.getFilters(), bl, button.getActionCommand());
                    gbc.gridx+=1;
                    add(b, gbc);
                }
            }
        }
    }

    //    /**
//     * An inner class that handles a panel which contains numeric data fields.
//     *
//     * @author Jessica S. Reuter
//     */
//    public static class FieldPanel
//            extends FieldPanel {
//        private static final long serialVersionUID = 20120802;
//
//        /**
//         * Constructor.
//         * <br><br>
//         * Create a panel of numeric data fields attached to labels.
//         *
//         * @param fields The data fields.
//         */
//        public FieldPanel(NumberField... fields) { super(fields); }
//    }
}
