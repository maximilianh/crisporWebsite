package ur_rna.RNAstructureUI.utilities;

import java.awt.event.ActionEvent;
import java.io.File;

/**
 * @author Richard M. Watson
 */
public class FileSelectedEvent extends ActionEvent {
    private final File file;
    public FileSelectedEvent(Object source, File file, String command) {
        super(source, ActionEvent.ACTION_PERFORMED, command, 0, 0);
        this.file = file;
    }
    public File getFile() { return file; }
}
