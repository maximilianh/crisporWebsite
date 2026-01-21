package ur_rna.GUITester;

import ur_rna.Utilities.Strings;

import java.io.Reader;
import java.io.StringReader;

class Script {
    public final String filePath;
    public final long startLine, endLine;
    public final Reader input;
    public String title;

    public Script(String filePath, String scriptText) {
        this(filePath, new StringReader(scriptText), 0, -1);
    }
    public Script(String filePath, Reader input) {
        this(filePath, input, 0, -1);
    }
    public Script(String filePath, Reader input, long start, long end) {
        this.filePath = filePath;
        this.input = input;
        this.startLine = start;
        this.endLine = end;
    }
    public String getFileInfo(boolean quotePath) {
        final String CMD = "(Command Line)";
        if (filePath == null)
            return CMD;
        return quotePath ? "\"" + filePath + "\"" : filePath;
    }
    public String getSegmentInfo() {
        return isSegment() ? "[lines " + startLine + " to " + endLine + "]" : "";
    }
    public boolean fromFile() {
        return filePath != null;
    }
    public boolean isSegment() {
        return startLine != 0 || endLine >= 0;
    }
    public Reader getInput() { return input; }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Script ");
        if (!Strings.isEmpty(title)) {
            sb.append('\"');
            sb.append(title);
            sb.append('\"');
            sb.append(' ');
        }
        sb.append('{');
        sb.append("file: ");
        sb.append(getFileInfo(true));
        if (isSegment()) {
            sb.append(' ');
            sb.append(getSegmentInfo());
            sb.append(' ');
        }
        sb.append('}');
        return sb.toString();
    }
}
