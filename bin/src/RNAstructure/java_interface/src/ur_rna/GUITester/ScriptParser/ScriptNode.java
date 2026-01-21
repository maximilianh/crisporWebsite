package ur_rna.GUITester.ScriptParser;

import ur_rna.Utilities.annotation.NotNull;

import java.io.PrintStream;

public class ScriptNode extends SimpleNode implements ScriptNodePropConstants {
  protected static final Node[] emptyChildArray = new Node[0];
  public int tokenKind;
  // public Token setToken;
  // public Node setNode;
  // protected List<Prop> props;
  
  public ScriptNode(int i) { super(i); }
  public ScriptNode(GuiTestScriptParser p, int i) { super(p, i); }
  
//  /** gets tmpToken and then sets it to null. */
  // public Token getToken() {
	  // Token t = setToken;
	  // setToken=null;
	  // return t;
  // }
  // /** gets tmpToken and then sets it to null. */
  // public Node tmpNode() {
	  // Node n = setNode;
	  // setNode=null;
	  // return n;
  // }
  
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("[");
    sb.append(getNodeName());
    writeProps(sb);
    sb.append("]");
	return sb.toString();
  }

  public void writeProps(StringBuilder sb) {
    if (value != null) sb.append(" Value: " + value);
    if (tokenKind != 0) sb.append(" Token: " + getTokenName());
  }

  public String propsToString() {
    StringBuilder sb = new StringBuilder();
    writeProps(sb);
    return sb.toString();
  }

  public ScriptNode child(int i) {
    return (ScriptNode)children[i];
  }

  private class DumpState {
    /*
      TODO: Make a nice tree dump
      ??Branch
      ??Branch
      ?  ??Limb
      ?  ??Limb
      ?  ?  ??Leaf
      ?  ?  ??Leaf
      ?  ??Limb
      ??Branch

     */
    // ? ?? ?
    public DumpState(PrintStream outputStream, String initialPrefix) {
      this.outputStream = outputStream;
      this.prefix = initialPrefix;
    }
    private final PrintStream outputStream;
    public String descendant_line = "\u2503"; // vertical line:   ┃
    public String branch_mid_line = "\u251d"; // sideways T:      ┝
    public String branch_end_line = "\u2515"; // short L-shape:   ┕
    public String branch_indent = "\u2500";   // wide horiz line: ─
    public String nonchild_indent = "  ";
    private String prefix = "";

    int level;
    StringBuilder sb;

    public void enterScope() {
      pushLevel(true);
    }
    public void exitScope() {
      pushLevel(false);
    }

    public String getCurrentNodePrefix(boolean isLastChild) {
      return prefix + (isLastChild ? branch_end_line : branch_mid_line) + branch_indent;
    }
    public void addNode(String description) {
      addNode(description, 1); //assume it's a middle-node
    }
    // positionFromEnd allows us to determine whether the node is:
    //   <= -1: The root node (use no branch prefix)
    //       0: The final node in the list of its siblings. (i.e. use branch_end_line)
    //   >= +1: A middle node in the list of its siblings. (i.e. use branch_mid_line)
    public void addNode(String description, int positionFromEnd) {
      print(prefix);
      if (positionFromEnd >= 0) {
        print(positionFromEnd==0 ? branch_end_line : branch_mid_line);
        print(branch_indent);
      }
      println(description);
    }
    public void print(String text) {
      outputStream.print(text);
    }
    public void println(String text) {
      outputStream.println(text);
    }
    public void println(String s1, String s2, String s3) {
      if (s1 != null) outputStream.print(s1);
      if (s2 != null) outputStream.print(s2);
      if (s3 != null) outputStream.print(s3);
      outputStream.println();
    }

    private void pushLevel(boolean push) {
      if (sb == null) { sb = new StringBuilder(); sb.append(prefix); }
      if (push) {
        level++;
        if (level > 1)
          sb.append(descendant_line+nonchild_indent);
      } else {
        level--;
        if (level > 0)
          //remove the last indent.
          sb.setLength(sb.length()-descendant_line.length()-nonchild_indent.length());
      }
      prefix = sb.toString();
    }
  }

  public String getTokenName() { return tokenKind == 0 ? "" : getTokenName(tokenKind); }
  public static String getTokenName(int kind) { return GuiTestScriptParserConstants.tokenImage[kind]; }
  public String getNodeName() { return getNodeName(this.id);}
  public static String getNodeName(int nodeId) { return GuiTestScriptParserTreeConstants.jjtNodeName[nodeId]; }


  @Override
  public void dump(String prefix) {
    dump(new DumpState(System.out, prefix), -1); //dump as if it is the root node
  }
  public void dump(String prefix, PrintStream out) {
    dump(new DumpState(out, prefix), -1); //dump as if it is the root node
  }
  private int getPositionFromEnd() {
    if (this.parent == null)
      return -1; //root
    else {
      Node[] siblings = ((ScriptNode)this.parent).getChildren();
      if (siblings == null) return -1; //should never be the case because if this has a parent, it should be in the children array.
      int pos = java.util.Arrays.asList(siblings).indexOf(this);
      return siblings.length - pos - 1;
    }
  }

  private void dump(DumpState state) { dump(state, getPositionFromEnd()); }
  private void dump(DumpState state, int nodePositionFromEnd) {
    state.addNode(this.toString(), nodePositionFromEnd);
    if (children != null) {
      state.enterScope();
      int len = children.length;
      for (int i = 0; i < len; i++) {
        Node n = children[i];
        if (n instanceof ScriptNode) {
          ((ScriptNode)n).dump(state, len-i-1);
        } else if (n != null) {
          ((SimpleNode)n).dump(state.prefix);
		}
      }
      state.exitScope();
    }
  }

  public void add(Node node) { this.jjtAddChild(node, children.length); }
  public ScriptNode get(int pos) { return (ScriptNode)this.jjtGetChild(pos); }

  public void setTokenKind(int kind) { this.tokenKind = kind; }
  public Object getTokenKind() { return tokenKind; }

  public int childCount() {    return children == null ? 0 : children.length;  }
  public boolean hasChildren() {
    return children != null && children.length > 0;
  }

  @NotNull
  public Node[] getChildren() {
    if (children == null) return ScriptNode.emptyChildArray;
    return children;
  }

  // public void setProp(int propID, Object propValue) {
	  // Prop found = findProp(propID)
	  // if (found == null) {
		  // if (props == null) props = new List<>();
		  // props.add(new Prop(propID, propValue));
	  // } else
		  // found.value = propValue;
  // }
  
  // public Object getProp(int propID) {
	  // Prop found = findProp(propID)
	  // if (found == null) return null;
	  // return found.value;
  // }
  
  // public Prop findProp(int propID) {
	  // if (props != null)
		  // for (Prop p : props)
			  // if (p.id == propID)
				  // return p;
	  // return null;
  // }
  
  //public void setToken(Token value) { this.t = value; }
  //public Token getToken() { return t; }  
  
  // public setLhs(Node value) { this.value = value; }
  // public Node getLhs() { return value; }
  
  // public setRhs(Node value) { this.value = value; }
  // public Node getRhs() { return value; }
}