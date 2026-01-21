package ur_rna.GUITester.GuiTools.Matchers;

 import abbot.finder.Matcher;

 public interface ComposableMatcher extends Matcher {
     ComposableMatcher and(Matcher other);
     ComposableMatcher or(Matcher other);
     ComposableMatcher xor(Matcher other);
     ComposableMatcher inverse();
     ComposableMatcher simplify();
     Matcher self();

     /**
      * Return an array of all Matchers that are immediate descendants of this one, or an empty array if
      * there are none.
      * @return A non-null array of Matchers.
      */
     Matcher[] getChildMatchers();
}
