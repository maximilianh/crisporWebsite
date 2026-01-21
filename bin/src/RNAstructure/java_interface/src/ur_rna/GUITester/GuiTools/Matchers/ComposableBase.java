package ur_rna.GUITester.GuiTools.Matchers;

import abbot.finder.Matcher;

import java.awt.*;
import java.util.Objects;

public abstract class ComposableBase implements ComposableMatcher {
    @Override
    public Matcher[] getChildMatchers() { return null; }

    public static final ComposableMatcher TRUE = new ComposableBase() {
        @Override
        public boolean matches(Component component) { return true; }

        @Override
        public ComposableMatcher inverse() {
            return FALSE;
        }

        @Override
        public ComposableMatcher and(Matcher other) {
            return identity(other);
        }

        @Override
        public ComposableMatcher or(Matcher other) {
            return TRUE;
        }

        @Override
        public ComposableMatcher xor(Matcher other) {
            return inverse(other);
        }

        @Override
        public ComposableMatcher simplify() {
            return this;
        }

    };
    public static final ComposableMatcher FALSE = new ComposableBase() {
        @Override
        public boolean matches(Component component) {
            return false;
        }

        @Override
        public ComposableMatcher inverse() {
            return TRUE;
        }

        @Override
        public ComposableMatcher and(Matcher other) {
            return FALSE;
        }

        @Override
        public ComposableMatcher or(Matcher other) {
            return identity(other);
        }

        @Override
        public ComposableMatcher xor(Matcher other) {
            return identity(other);
        }

        @Override
        public ComposableMatcher simplify() {
            return this;
        }

    };

    public static ComposableMatcher composable(Matcher m) {
        return identity(m);
    }

    public static ComposableMatcher inverse(Matcher operand) {
        if (operand instanceof Inverse)
            return identity(((Inverse) operand).getReference());
        if (operand instanceof Identity)
            return inverse(((Identity) operand).self());
        if (operand == TRUE)
            return FALSE;
        if (operand == FALSE || operand == null)
            return TRUE;
        return new Inverse(operand);
    }

    protected static ComposableMatcher identity(Matcher operand) {
        if (operand instanceof ComposableMatcher)
            return (ComposableMatcher) operand;
        return new Identity(operand);
    }

    protected static ComposableMatcher simplify(Matcher m) {
        if (m == null) return null;
        if (m instanceof ComposableMatcher)
            return ((ComposableMatcher) m).simplify();
        return new Identity(m);
    }

//    public static Matcher simplifyFull(Matcher m) {
//        if (m instanceof ComposableMatcher)
//            return ((ComposableMatcher) m).simplify().self();
//        return m;
//    }

    public ComposableMatcher and(Matcher other) {
        return BinaryMatcher.and(this, other);
    }

    public ComposableMatcher or(Matcher other) {
        return BinaryMatcher.or(this, other);
    }

    public ComposableMatcher xor(Matcher other) {
        return BinaryMatcher.xor(this, other);
    }

    public ComposableMatcher inverse() {
        return ComposableBase.inverse(this);
    }

    public ComposableMatcher simplify() { return this; } //override this if this matcher can be simplified or reduced.

    public Matcher self() {
        return this; //override this if there is a more direct representation of this Matcher.
    }

    public ComposableMatcher withText(String text) {
        return this.and(new DescriptionMatcher(text));
    }
    public ComposableMatcher withText(String text, int searchScope) {
        return this.and(new DescriptionMatcher(text, searchScope));
    }
    public ComposableMatcher withText(String text, boolean searchCaption, boolean searchText, boolean searchName) {
        int scope = 0;
        if (searchCaption)scope |= DescriptionMatcher.SEARCH_CAPTION;
        if (searchText)scope |= DescriptionMatcher.SEARCH_TEXT;
        if (searchName)scope |= DescriptionMatcher.SEARCH_NAME;
        return withText(text, scope);
    }

    @Override
    public String toString() {
        String type = this.getClass().getSimpleName();
        if (!type.toLowerCase().contains("matcher"))
            type = type + " Matcher";
        return "{ " + type + " }";
    }

    protected static class Inverse extends ComposableBase {
        private Matcher _m;

        public Inverse(Matcher reference) {
            _m = reference;
        }

        public Matcher getReference() {
            return _m;
        }

        @Override
        public boolean matches(Component component) {
            return !_m.matches(component);
        }

        @Override
        public ComposableMatcher simplify() {
            return inverse(simplify(_m));
        }

        @Override
        public String toString() {
            return "NOT " + Objects.toString(_m);
        }

        @Override
        public Matcher[] getChildMatchers() { return new Matcher[] { _m };  }

    }

    protected static class Identity extends ComposableBase {
        private Matcher _m;

        public Identity(Matcher reference) {
            _m = reference;
        }

        public Matcher getReference() {
            return _m;
        }

        @Override
        public boolean matches(Component component) {
            return _m.matches(component);
        }

        @Override
        public ComposableMatcher simplify() {
            if (_m instanceof ComposableMatcher)
                return ((ComposableMatcher) _m).simplify(); // implicitly performs sub-IdentityMatcher flattening.
            return this;
        }

        @Override
        public Matcher self() {
            return _m;
        }

        @Override
        public String toString() {
            return Objects.toString(_m);
        }

        @Override
        public Matcher[] getChildMatchers() { return new Matcher[] { _m };  }
    }

    protected static abstract class BinaryMatcher extends ComposableBase {
        protected BinaryMatcher(Matcher a, Matcher b) {
            if (a instanceof ComposableMatcher) a = ((ComposableMatcher) a).self();
            if (b instanceof ComposableMatcher) b = ((ComposableMatcher) b).self();
            _a = a;
            _b = b;
        }

        protected Matcher _a;
        protected Matcher _b;

        public void setRhs(Matcher m) {
            _b = m;
        }

        public Matcher getRhs() {
            return _b;
        }

        public void setLhs(Matcher m) {
            _a = m;
        }

        public Matcher getLhs() {
            return _a;
        }

        @Override
        public String toString() {
            return String.format("(%s %s %s)", _a, getOperandName().toUpperCase(), _b);
        }

        public String getOperandName() {
            return this.getClass().getSimpleName().toLowerCase();
        }

        public static ComposableMatcher and(Matcher a, Matcher b) {
            if (a == TRUE || a == null)
                return identity(b);
            if (b == TRUE || b == null)
                return identity(a);
            if (a == FALSE || b == FALSE)
                return FALSE;
            return new And(a, b);
        }

        public static ComposableMatcher or(Matcher a, Matcher b) {
            if (a == FALSE || a == null)
                return identity(b);
            if (b == FALSE || b == null)
                return identity(a);
            if (a == TRUE || b == TRUE)
                return TRUE;
            return new Or(a, b);
        }

        public static ComposableMatcher xor(Matcher a, Matcher b) {
            if (a == FALSE || a == null)
                return identity(b);
            if (b == FALSE || b == null)
                return identity(a);
            if (a == TRUE)
                return inverse(b);
            if (b == TRUE)
                return inverse(a);
            return new Xor(a, b);
        }

        @Override
        public Matcher[] getChildMatchers() { return new Matcher[] { _a, _b }; }

        public static class And extends BinaryMatcher {
            public And(Matcher a, Matcher b) {
                super(a, b);
            }

            @Override
            public boolean matches(Component c) {
                return (_a == null || _a.matches(c)) && (_b == null || _b.matches(c));
            }

            @Override
            public ComposableMatcher simplify() {
                return BinaryMatcher.and(simplify(_a), simplify(_b));
            }
        }

        public static class Or extends BinaryMatcher {
            public Or(Matcher a, Matcher b) {
                super(a, b);
            }

            @Override
            public boolean matches(Component c) {
                return (_a != null && _a.matches(c)) || (_b != null && _b.matches(c));
            }

            @Override
            public ComposableMatcher simplify() {
                return BinaryMatcher.or(simplify(_a), simplify(_b));
            }
        }

        public static class Xor extends BinaryMatcher {
            public Xor(Matcher a, Matcher b) {
                super(a, b);
            }

            @Override
            public boolean matches(Component c) {
                return (_a != null && _a.matches(c)) ^ (_b != null && _b.matches(c));
            }

            @Override
            public ComposableMatcher simplify() {
                return BinaryMatcher.xor(simplify(_a), simplify(_b));
            }
        }
    }
}
