package ur_rna.Utilities;

import ur_rna.Utilities.annotation.NotNull;

import java.util.*;

public class CommandLineParser {
    // private final ArrayList<Argument> args = new ArrayList<>();
    // private final ArrayList<FlagDef> flagDefMap = new ArrayList<>();
    private final HashMap<String, FlagDef> flagDefMap = new HashMap<>();
    public final boolean ignoreFlagCase;

    public boolean defaultAllowRepeatFlags; // default false;
    public ParamAffinity defaultParamAffinity; //default ParamAffinity.FORBIDDEN;
    /**
     * Whether or not to allow undefined flags when parsing the command line.
     * If this is true and an undefined flag is encountered, a new FlagDef will automatically be created. The new
     * FlagDef will have the found flag as its {@link FlagDef#name name}
     * value and it will have no {@link FlagDef#aliases aliases}. The values for
     * {@link FlagDef#allowRepeats allowRepeats} and {@link FlagDef#paramAffinity paramAffinity} will be
     * {@link #defaultAllowRepeatFlags} and {@link #defaultParamAffinity} respectively.
     */
    public boolean allowUndefinedFlags;

    //public boolean allowUndefinedFlags() { return allowUndefinedFlags; }
    //public void allowUndefinedFlags(final boolean allowUndefinedFlags) { this.allowUndefinedFlags = allowUndefinedFlags; }

    public CommandLineParser() { this(false); }
    public CommandLineParser(boolean ignoreFlagCase) { this(ignoreFlagCase, ParamAffinity.FORBIDDEN, false); }
    public CommandLineParser(boolean ignoreFlagCase, ParamAffinity defaultParamAffinity, boolean defaultAllowRepeatFlags) {
        this.ignoreFlagCase = ignoreFlagCase;
        this.defaultParamAffinity = defaultParamAffinity;
        this.defaultAllowRepeatFlags = defaultAllowRepeatFlags;
    }

    /**
     * Add a flag and specify whether or not it requires a following parameter.
     *
     * @param flagName          The name of the flag. If <tt>ignoreFlagCase</tt>, it will be converted to lower case.
     * @param requiresParameter Whether or not the flag requires a parameter. (The parameter is never optional.)
     */
    public FlagDef addFlag(String flagName, boolean requiresParameter) {
        ParamAffinity defaultAffinity = defaultParamAffinity == ParamAffinity.FORBIDDEN ?
                ParamAffinity.FORBIDDEN :
                ParamAffinity.OPTIONAL;
        return addFlag(flagName, requiresParameter ? ParamAffinity.REQUIRED : defaultAffinity, defaultAllowRepeatFlags);
    }
    public FlagDef addFlag(String... aliases) {
        return addFlag(aliases[0], defaultParamAffinity, defaultAllowRepeatFlags, Arrays.copyOfRange(aliases, 1, aliases.length));
    }
    public FlagDef addFlag(String name, ParamAffinity canHaveParams, boolean allowRepeats, String... aliases) {
        FlagDef f = new FlagDef(name, canHaveParams, allowRepeats, aliases);
        addFlag(f);
        return f;
    }

    public void addFlag(FlagDef f) {
        if (f.name != null) addFlagAlias(f, f.name);
        if (f.aliases != null)
            for (String alias : f.aliases)
                addFlagAlias(f, alias);
    }

    public void addFlagAlias(final String existingName, final String alias) {
        FlagDef f = findFlagDef(existingName, true);
        if (existingName == null)
            throw new IllegalArgumentException(Strings.fmt("No flag has been defined with the name '%s'.", existingName));
        addFlagAlias(f, alias);
    }
    public void addFlagAlias(final FlagDef f, final String alias) {
        FlagDef other = flagDefMap.putIfAbsent(alias, f);
        if (other != null && other != f)
            throw new KeyAlreadyExistsException(Strings.fmt("The alias '%s' is already associated with the %s flag.", alias, other.name))
                    .initDetails(alias, other);
    }

    /**
     * Define flags that do NOT have parameters.
     *
     * @param flagNames An array or list of flag names.
     */
    public void addFlags(String... flagNames) {
        for (String flagName : flagNames)
            addFlag(flagName, false);
    }

    /**
     * Define flags that have (and require) parameters.
     *
     * @param flagNames An array or list of flag names.
     */
    public void addFlagsWithParams(String... flagNames) {
        for (String flagName : flagNames)
            addFlag(flagName, true);
    }


    /**
     * Finds the FlagDef associated with the given name.
     *
     * @param nameOrAlias The name or alias to search for. This can be either the
     *                    {@link CommandLineParser.FlagDef#name name} of the FlagDef or one of its
     *                    {@link CommandLineParser.FlagDef#aliases aliases}.
     * @return The FlagDef associated with the given name or alias, or null if no matching FlagDef has been defined.
     */
    public FlagDef findFlagDef(String nameOrAlias) {
        return flagDefMap.get(normalizeFlagName(nameOrAlias));
    }

    /**
     * /**
     * Finds the FlagDef associated with the given name.
     *
     * @param nameOrAlias     The name or alias to search for. This can be either the
     *                        {@link CommandLineParser.FlagDef#name name} of the FlagDef or one of its
     *                        {@link CommandLineParser.FlagDef#aliases aliases}.
     * @param required        If true and the FlagDef is not defined, an exception will be thrown, unless
     *                        {@link #allowUndefinedFlags } is true, in which chase the exception will not be thrown.
     * @return The FlagDef associated with the given name or alias, or null if no matching FlagDef has been defined.
     * @throws IllegalArgumentException if the FlagDef is not defined and {@code required} is true and
     *                                  {@link #allowUndefinedFlags } is false.
     */
    public FlagDef findFlagDef(String nameOrAlias, boolean required) {
        FlagDef f = flagDefMap.get(normalizeFlagName(nameOrAlias));
        if (f == null && required && !allowUndefinedFlags)
                throw new IllegalArgumentException("The flag '" + nameOrAlias + "' has not been defined.");
        return f;
    }

    public Set<FlagDef> getFlagDefinitions() {
        return new HashSet<>(flagDefMap.values());
    }

    public List<FlagDef> parseFlagDefs(String flagDefSpecifier)
            throws SyntaxErrorException { return parseFlagDefs(flagDefSpecifier, true); }
    public List<FlagDef> parseFlagDefs(String flagDefSpecifier, boolean addResults) throws SyntaxErrorException {
        //  verbose,v:R  file,f:R*  pattern,p:A*  trace,t:D

        String[] flags = flagDefSpecifier.split(" +"); //split on whitespace
        if (flags.length == 0) return Collections.emptyList();
        List<FlagDef> results = new ArrayList<>(flags.length);

        for (String flagInfo : flags) {
            String aliasList, options = "";
            ParamAffinity params = defaultParamAffinity;
            boolean multiple = defaultAllowRepeatFlags;

            int pos = flagInfo.indexOf(':');
            if (pos == -1)
                aliasList = flagInfo;
            else {
                aliasList = flagInfo.substring(0, pos);
                if (++pos < flagInfo.length())
                    options = flagInfo.substring(pos);

                for (char opt : options.toUpperCase().toCharArray()) {
                    switch (opt) {
                        case '?':
                        case 'O': params = ParamAffinity.OPTIONAL; break;
                        case '!':
                        case 'R': params = ParamAffinity.REQUIRED; break;
                        case 'X':
                        case '/':
                        case 'F': params = ParamAffinity.FORBIDDEN; break;
                        case '*':
                        case '+':
                        case 'M': multiple = true; break;
                        case '1': multiple = false; break;
                        default:
                            throw new SyntaxErrorException("Invalid flag definition option.", Character.toString(opt), "in the definition of " + aliasList);
                    }
                }
            }
            String[] aliases = aliasList.split(",");
            if (aliases.length == 0 || aliases[0].trim().length() == 0)
                throw new SyntaxErrorException("No aliases specified for flag.", flagInfo, "in the definition of flag " + results.size() + 1);
            String name = aliases[0];
            aliases = Arrays.copyOfRange(aliases, 1, aliases.length);
            results.add(new FlagDef(name, params, multiple, aliases));
        }
        if (addResults)
            for (FlagDef fd : results)
                addFlag(fd);
        return results;
    }

    public ParseResults parse(String commandLine) throws SyntaxErrorException {
        return parse(CommandLineSplitter.split(commandLine));
    }
    public ParseResults parse(String[] commandArgs) throws SyntaxErrorException {
        List<Argument> args = new ArrayList<>(commandArgs.length);
        int argCount = 0;
        for (int i = 0; i < commandArgs.length; i++) {
            String fullText;
            String arg = fullText = commandArgs[i];
            if (arg.startsWith("-")) {
                Pair<String, String> fv = splitFlagValue(removeFlagPrefix(arg));
                FlagDef fd = findFlagDef(fv.key);

                if (fd == null) {
                    if (allowUndefinedFlags)
                        fd = addFlag(fv.key, defaultParamAffinity, defaultAllowRepeatFlags);
                    else
                        throw new SyntaxErrorException("Invalid flag or option.", fv.key, "at argument " + (i+1));
                }

                String value = fv.value;
                if (fd.paramAffinity == ParamAffinity.REQUIRED) {
                    if (value == null) {
                        if (i + 1 >= commandArgs.length)
                            throw new SyntaxErrorException(String.format("The flag '%s' requires a parameter.", fd.name), arg, "at the end of the argument list.");
                        value = commandArgs[++i];
                        fullText += " " + value;
                    }
                } else if (fd.paramAffinity == ParamAffinity.FORBIDDEN) {
                    if (value != null)
                        throw new SyntaxErrorException(String.format("The flag '%s' does not accept a parameter value.", fd.name), arg, "at argument " + (i+1));
                }
                args.add(new FlagArg(fullText, argCount++, fd, fv.key, value, value));
            } else
                args.add(new UnnamedArg(fullText, argCount++));
        }
        return new ParseResults(this, args);
    }

    private Pair<String, String> splitFlagValue(String flagNameAndValue) {
        int pos1 = flagNameAndValue.indexOf('=');
        int pos2 = flagNameAndValue.indexOf(':');
        int minPos = Math.min(pos1, pos2);
        if (minPos == -1)
            return new Pair<String, String>(flagNameAndValue, null);
        return new Pair<String, String>(
                flagNameAndValue.substring(0, minPos),
                flagNameAndValue.substring(minPos + 1, flagNameAndValue.length()));
    }

    public String normalizeFlagName(String flagName) {
        String s = removeFlagPrefix(flagName);
        return ignoreFlagCase ? s.toLowerCase() : s;
    }
    public String removeFlagPrefix(String s) {
        return s.startsWith("-") ? s.substring(1) : s;
    }
    public String addFlagPrefix(String s) {
        return s.startsWith("-") ? s : "-" + s;
    }
    public String getUsageMessage() {
        StringBuilderEx sb = new StringBuilderEx();
        Set<FlagDef> fds = getFlagDefinitions();
        Set<String> names = new HashSet<>(4);
        for (FlagDef f : fds) {
            names.clear();
            names.add(f.name);
            names.addAll(Arrays.asList(f.aliases));
            // get any additional names added after the FlagDef was created.
            for (Map.Entry<String,FlagDef> e : flagDefMap.entrySet())
                if (e.getValue() == f) names.add(e.getKey());
            sb.append("Flag Names: ");
            for (String name : names) {
                sb.append(addFlagPrefix(name));
                sb.append(", ");
            }
            sb.truncate(2); //remove final ", "
            sb.append(";\tParameter: " + f.paramAffinity.toString());
            sb.append(";\tRepeats: " + (f.allowRepeats ? "yes" : "no"));
            sb.append('\n');
        }
        return sb.toString();
    }

/* ------------------------------------------------------------------------------ */
/* ---------------------- Begin inner/nested classes ---------------------------- */
/* ------------------------------------------------------------------------------ */

    /**
     * Describes whether a parameter value is required, optional (allowed),
     * or forbidden (invalid) for a given flag.
     */
    public enum ParamAffinity {
        /** The flag requires a parameter. */
        REQUIRED,
        /** The flag can have a parameter, but it is optional. */
        OPTIONAL,
        /** Parameters are not allowed for this flag. */
        FORBIDDEN
    }

    /**
     * Arguments are the results of parsing a commandline.
     * There are two types of arguments: flags (aka <em>named</em> arguments), represented by {@link FlagArg } objects
     * and <em>unnamed</em> arguments, represented by  {@link UnnamedArg } objects.
     */
    public static abstract class Argument {
        /** The full original text of this Argument in the command line. */
        public final String fullText;

        /** The position of this Argument in the parsed command line. */
        public final int position;

        protected Argument(String text, int position) {
            this.fullText = text;
            this.position = position;
        }

//        public abstract boolean isFlag();
//        public abstract FlagDef associatedFlag();
//        public abstract String getValue();
    }

    /**
     * Flags arguments are indicated by a preceding hyphen (-) or slash (/).
     * Flag arguments have a {@link FlagArg#flag flag} field which indicates the {@link FlagDef flag definition}
     * that provided the rules for parsing and producing this FlagArg object.
     * Some flags have parameters (which can be optional, required, or forbidden/not allowed).
     * If a parameter value is encountered during parsing, it is stored in the
     * {@link FlagArg#valueText valueText} field, which is {@code null} in Flags for which no parameter value
     * was encountered.
     */
    public static class FlagArg extends Argument {
        /**
         * Gets the actual flag name from the command-line, without any parameter value or separator that
         * might have been present.
         * <p>
         * When a FlagArg is parsed, the flag name (e.g. "-f" ) is used to associate it with a FlagDef.
         * From that point on, the FlagArg can be retrieved using the FlagDef name (e.g. "file")
         * or any of its aliases.
         * Usually the actual flag name (-f vs -file) is not relevant, but the name is stored in case some
         * implementation requires it.
         * </p>
         */
        public final String flagText;

        /**
         * The Flag associated with this Argument, which specifies the name (and aliases) of the
         * flag as well as the rules used to parse it.
         */
        public final FlagDef flag;
        /** Gets the value of this flag, if one was specified, or null otherwise. */
        public final String valueText;
        /**
         * Gets the value of this flag once it has been converted
         * to the type specified in the associated Flag
         */
        public final Object parsedValue;

        public FlagArg(String fullText, int position, FlagDef flag, String flagText, String valueText, Object value) {
            super(fullText, position);
            this.flag = flag;
            this.flagText = flagText;
            this.valueText = valueText;
            this.parsedValue = value;
        }
    }

    /**
     * Unnamed Arguments are those that stand alone on the command line. I.e. they are arguments that are
     * neither a flag nor a parameter value of a flag.
     */
    public static class UnnamedArg extends Argument {
        public UnnamedArg(String text, int index) {
            super(text, index);
        }
    }

    /**
     * Contains information about flags expected on the command-line. This information includes
     * the flag name (and any aliases it may have).
     * whether the flag can be repeated,
     */
    public static class FlagDef {
        public static final String[] NO_ALIASES = {};
        public final String name;
        public final String[] aliases;
        public final ParamAffinity paramAffinity;
        public final boolean allowRepeats;
        public FlagDef(String name, ParamAffinity canHaveParams, boolean allowRepeats) {
            this(name, canHaveParams, allowRepeats, NO_ALIASES);
        }
        public FlagDef(String name, ParamAffinity canHaveParams, boolean allowRepeats, String... aliases) {
            this.name = name;
            this.paramAffinity = canHaveParams;
            this.allowRepeats = allowRepeats;
            this.aliases = (aliases == null || aliases.length == 0 ? NO_ALIASES : aliases);
        }
    }

    public static class ParseResults {
        public final CommandLineParser parser;
        public final List<Argument> arguments;
        public final List<FlagArg> flags;
        public final List<UnnamedArg> unnamed;

        public ParseResults(CommandLineParser parser, List<Argument> args) {
            this.parser = parser; // explicit reference to enclosing class.
            this.arguments = args;
            this.flags = new ArrayList<>();
            this.unnamed = new ArrayList<>();
            for (Argument a : arguments) {
                if (a instanceof FlagArg)
                    flags.add((FlagArg) a);
                else if (a instanceof UnnamedArg)
                    unnamed.add((UnnamedArg) a);
            }
        }

//        /**
//         * Gets a list of all values specified for a specific flag.
//         *
//         * @param flagNameOrAlias The name of the flag for which to retrieve values.
//         * @return A list of all values specified for a specific flag. The return value will never be null, but may be empty (if the flag was not found) or it may have a single value (if the flag was specified only once.)
//         */
//        public List<String> getFlagValues(String flagNameOrAlias) {
//            FlagDef f = parser.findFlagDef(flagNameOrAlias, true, true);
//            if (f == null)
//            return getFlagValues(f);
//        }

//        /**
//         * Gets a list of all values specified for a specific flag.
//         *
//         * @param flag The Flag definition for which to retrieve values.
//         * @return A list of all values specified for a specific flag. The return value will never be null, but may be empty (if the flag was not found) or it may have a single value (if the flag was specified only once.)
//         */
//        public List<String> getFlagValues(FlagDef flag) {
//            List<String> found = Collections.emptyList();
//            for (Argument a : arguments) {
//                if (flag == a.flag && a.value != null) {
//                    if (found.size() == 0)
//                        found = Collections.singletonList(a.value);
//                    else {
//                        if (found.size() == 1)
//                            found = new ArrayList<>(found);
//                        found.add(a.value);
//                    }
//                }
//            }
//            return found;
//        }

        public FlagArg getFlag(String flagNameOrAlias) {
            FlagDef f = parser.findFlagDef(flagNameOrAlias, true);
            if (f == null)
                return null;
            for (FlagArg a : flags)
                if (a.flag == f)
                    return a;
            return null;
        }

        /**
         * Returns the value of the flag that matches the given name or alias.
         * <p>
         * If {@link CommandLineParser.FlagDef#allowRepeats allowRepeats } is true for the flag,
         * the value returned will be that of the <em>first</em> instance of the flag.
         * </p>
         * <p>
         * If a parameter value is optional for this flag (i.e. if the
         * { CommandLineParser.FlagDef#paramAffinity  paramAffinity} of this flag is
         * { CommandLineParser.ParamAffinity#OPTIONAL OPTIONAL}), you will not be able to use the
         * return value to distinguish between a missing flag and a flag that is present but without a value,
         * because {@code null } is returned in both cases.
         * However, you can call {@link #hasFlag(String)} to determine whether the flag was present or not.
         * </p>
         *
         * @param flagNameOrAlias The name of the flag to search for.
         * @return The value of the flag that matches the given name, or null if no match is found.
         */
        public String getFlagValue(String flagNameOrAlias) {
            FlagArg a = getFlag(flagNameOrAlias);
            if (a == null)
                return null;
            return a.valueText;
        }

//        /**
//         * Returns the value of the first flag that matches a name in the list.
//         *
//         * @param flagNames The list of flags to search for.
//         * @return The value of the first flag that matches a name in the list, or null if no name is found. Use @code hasFlag to distinguish between a null valued flag and a missing flag.
//         */
//        public String getFlagValue(String... flagNames) {
//            for (String name : flagNames) {
//                name = normalizeFlagName(name);
//                Flag f = flags.get(name);
//
//            }
//            return null;
//        }

        /**
         * Returns true if the specified flag was found.
         *
         * @param flagNameOrAlias The name or one of the aliases the flag.
         * @return True if the flag was found on the command-line or false otherwise.
         */
        public boolean hasFlag(String flagNameOrAlias) {
            return getFlag(flagNameOrAlias) != null;
        }

        /**
         * Returns true if any of the listed flags were found.
         *
         * @param flagNameList The list of flags to search for.
         * @return True if any of the listed flags were found or false otherwise.
         */
        public boolean hasAnyFlag(String... flagNameList) {
            for (String name : flagNameList)
                if (hasFlag(name))
                    return true;
            return false;
        }
        public boolean foundAllFlags(String... flagNameList) {
            for (String name : flagNameList)
                if (!hasFlag(name))
                    return false;
            return true;
        }

        /**
         * Get the list of Unnamed Arguments (i.e. non-flag values) in the order they were listed on the command-line.
         *
         * @return A non-null list of ordered values.
         */
        @NotNull
        public List<String> getUnnamedValues() {
            if (unnamedValues == null) {
                unnamedValues = new ArrayList<String>(unnamed.size());
                for (UnnamedArg a : unnamed)
                    unnamedValues.add(a.fullText);
            }
            return unnamedValues;
        } private List<String> unnamedValues;

    }

}
