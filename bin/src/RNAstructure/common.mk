##############################################################################
#  This file is meant to be included by other Makefiles. 
#  It includes general-purpose variable defintions and functions.
#  Additional functions useful for inspecting variables and debugging
#  makefiles are also included.
##############################################################################

#### Note on Terminology: 
#   In Makefiles, a user-defined "function" is simply a variable 
#   that is meant to be evaluated using the `call` built-in function.
#   To differentiate these from normal variables, the prefix "fn" is 
#   added to the names of all functions that are meant to be accessed 
#   from external Makefiles.
#
#   Example of function usage:
#     # `fnToUpper` converts its first argument to upper case.
#     # After the following line is parsed, MyVar evaluates to "HELLO WORLD"
#     MyVar = $(call fnToUpper,Hello World)

#### Notes on variable naming:
#  (*)  As mentioned above, public functions are prefixed with "fn".
#  (*)  Variables and functions that are meant to be INTERNAL (i.e.
#         not used outside of THIS file) are prefixed with a tilde (~).
#       Example: ~fnInternalFunction or ~INTERNAL_VAR
#  (*)  Variables that are temporary (i.e. redefined in functions or 
#       loops etc) are prefixed with double-tilde (~~).
#       Example: ~~TEMP_VAR
#  (*)  Variables that cause side-effects (such as file-system changes)
#       when evaluated are prefixed with an exclamation point (!).
#       Example: !MAKE_OUTDIR  # Creates the parent directory of a make target.
#
##############################################################################
#  Contents:
#  Section 1. Variable definitions for special characters.
#  Section 2. Global special variable definitions
#  Section 3. General-Purpose Text-Processing Functions 
#  Section 4. Other General-Purpose Functions 
#  Section 5: Variable Inpection and other Debugging functions
#  Section 6. Operating System and Architecture Detection
#  Section 7. Common, Useful Recipe Macros
#  Section 8. Common, Useful Make Targets (for debugging etc)
#
##############################################################################
# The first target in a Makefile becomes the default target. But
# we don't want any of the special targets below to become
# the default.
# .DEFAULT_GOAL is a built-in variable that stores the default makefile target.
# It is probably blank here, unless the client Makefile already listed a target.
# Either way, we will restore it at the end of this (included) makefile.
~~PRESERVE_DEFAULT_GOAL := $(.DEFAULT_GOAL)

##############################################################################
#  Section 1. Variable definitions for special characters.
##############################################################################
#newline (linefeed) character
define C_NEWLINE


endef
|=# define $| to mean empty ("") for convienience in the following defintions.
EMPTY:=
C_COMMA:= ,
C_SPACE:=$| #
C_TAB:=$|	#
# The following symbol definitions are sometimes useful in places where the
# actual symbol cannot be used because it has special meaning.
C_HASH:=\#
C_BKSLASH:=\$|
C_COLON:=:
C_EXCL:=!
C_SEMI=;
C_BAR=|
t:=$(C_TAB)
n:=$(C_NEWLINE)
s:=$(C_SPACE)
k:=$(C_BKSLASH)
h:=$(C_HASH)
,:=$(C_COMMA)#  define $, to expand to a comma
\:=$(C_BKSLASH)#  define $\ to expand to a backslash
$(C_COLON) := $(C_COLON)#  define $: to expand to a colon

#### Definitions for formatted terminal output ###
# Use these in terminal output for bold, underlined, colored text etc.
# Set the following ONLY if make output is being redirected to a terminal.
ifneq ($(MAKE_TERMOUT),)
  C_ESC:=$(shell echo $$'\033')
  E:=$(C_ESC)# Assign $E to be the escape character for the following definitions.
  
  # Normal/Dark colors
  TFMT_T:=$E[39m#  Default Text Color
  TFMT_K:=$E[30m#  Black
  TFMT_R:=$E[31m#  Red
  TFMT_G:=$E[32m#  Green
  TFMT_Y:=$E[33m#  Yellow
  TFMT_B:=$E[34m#  Blue
  TFMT_M:=$E[35m#  Magenta
  TFMT_C:=$E[36m#  Cyan
  TFMT_GY:=$E[90m# Dark Gray
  TFMT_W:=$E[97m#  White
  
  # Light colors
  TFMT_LR:=$E[91m#  Light Red
  TFMT_LG:=$E[92m#  Light Green
  TFMT_LY:=$E[93m#  Light Yellow
  TFMT_LB:=$E[94m#  Light Blue
  TFMT_LM:=$E[95m#  Light Magenta
  TFMT_LC:=$E[96m#  Light Cyan
  TFMT_LGY:=$E[37m# Light Gray
  
  # Bold/Bright colors
  TFMT_BR:=$E[1;31m# Bold/Bright Red
  TFMT_BG:=$E[1;32m# Bold/Bright Green
  TFMT_BY:=$E[1;33m# Bold/Bright Yellow
  TFMT_BB:=$E[1;34m# Bold/Bright Blue
  TFMT_BM:=$E[1;35m# Bold/Bright Magenta
  TFMT_BC:=$E[1;36m# Bold/Bright Cyan
  
  # Other formatting
  TFMT_RST:=$E[0m#  Reset all attributes (formatting and color)
  TFMT_BLD:=$E[1m#  Bold
  TFMT_DIM:=$E[2m#   Dim
  TFMT_UND:=$E[4m#   Underlined
  TFMT_HID:=$E[8m#   Hidden (for passwords etc)
  
  TFMT_RBLD:=$E[21m#  Reset Bold
  TFMT_RDIM:=$E[22m#  Reset Dim
  TFMT_RUND:=$E[24m#  Reset Underline
  TFMT_RHID:=$E[28m#  Reset Hidden
  
  TFMT0 := $(TFMT_RST)# Alias for Reset All
  TFMT! := $(TFMT_BLD)# Alias for Bold
  TFMT_ := $(TFMT_UND)# Alias for Underline
  
  #### Single-character formatting Aliases for easier use in messages: ####
  $(C_EXCL):=$(TFMT_BLD)#  Bold Text with "$!"
  $(C_SEMI):=$(TFMT_RST)#  Reset format with "$;"
  _:=$(TFMT_UND)#          Unerline with $_
  # R:=$(TFMT_BR)# $R = Bold Red
  # G:=$(TFMT_GG)# $G = Bold Green
  # Y:=$(TFMT_BY)# $Y = Bold Yellow
  # B:=$(TFMT_BB)# $B = Bold Blue
  #$(info $!BoldText$; $(TFMT_BR)Bold Red Text$; $(TFMT_R) Red Text$;)
endif


##############################################################################
#  Section 2. Global special variable definitions
##############################################################################

# Standard Exclusions when using `shell find`
# -- Exclude directories that start with '.' or end with .bak, .old, _bak or _old.
# -- Exclude files that end with .bak or .old and files of the form  '<FILENAME>_old<NUMBER>.<EXT>' e.g. 'structure_old2.cpp'
SHELL_FIND_EXCLUDED_DIRS:=\( -type d \( -name '.*' -o -iname '*[._]bak' -o -iname '*[._]old' \) -prune -false \) -o  # Note that '-false' here is used to disable automatic printing after -prune.
SHELL_FIND_EXCLUDED_FILES:=\( -not \( -type f \( -iname '*.bak' -o -iname '*.old' -o -iname '*_old[0-9].*'  \) \) \)
SHELL_FIND_EXCLUDE=$(SHELL_FIND_EXCLUDED_DIRS) $(SHELL_FIND_EXCLUDED_FILES)


##############################################################################
#  Section 3. General-Purpose Text-Processing Functions 
##############################################################################

####  fnToLower -- Convert text to lower case. ####
#  Usage: $(call fnToLower,SOME_TEXT)
#  Example: $(call fnToLower,Hello World) --> hello world
fnToLower = $(subst A,a,$(subst B,b,$(subst C,c,$(subst D,d,$(subst E,e,$(subst F,f,$(subst G,g,$(subst H,h,$(subst I,i,$(subst J,j,$(subst K,k,$(subst L,l,$(subst M,m,$(subst N,n,$(subst O,o,$(subst P,p,$(subst Q,q,$(subst R,r,$(subst S,s,$(subst T,t,$(subst U,u,$(subst V,v,$(subst W,w,$(subst X,x,$(subst Y,y,$(subst Z,z,$1))))))))))))))))))))))))))
~TEST_fnToLower = $(call fnToLower,ABCDEFGHIJKLMNOPQRSTUVWXYZ)
~RESULT_fnToLower = abcdefghijklmnopqrstuvwxyz

####  fnToUpper -- Convert text to upper case. ####
#  Usage: $(call fnToUpper,SOME_TEXT)
#  Example: $(call fnToUpper,Hello World) --> HELLO WORLD
fnToUpper = $(subst a,A,$(subst b,B,$(subst c,C,$(subst d,D,$(subst e,E,$(subst f,F,$(subst g,G,$(subst h,H,$(subst i,I,$(subst j,J,$(subst k,K,$(subst l,L,$(subst m,M,$(subst n,N,$(subst o,O,$(subst p,P,$(subst q,Q,$(subst r,R,$(subst s,S,$(subst t,T,$(subst u,U,$(subst v,V,$(subst w,W,$(subst x,X,$(subst y,Y,$(subst z,Z,$1))))))))))))))))))))))))))
~TEST_fnToUpper = $(call fnToUpper,abcdefghijklmnopqrstuvwxyz)
~RESULT_fnToUpper=ABCDEFGHIJKLMNOPQRSTUVWXYZ

# Converts a list of space-separated values into a newline-separated list (i.e. one value per line)
# Usage: $(call fnSpacesToLines,WORD_LIST)
fnSpacesToLines=$(subst $s,$n,$(strip $1))
~TEST_fnSpacesToLines = $(call fnSpacesToLines,apple  banana cherry   dog fred)
define ~RESULT_fnSpacesToLines
apple
banana
cherry
dog
fred
endef

fnEscapeChars=$(subst $(C_NEWLINE),[LF],$(subst $(C_ESC),[ESC],$(subst $(C_TAB),[TAB],$1)))
#$(subst $(subst $(subst $(subst $(subst $(subst )))))))))

#### fnIsTruthy -- Return T for a truthy value or empty otherwise. ####
# Usage: $(call fnIsTruthy,TEXT,[T_VAL],[F_VAL])
# Description: If TEXT is "truthy" (yes true on 1 Y or T) then T_VAL is returned.
#    Otherwise F_VAL is returned. T_VAL defaults to 'T' while F_VAL defaults to 
#    empty ('').
fnIsTruthy=$(if $(filter yes true on 1 y t,$(call fnToLower,$1)),$(or $2,T),$3)

#### fnIsFalsey -- Return F for a falsey value or empty otherwise. ####
# Usage: $(call fnIsFalsey,TEXT,[F_VAL],[T_VAL])
# Description: If TEXT is "falsey" (no false off 0 N or F) then F_VAL is returned.
#    Otherwise T_VAL is returned. F_VAL defaults to 'F' while T_VAL defaults to 
#    empty ('').
fnIsFalsey=$(if $(filter no false off 0 n f,$(call fnToLower,$1)),$(or $2,F),$3)

#### fnToBool -- Convert a text value into a boolean (T/F) ####
# $(call fnTriBool TEXT,[T_VAL],[F_VAL],[OTHER_VAL],[EMPTY_VAL])
#  TEXT  -- the value to be tested.
#  T_VAL -- the value to return if TEXT is "truthy" (yes, true, 1, T, on).   (default is 'T')
#  F_VAL -- the value to return if TEXT is "falsey" (no, false, 0, F, off).  (default is 'F')
#  OTHER_VAL -- the value to return if TEXT is neither truthy nor falsey nor blank. (default is 'T')
#  EMPTY_VAL -- the value to return if TEXT is empty/blank ('') (default is 'F')
# Examples:
#   $(call fnToBool,on)    --> T
#   $(call fnToBool,off)   --> F
#   $(call fnToBool,hello) --> T
#   $(call fnToBool,)      --> F
#   $(call fnToBool,on,Great,Terrible) --> Great
#   $(call fnToBool,F,Great,Terrible)  --> Terrible
#   $(call fnToBool,,Great,Terrible)  --> F  #(because EMPTY_VAL was not specified and defaults to 'F')
#   $(call fnToBool,,Great,Terrible,Wow,Nothing)  --> Nothing
fnToBool=$(call fnToBoolEx,$1,$(or $2,T),$(or $3,F),$(or $4,T),$(or $5,F))

#### fnToBoolEx ####
# Description: Converts a text value into one of four values, depending on whether 
#  it is truthy, falsey, other, or blank.
# $(call fnToBoolEx TEXT,T_VAL,F_VAL,OTHER_VAL,EMPTY_VAL)
#  TEXT  -- the value to be tested.
#  T_VAL -- the value to return if TEXT is "truthy" (yes, true, 1, T, on).
#  F_VAL -- the value to return if TEXT is "falsey" (no, false, 0, F, off).
#  OTHER_VAL -- the value to return if TEXT is neither truthy nor falsey nor blank.
#  EMPTY_VAL -- the value to return if TEXT is empty/blank ('')
# Unlike fnToBool, the default for any of the return values (including T_VAL) is blank! So they must be specified!
# Examples:
#   $(call fnToBool,on)        --> '' (empty -- because no T_VAL was specified)
#   $(call fnToBool,off,T,F)   --> F
#   $(call fnToBool,hello,T,F) --> ''  (empty -- because no OTHER_VAL was specified)
#   $(call fnToBool,,T,F)      --> ''  (empty -- because no EMPTY_VAL was specified)
#   $(call fnToBool,hello,T,F,O,E) --> O
fnToBoolEx=$(if $(call fnIsTruthy,$1),$2,$(if $(call fnIsFalsey,$1),$3,$(if $1,$4,$5)))

#### fnNot ####
# Description: Converts a truthy value to F or a falsey value to T. 
#              Returns blank for non-boolean values.
# Usage: $(call fnNot,TEST)
fnNot=$(call fnToBoolEx,$1,F,T)

#### fnIf ####
# Description: Converts a truthy value to T_VAL or a falsey value to F_VAL. 
#   Returns blank for non-boolean values.
# Note that this differs from the built-in $(if) function in a number of cases:
#    fnIf                          |  Built-in $(if)        | Comparison
#    ------------------------------|------------------------|------------
# 1. $(call fnIf,yes,A,B)   ==> A  |  $(if yes,A,B)   ==> A | same -- Both functions return T_VAL for truthy (non-empty) values.
# 2. $(call fnIf,pizza,A,B) ==> '' |  $(if pizza,A,B) ==> A | diff -- fnIf returns empty for non-boolean values.
# 3. $(call fnIf,no,A,B)    ==> B  |  $(if no,A,B)    ==> A | diff -- fnIf returns F_VAL for falsey (non-empty) values.
# 4. $(call fnIf,,A,B)      ==> '' |  $(if pizza,A,B) ==> B | diff -- fnIf returns empty for empty values.
#
# Usage: $(call fnIf,TEST,T_VAL,F_VAL)
fnIf= $(call fnToBoolEx,$1,$2,$3)

##############################################################################
#  Section 4. Other General-Purpose Functions 
##############################################################################

####  fnAddAlias  ####
#  Description: fnAddAlias -- adds an alias for a target, but only 
#               if the alias is NOT identical to the target.
#  Usage: $(call fnAddAlias,ALIAS,TARGET)
#  Example:  $(call fnAddAlias,oligowalk,OligoWalk) --> equivalent to rule "oligowalk: OligoWalk"
#  Example:  $(call fnAddAlias,bifold,bifold) --> No effect (since alias and target are identical)
fnAddDep = $(if $(filter-out $1,$2),$(eval $1: $2))

####  fnMkdir ####
# Description: Produces shell code that checks whether a directory exists
#   and creates it if it doesn't exist. This uses mkdir -p, so an entire 
#   directory tree can be created.
# Usage: $(call fnMkdir,PATH_TO_DIR)
fnMkdir = @[ -z "$1" -o -d "$1" ] || mkdir -p "$1"

####  fnToColumns  ####
# Description: This function formats whitespace-separated words into columns 
# (of equal character-width).
# Usage: $(call fnToColumns,WORDS,[LEFT_INDENT],[COLUMN_WIDTH],[MAX_LINE_WIDTH])
# Parameters: 
#   WORDS -- A list of words that should be split into columns. The words are 
#     white-space separated.
#   LEFT_INDENT -- Text that should be prepended to each line in the table to 
#     create a left margin. Typically one of: "" or "\t" or "    "  etc.
#   COLUMN_WIDTH -- The length (in characters) of each column. Words that exceed
#     this length will occupy more than one column, pushing subsequent words into
#     the next column. Typically 12 to 18 etc.
#   MAX_LINE_WIDTH -- The maximum length (in characters) of any line in the table.
#     Typically 80.
# Example: $(call fnToColumns,$(PROGRAM_LIST),\t,18,80)
fnToColumns=$(shell $(call ~TO_COLUMNS_BASH_CODE,$1,\n$(or $2,$s$s$s$s),$(or $3,18),$(or $4,80)))
define ~TO_COLUMNS_BASH_CODE
  P=$3 T=0; for W in $(patsubst %,'%',$1); do
    ((L=$${#W}/P*P+P)); (((T+=L)>$4)) && { printf '\n%s' '$2'; T=$$L; };
    printf "%-$${L}s" "$$W";
  done
endef

####  fnCacheVar  ####
# Description: Use to create a macro that gets redefined after it is first 
#    expanded, so that its result is cached.
# Usage:  $(call fnCacheVar,VAR_NAME,ESCAPED_EXPENSIVE_CODE)
# Parameters: 
#    VAR_NAME -- The name of the variable to define.
#    ESCAPED_EXPENSIVE_CODE -- Makefile code (a function etc) that 
#      will be evaluated at most once to determine the final true value
#      of the variable. Dollar signs should in most cases be escaped ($$).
# Example:  $(call fnCacheVar,OSNAME,$$(shell uname))
#   This way uname will not be invoked unless Make expands OSNAME. But
#   when OSNAME is expanded for the first time, it is redefined to be
#   the result of the shell code, so that code is never evaluated more 
#   than once.
#
# Discussion:
#     Assume you have a macro that invokes a shell program:
#       (1) MyMacro = $(shell /bin/SomeProgram)
#       (2) MyMacro := $(shell /bin/SomeProgram)
#   Assignment #1 is inefficient when MyMacro is expanded multiple times,
#   because the shell code is re-executed each time MyMacro is evaluated. 
#   Assignment #2 is more efficient, but the downside of it is that the 
#   assignment happens immediately and it happens every time the makefile 
#   is parsed, even if MyMacro is never evaluated.
#   The best approach is a combination of the two: MyMacro is initially 
#   defined using "=" (aka recursively-evaluated) but the first time it 
#   is evaluated, it redefines itself using := so that future evaluations
#   are expanded to a simple scalar value instead of re-invoking the shell 
#   code.
#   This can be achieved as follows:
#       (3) MyMacro = $(eval MyMacro := $$(shell /bin/SomeProgram))$(MyMacro)
#   This way, MyMacro is NOT immediately evaluated. But the first time it is 
#   evaluated, it redefines itself. The final $(MyMacro) at the end causes the
#   result to be returned the time the macro is evaluated. (Because eval does 
#   not return anything itself.)
fnCacheVar=$(eval $1 = $$(eval $1 := $2)$$($1))

#### fnGetOption ####
# Description: Expands the variable name passed in as the first argument, 
#   but ONLY if it was specified as an option on the command-line to `make`
# Usage:  $(call fnGetOption,OPTION_NAME,DEFAULT_VALUE)
# Parameters:
#   OPTION_NAME   -- the literal name of the option (e.g. 'debug')
#   DEFAULT_VALUE -- (optional) the value to return if the option was 
#                    not specified. Typically blank.
# Example:  `make all debug=on`        # Command-line to make
#           $(call fnGetOption,debug)  # expands to "on"
#           $(call fnGetOption,CXX)    # expands to "" (even if CXX is defined)
#           $(call fnGetOption,display,verbose)    # expands to "verbose"
fnGetOption=$(if $(findstring command line,$(origin $1)),$($1),$2)

#### fnIsOptionON ####
# Description: Expands to a non-empty value ONLY if the specified 
#   option was passed on the command-line to `make` AND the value of the option
#   is "1", "yes", "true", or "on" (case insensitive)
# Usage:  $(call fnIsOptionON,<OPTION_NAME>)
# Example:  
#   `make all debug=on optimize=off`       # Command-line to make
#   $(call fnIsOptionON,debug)    # expands to "on"
#   $(call fnIsOptionON,optimize) # expands to ""
#   $(call fnIsOptionON,CXX)      # expands to "" (even if CXX is defined)
fnIsOptionON= $(call fnIsTruthy,$(call fnGetOption,$1),on)

#### fnIsOptionOFF ####
# Description: Expands to a non-empty value ONLY if the specified 
#   option was passed on the command-line to `make` AND the value of the option
#   is "0", "no", "false", or "off" (case insensitive)
# Usage:  $(call fnIsOptionOFF,<OPTION_NAME>)
# Example:  
#   `make all debug=on optimize=off`       # Command-line to make
#   $(call fnIsOptionOFF,debug)    # expands to ""
#   $(call fnIsOptionOFF,optimize) # expands to "off"
#   $(call fnIsOptionOFF,CXX)      # expands to "" (even if CXX is defined)
fnIsOptionOFF=$(call fnIsFalsey,$(call fnGetOption,$1),off)

#### fnHandleOption ####
# Description: Useful to handle settings that can have default values of 
#   ON or OFF (or True/False etc) that can also be forced ON or OFF by
#   a command-line option (flag).
# Usage: 
#   $(call fnHandleOption,[DEFAULT_VALUE],[FLAG_ENABLE],[FLAG_DISABLE],[ON_VALUE=on],[OFF_VALUE])
# Parameters:
#   DEFAULT_VALUE -- (optional) The default value of the setting. Any NON-truthy value is 
#                    interpreted as OFF.
#   FLAG_ENABLE   -- (optional) The name of the command-line flag that controls
#      this setting in a direct/positive way -- i.e. the setting will be ON
#      if the flag has a truthy value and OFF if the flag has a falsey value.
#      (The flag will be ignored if it does not have a boolean-type value.)
#   FLAG_DIABLE   -- (optional) The name of the command-line flag that controls 
#      this setting in an inverted/opposite way -- i.e. the setting will be OFF
#      if the flag has a truthy value and ON if the flag has a falsey value.
#      E.g. setting the flag 'dynamic=off' turns the COMPILE_STATIC setting ON.
#      (The flag will be ignored if it does not have a boolean-type value.)
#    ON_VALUE     -- (optional) The value that should be returned if the setting 
#                    is enabled/ON. The default is 'on'.
#    OFF_VALUE    -- (optional) The value that should be returned if the setting 
#                    is disabled/OFF. The default is empty/blank ('').
# Examples:
#   $(call fnHandleOption,$(DEFAULT_DEBUG),debug)
#   $(call fnHandleOption,$(DEFAULT_STATIC),static,dynamic)
#   VERBOSITY:=$(call fnHandleOption,$(DEFAULT_VERBOSE),verbose,,-V,-Q)
fnHandleOption=$(call fnIf,$(or \
  $(and $2,$(call fnToBoolEx,$(call fnGetOption,$2),T,F)),\
  $(and $3,$(call fnToBoolEx,$(call fnGetOption,$3),F,T)),\
  $(call fnToBool,$1)\
),$(or $4,on),$5)

# Logic schematic for fnHandleOption 
# Note: T* is any truthy value, F* is any (non-empty) falsey value. (See fnIsTruthy, fnIsFalsey)
#  if FLAG_ENABLE  is T*   return ON_VALUE
#  if FLAG_ENABLE  is F*   return OFF_VALUE
#  if FLAG_DISABLE is T*   return OFF_VALUE
#  if FLAG_DISABLE is F*   return ON_VALUE
#  if DEFAULT      is T*   return ON_VALUE
#  else                    return OFF_VALUE

# fnEvalForeach -- Evaluate an expression once for each word in the list, replacing the text "{?}" with each word.
# Usage: $(call fnEvalForeach,<WORD_LIST>,<EXPRESSION>)
# This will call $(eval <EXPRESSION>) once for each word in <WORD_LIST> with the text "{?}" replaced with that word.
# IMPORTANT: The <EXPRESSION> argument is expanded twice; first by the eval function, then the results of that expansion 
#   are expanded AGAIN when they are parsed as makefile syntax. This means you may need to provide extra levels 
#   of escaping for “$” characters when using eval. The value function can sometimes be useful in these situations, 
#   to circumvent unwanted expansions.
fnEvalForeach=$(foreach ~~tmp_var,$1,$(eval $(subst {?},$(~~tmp_var),$2)))

# fnFind -- Run the shell `find` command, excluding specified files and dirs.
# Usage: $(call fnFind,<DIRS>,<FIND_EXPR>,<EXCLUDE_DIRS_EXPR>,<EXCLUDE_FILES_EXPR>)
#   DIRS - A directory or space-separated list of directories to search. 
#          These should be quoted if they contain spaces or other special chars.
#   FIND_EXPR - a `find` expression for the files or directories to match.
#          for example: "-type f -iname '*.cpp'"     or     "-path '*/proxy'"
#   EXCLUDE_DIRS_EXPR - (optional) A `find` expression for directories to prune.
#          for example: "-name bak -o -name '*_old'"
#          (note that "-type d" and "-prune" flags are implied and should NOT be included in the expression)
#   EXCLUDE_FILES_EXPR - (optional) A `find` expression for files or directories to exclude.
#          for example: "-iname readme"   or   "-name '*.o' -o -path 'exe/*.so"
fnFind=$(shell find $1 $(SHELL_FIND_EXCLUDE) \( -type d \( $(or $3,-false) \) -prune -false \) -o \! \( -type f \( $(or $4,-false) \) \) $2 )
~TEST_fnFind: # goal to test the fnFind function.
	@echo "Find cpp files in src. Exclude phmm and dynalign*"
	@printf '%s\n'  $(call fnFind,src,-iname '*.cpp',-path src/phmm,-iname 'dynalign*')

####################################################
####  Recipe Functions
####################################################
# fnCpDir -- copy a directory recursively in a reliable cross-platform way.
#   Assume we want to recursively copy dirA into dir2 and name it "B": dirA ==> dir2/B
#  	    `cp -a dirA dir2/B` will not work reliably, because its behavior depends on the 
#   existance of dir2/B: If dir2/B does NOT exist, dirA will be copied to dir2/B (as expected)
#   But if dir2/B already exists, dirA will be will copied *INTO* it: dir2/B/dirA (not what we want)
#   Appending "/." to the source directory will do what we want:
#       `cp -a dirA/. dir2/B` will reliably copy dirA to dir2/B whether dir2/B exists or not.
#   Flag -a means copy recursively, preserving file types and permissions/attributes when possible.
#   Minor Notes:
#   - Do not append a dot to the destination directory. `cp -a dirA/. dir2/B/.` will fail if 
#     dir2/B doesn't exist.
#   - BSD `cp` (MacOSX) differs from GNU (Linux|Cygwin) in behavior and flags.
#     E.g. `cp -a dirA/ dir2/B` (no dot) will work on BSD cp but not GNU.
#     On MacOSX, the `cp` command does not support -u which means 'update'
#      -- i.e. do not re-copy existing, unchanged files. So rsync is preferred on Mac if it exists.
#
# Usage: $(call fnCpDir,SRC_DIR,DEST_DIR) -- produces shell code to copy the directory SRC_DIR to DEST_DIR
fnCpDir=$(!CP) "$1/." "$2"

# fnIsExe -- Determine if a program (or shell builtin) can be run/executed.
# Outputs 1 if the given argment is a command or executable.
# Output is blank if the argment is a normal file or does not exist.
# Usage: $(call fnIsCmd,COMMAND)
# Example: COPY_CMD=$(if $(call fnIsCmd,rsync),rsync -au,cp -au) # use `rsync` to copy if it exists, otherwise use `cp`
fnIsCmd=$(shell type -p "$1" >/dev/null && echo 1)


##############################################################################
#  Section 5: Variable Inpection and other Debugging functions
#############################################################################

#### fnShowVar - Show information about a single variable ####
# Usage: $(call fnShowVar,VAR_NAME,[FIELD_SEP],[LABEL_SEP])
# Parameters: VAR_NAME is the name of the variable to show information for.
#             FIELD_SEP is text that separates the fields (Source,Define,Value) 
#             LABEL_SEP is text that comes after the labels ("Source","Define","Value")
# Output: VAR_NAME>>>>Source:--VAR_SOURCE>>>>Define:--VAR_DEF>>>>Value :--VAR_VALUE  ### where FIELD_SEP=">>>>" and LABEL_SEP="--"
fnShowVar=$(call ~SHOW_VAR_CODE,$1,$(or $2,$(C_NEWLINE)$(C_SPACE)$(C_SPACE)),$(or $3,$2,$(C_SPACE))) # Specify defaults for field separators ($2 and $3)
 ~SHOW_VAR_CODE=$1$2Source=$3$(origin $1)$2Define=$3\
	$(call fnEscapeChars,$(value $1))$2Value =$3\
	$(eval ~~DISABLE_INFO := 1)$(call fnEscapeChars,$($1))$(eval ~~DISABLE_INFO := )

fnInfo=$(if $(~~DISABLE_INFO),$$(info $1),$(info $1))# Do not show info when displaying the value of all variables.

# fnShowVarList - Show information about each variable in a list of whitespace-separated names.
# Usage: $(call fnShowVarList,VARIABLE_NAMES)
fnShowVarList=$(foreach var,$1,$(C_NEWLINE)  $(call fnShowVar,$(var),$(C_NEWLINE)    ,$(C_SPACE)))

# fnFilterVarsByOrigin - Filter a list of variables, including only those with the specified origins.
# Note that spaces in the origin are replaced by underscore: e.g.: command_line and environment_override
fnFilterVarsByOrigin=$(foreach var,$2,$(if $(filter $1,$(subst $s,_,$(origin $(var)))),$(var)))
# fnFilterVarsByOrigin - Filter a list of variables, excluding those with the specified origins.
fnFilterOutVarsByOrigin=$(foreach var,$2,$(if $(filter-out $1,$(subst $s,_,$(origin $(var)))),$(var)))

##############################################################################
#  Section 6. Operating System and Architecture Detection
##############################################################################

# Call `uname` to get information about the OS/kernel
SHELL_UNAME:=     $(shell uname -srm)
OS_NAME:=   $(word 1,$(SHELL_UNAME))
OS_VERSION:=$(word 2,$(SHELL_UNAME))
OS_ARCH:=   $(word 3,$(SHELL_UNAME))
#SHELL_UNAME=$(call fnCacheVar,SHELL_UNAME,$(shell uname -srm))
#SHELL_OS_NAME=$(word 1,$(SHELL_UNAME))
#SHELL_OS_VERSION=$(word 2,$(SHELL_UNAME))
#OS_ARCH=$(word 3,$(SHELL_UNAME))

# Do some replacements to standardize the OS name.
OS_NAME := $(OS_NAME:CYGWIN%=Windows)
OS_NAME := $(OS_NAME:MSYS%=Windows)
OS_NAME := $(OS_NAME:Darwin=Mac)
OS_NAME := $(OS_NAME:GNU%=Linux)

# Do some replacements to standardize the Architecture (x86 or x86_64)
OS_ARCH := $(OS_ARCH:i%86=x86)

# Get architecture bits (32 or 64)
OS_BITS := $(if $(findstring 64,$(OS_ARCH)),64,32)
ARCH_BITS := $(OS_BITS)# Alias for backwards compatibility.

#Debug:
# $(info OS -- $(OS_NAME) $(OS_ARCH)  $(OS_BITS))

# OS_PATH_SEP is used to separate directory paths.
# On Mac and Linux, this is ":" but on windows it is ";"
# Specifically, when invoking java with a -cp (classpath), the paths must be 
# separated by ";" on windows (even if run from Cygwin which itself uses ":")
# This is modified for Windows in the OS-specific section further below.
ifeq ($(OS_NAME),Windows)
  OS_PATH_SEP:=;
else
  OS_PATH_SEP:=:
endif

# OS_LIB_EXT is the extension of shared object files (libraries)
ifeq ($(OS_NAME),Windows)
  OS_LIB_EXT=dll
else ifeq ($(OS_NAME),Linux)
  OS_LIB_EXT=so
else ifeq ($(OS_NAME),Mac)
  OS_LIB_EXT=dylib
endif

MAKE_CD=$(MAKE) --no-print-directory -C 

# define GET_OS_NAME_BASH_CODE
#   # First output OS Name
#   case $$OSTYPE in
#     cygwin|msys)    echo Windows ;;
#     linux*)         echo Linux   ;;
#     darwin*)        echo Mac     ;;
#     openbsd*|FreeBSD)  echo BSD  ;;
#     *) #OSTYPE is undefined. Use uname instead.
#        case $$(uname -s) in
#          CYGWIN*|MSYS*) echo Windows ;;
#          Linux)   echo Linux ;;
#          Darwin)  echo Mac   ;;
#          *BSD)    echo BSD   ;;
#          *) echo UNKNOWN     ;;
#        esac ;; # end of empty OSTYPE
#   esac
#   # Next output architecture
#   case $$OSTYPE in
#       else 
#         # PREVIOUS method was to use Use the `uname` command to determine the OS.
#         # This might still be required for some operating systems
#         #   OPSYSTEM:=$(shell uname -s 2>/dev/null || echo UNKNOWN)#   Replace 'UNKNOWN' with default OS if desired. 
#         # Perform some replacements to normalize the output of uname on various systems.
#         #OPSYSTEM := $(OSTYPE:cygwin=Windows)
#         #OPSYSTEM := $(OSTYPE:CYGWIN%=Windows)
#         #OPSYSTEM := $(OPSYSTEM:MSYS%=Windows)
#         #OPSYSTEM := $(OPSYSTEM:Darwin=Mac)
#         #OPSYSTEM := $(OPSYSTEM:GNU%=Linux)
#       endif
#       $(if $(DEBUG), $(info Make: Operating System: $(OPSYSTEM)))
#       export OPSYSTEM #make it available for recursive calls to make, so auto-detection is performed only once.
#     endif


##############################################################################
#  Section 7. Common, Useful Recipe Macros
##############################################################################

# !MAKE_OUTDIR: A macro to use in recipes to automatically
#    create the parent directory of the target file.
#    Usage: $(!MAKE_OUTDIR)  # in a recipe
 !MAKE_OUTDIR = $(call fnMkdir,$(@D))

# Recipe command to recursively copy files or directories
# rsync is used on Mac because the `cp` command does not
# have the same features as GNU cp.
ifneq ($(call fnIsCmd,rsync),)
  !CP:=rsync -au
else ifeq ($(OS_NAME),Mac)
  # Use `cp` command to copy in the absence of rsync
  # OSX/BSD `cp` command does not have -u flag.
  !CP:=cp -a
else
  # Include -u option to only update changed files.
  !CP:=cp -au
endif

##############################################################################
#  Section 8. Common, Useful Make Targets (for debugging etc)
##############################################################################

# Filter out these variables from the showallvars target.
~FILTER_ALL_VARS= $(.HIDDEN_VARS) fn% ~fn% ~% ~~% !% TFMT_%

# Filter out these origins from the showallvars target.
~FILTER_ALL_VARS_ORIGINS = automatic default environment
# Possible Variable origins: (from the $(origin) function)
# undefined     -- never defined.
# default       -- has a default definition, as is usual with CC and so on. See Variables Used by Implicit Rules. Note that if you have redefined a default variable, the origin function will return the origin of the later definition.
# environment   -- inherited from the environment provided to make.
# environment override  -- inherited from the environment provided to make, and is overriding a setting for variable in the makefile as a result of the ‘-e’ option (see Summary of Options).
# file          -- defined in a makefile.
# command line  -- defined on the command line.
# override      -- defined with an override directive in a makefile 
# automatic     -- an automatic variable defined for the execution of the recipe for each rule 
# Note: When fnShowVarList filters results, it places an underscore in command_line and environment_override.

~~SHOW_VARS_END=$n

.PHONY: show-vars
show-vars: ; @:
	$(info Make Variables: $n(Excluding hidden, internal, temporary, $(subst $s,$(comma) ,$(strip $(~FILTER_ALL_VARS_ORIGINS))), and function variables.) \
		$(call fnShowVarList,$(sort $(filter-out $(~FILTER_ALL_VARS), \
			$(call fnFilterOutVarsByOrigin,$(~FILTER_ALL_VARS_ORIGINS),$(.VARIABLES))\
		)))\
	$(~~SHOW_VARS_END))

.PHONY: show-env-vars
show-env-vars: ; @:
	$(info Environment Variables: \
		$(call fnShowVarList,$(sort $(filter-out $(~FILTER_ALL_VARS), \
			$(call fnFilterVarsByOrigin,environment environment_override,$(.VARIABLES))\
		)))\
	$(~~SHOW_VARS_END))

.PHONY: show-make-vars
show-make-vars:  ; @:
	$(info Default (Built-In) Make Variables: \
		$(call fnShowVarList,$(sort $(filter-out $(~FILTER_ALL_VARS), \
			$(call fnFilterVarsByOrigin,default,$(.VARIABLES))\
		)))\
	$(~~SHOW_VARS_END))

.PHONY: show-all-vars
show-all-vars: ~~SHOW_VARS_END+=$n----------------------------$n
show-all-vars: show-vars show-env-vars show-make-vars ; @:

#### fnRunTest - Run a Unit Test ###
# TEST_<NAME>
# <FUNCTION>===<RESULT>
# 
# ~fnRunUnitTest=$(eval $(call RUN_TEST_CODE,$(patsubst TEST_%,RESULT_%,$1),$1))
 ~fnRunUnitTest=$(eval $(call ~RUN_UNIT_TEST_CODE,$(subst ~TEST_,,$1),$1,$(subst ~TEST_,~RESULT_,$1))) #$(eval )

# $(call RUN_TEST_CODE,ACTUAL_RESULT,EXPECTED_RESULT,TEST_NAME)
define ~RUN_UNIT_TEST_CODE
  ifeq ($$($2),$$($3))
    $$(info PASSED Test $1)
  else
  	$$(warning FAILED Test $1:$$n$$tExpected: "$$($3)"$$n$$tActual: "$$($2)")
  endif
endef


# Run all Unit-Tests (variables of the form ~TEST_*)
run-unit-tests: ; :
	$(foreach var,$(filter ~TEST_fn%,$(.VARIABLES)),$(call ~fnRunUnitTest,$(var)))


#Restore the default goal if it was set to one of the targets in this file. (It should be set in the Makefile itself.)
.DEFAULT_GOAL := $(~~PRESERVE_DEFAULT_GOAL)


# $(info banana = $(call fnGetOption,banana))
# $(info banana ON = $(call fnIsOptionON,banana))
# $(info banana OFF = $(call fnIsOptionOFF,banana))
# $(info banana Truthy = $(call fnIsTruthy,$(banana)))
# $(info banana Falsey = $(call fnIsFalsey,$(banana)))
# $(info banana Bool = $(call fnToBool,$(banana)))
# $(info banana Bool = $(call fnToBoolEx,$(banana),T,F,O,E))
# $(info banana Not = $(call fnNot,$(banana)))
# $(info banana IF = $(call fnIf,$(banana),AAA,BBB))

# $(info orange = $(call fnGetOption,orange))
# $(info orange ON = $(call fnIsOptionON,orange))
# $(info orange OFF = $(call fnIsOptionOFF,orange))
# $(info orange Truthy = $(call fnIsTruthy,$(orange)))
# $(info orange Falsey = $(call fnIsFalsey,$(orange)))
# $(info orange Bool = $(call fnToBool,$(orange)))
# $(info orange Bool = $(call fnToBoolEx,$(orange),T,F,O,E))
# $(info orange Not = $(call fnNot,$(orange)))
# $(info orange IF = $(call fnIf,$(orange),AAA,BBB))

# $(info pizza = $(call fnGetOption,pizza))
# $(info pizza ON = $(call fnIsOptionON,pizza))
# $(info pizza OFF = $(call fnIsOptionOFF,pizza))
# $(info pizza Truthy = $(call fnIsTruthy,$(pizza)))
# $(info pizza Falsey = $(call fnIsFalsey,$(pizza)))
# $(info pizza Bool = $(call fnToBool,$(pizza)))
# $(info pizza Bool = $(call fnToBoolEx,$(pizza),T,F,O,E))
# $(info pizza Not = $(call fnNot,$(pizza)))
# $(info pizza IF = $(call fnIf,$(pizza),AAA,BBB))

#$(info $(call fnHandleOption,$(pizza),banana,orange,1,0))