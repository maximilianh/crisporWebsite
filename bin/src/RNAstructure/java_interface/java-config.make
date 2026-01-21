################################################################################
# In order to compile Java programs (and the required C++ JNI library), `make`
# has to be able to locate the Java JDK. (The actual folder path is required,
# not just access to the binary executables on the PATH)
#  
# The location of the JDK differs drastically between platforms and even  
# betweendifferent versions of Java and implementations of it (e.g. Oracle
# vs OpenJDK)
#
# Fortunately this Makefile can auto-detect the location in most cases.


########################### JDK_HOME ###########################################
# In most modern Java installations, Java is installed into a main directory 
# (referred to as 'JDK_HOME' that contains subdirectories 'bin' and 
# 'include' (among others).  
#
# If this is the case on your system, it is only necessary to set the location 
# of JDK_HOME. Otherwise it might be necessary to set both JDK_BIN and 
# JDK_INCLUDE variables, which are described later.
#
# Uncomment and set JDK_HOME if your Java installation contains both 'bin' and 
# 'include' subfolders. 
# (This is only necessary if auto-detection of Java fails on your system.)
#
#JDK_HOME  =


##################### JDK_BIN and JDK_INCLUDE ##################################
#  If your Java installation does NOT place both 'bin' and 'include' 
#  sub-directories together in the same parent directory, then it is necessary 
#  to define *both* of the following variables individually:
#
#  1) JDK_BIN     - The Java JDK 'bin' folder, which contains the Java 
#                   compiler (javac) and archiver (jar).
#  2) JDK_INCLUDE - The JNI 'inlcude' folder, which contains C++ headers 
#                   necessary for compilation of the JNI interface.
#
#  The defaults for JDK_BIN and JDK_INCLUDE (if they are not set below) are
#  respectively:  $(JDK_HOME)/bin   and   $(JDK_HOME)/include
#
# Uncomment and set both JDK_BIN and JDK_INCLUDE if your Java installation 
# places 'bin' in a different location than 'include'.
#
#JDK_BIN    =# JDK_BIN is the directory containing 'javac'. 
#JDK_INCLUDE=# JDK_INCLUDE is the directory containing 'jni.h'. Default: $(JDK_HOME)/include

####################### DISABLE_JAVA_AUTODETECT ################################
# If DISABLE_JAVA_AUTODETECT is set (to anything other than blank) 
# auto-detection of JAVA_HOME, JDK_BIN, and JDK_INCLUDE  will be disabled. 
#
# It is NOT RECOMMENDED that you do this, even if you set the correct path(s) 
# above, because auto-detection also performs additional verifications 
# (e.g. correct java version, existance of required files, etc.).
# So set this ONLY if you have already verified that compilation proceeds 
# correctly for the path(s) defined above.
#
#DISABLE_JAVA_AUTODETECT=DISABLED