RNAstructure/java_interface/SWIG ReadMe

There are two types of files in this folder, .cxx and .i. Each file type has a
different use when building the functional RNAstructure program.

1.   I Files
     Files with a .i extension are SWIG interface files. These files are needed
     by SWIG to generate the appropriate proxy files to link the Java interface
     with the C++ back end.

     In the Makefile located in RNAstructure/java_interface, these
     files are used in the swigFiles target when running SWIG. If a user wants
     to build all SWIG-generated files from scratch, these .i files enable that
     capability. Unless SWIG is installed on the current machine, however, it's
     recommended that one not change or use these files.

     Use the following command in the aforementioned Makefile to build
     the program including a SWIG file rebuild:
     make all

2.   CXX Files
     Files with a .cxx extension are SWIG-generated C++ files. These files are
     proxy files which can be linked into an object library using the buildC
     target in the Makefile found in RNAstructure/java_interface.
     These files allow RNAstructure to be built without the use of SWIG, since
     all the other required files (Java files are located in the
     RNAstructure/java_interface/source directory) are present and
     no additional build steps are needed.

     These files should not be changed in any way. If an issue arises in them,
     a bug should be reported to the RNAstructure developers.

     Use the following make command in the aforementioned Makefile to build the
     program without a SWIG file rebuild:
     make GUI
