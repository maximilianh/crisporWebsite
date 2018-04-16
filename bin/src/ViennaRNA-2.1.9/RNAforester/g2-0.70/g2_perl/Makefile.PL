use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	=> 'G2',
    'VERSION_FROM' => 'G2.pm', # finds $VERSION
    'LIBS'	=> ['-L./.. -lg2   -L/usr/local/lib -lm'],  # e.g., '-lm' 
    'DEFINE'	=> '-DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DLINUX=1 -DX_DISPLAY_MISSING=1 -DDO_PS=1 -DDO_FIG=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1',                   # e.g., '-DHAVE_SOMETHING' 
    'INC'	=> '-I./../include',           # e.g., '-I/usr/local/include' 
);
