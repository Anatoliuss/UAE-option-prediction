# CMake generated Testfile for 
# Source directory: C:/Users/achuv/UAE-option-prediction
# Build directory: C:/Users/achuv/UAE-option-prediction/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  add_test([=[all_unit_tests]=] "C:/Users/achuv/UAE-option-prediction/build/Debug/tests.exe")
  set_tests_properties([=[all_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/achuv/UAE-option-prediction/CMakeLists.txt;32;add_test;C:/Users/achuv/UAE-option-prediction/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  add_test([=[all_unit_tests]=] "C:/Users/achuv/UAE-option-prediction/build/Release/tests.exe")
  set_tests_properties([=[all_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/achuv/UAE-option-prediction/CMakeLists.txt;32;add_test;C:/Users/achuv/UAE-option-prediction/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
  add_test([=[all_unit_tests]=] "C:/Users/achuv/UAE-option-prediction/build/MinSizeRel/tests.exe")
  set_tests_properties([=[all_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/achuv/UAE-option-prediction/CMakeLists.txt;32;add_test;C:/Users/achuv/UAE-option-prediction/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
  add_test([=[all_unit_tests]=] "C:/Users/achuv/UAE-option-prediction/build/RelWithDebInfo/tests.exe")
  set_tests_properties([=[all_unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/achuv/UAE-option-prediction/CMakeLists.txt;32;add_test;C:/Users/achuv/UAE-option-prediction/CMakeLists.txt;0;")
else()
  add_test([=[all_unit_tests]=] NOT_AVAILABLE)
endif()
subdirs("_deps/catch2-build")
