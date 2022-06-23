# Module that checks whether the compiler supports the
# abi::__cxa_demangle function required to
# make the type names returned by typeid() human-readable
#
# Sets the following variable:
# HAVE_CXA_DEMANGLE
#
# perform tests
include(CheckCXXSourceCompiles)
include(AppendCompilerFlags)

CHECK_CXX_SOURCE_COMPILES("#include <cxxabi.h>
int main(void){
    int foobar = 0;
    const char *foo = typeid(foobar).name();
    int status;
    char *demangled = abi::__cxa_demangle( foo, 0, 0, &status );
}" HAVE_CXA_DEMANGLE)

include (FindPackageHandleStandardArgs)
# prevent useless message from being displayed
set (FIND_PACKAGE_MESSAGE_DETAILS_CxaDemangle "[1][v()]"
  CACHE INTERNAL "Details about finding CxaDemangle")
find_package_handle_standard_args (CxaDemangle DEFAULT_MSG HAVE_CXA_DEMANGLE)

if( HAVE_CXA_DEMANGLE )
    message(STATUS "${Green}Compiler support DEMANGLE${ColourReset}")
else()
    set(HAVE_CXA_DEMANGLE 0)
    message(STATUS "${Red}Compiler NOT support DEMANGLE${ColourReset}")
endif()