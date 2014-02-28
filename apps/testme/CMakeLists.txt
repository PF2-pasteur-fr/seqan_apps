# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/seqan_ngs_apps/apps/testme/CMakeLists.txt
#
# CMakeLists.txt file for testme.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_seqan_ngs_apps_apps_testme)
message (STATUS "Configuring sandbox/seqan_ngs_apps/apps/testme")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

if (NOT ZLIB_FOUND)
  message (STATUS "  ZLIB not found: RABEMA will not be built.")
  return ()
endif (NOT ZLIB_FOUND)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (testme testme.cpp testme.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (testme ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (testme)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install testme in ${PREFIX}/bin directory
install (TARGETS testme
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/testme for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (testme)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:testme")
  set (CPACK_PACKAGE_NAME "testme")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "testme")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Bernd Jagla <bernd.jagla@pasteur.fr>")
  set (CPACK_PACKAGE_VENDOR "Bernd Jagla <bernd.jagla@pasteur.fr>")

  seqan_configure_cpack_app (testme "testme")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:testme")

