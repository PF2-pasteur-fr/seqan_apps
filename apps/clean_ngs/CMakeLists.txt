# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/seqan_ngs_apps/apps/clean_ngs/CMakeLists.txt
#
# CMakeLists.txt file for clean_ngs.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_seqan_ngs_apps_apps_clean_ngs)
message (STATUS "Configuring sandbox/seqan_ngs_apps/apps/clean_ngs")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (clean_ngs
                clean_ngs.cpp
                clean_ngs.h
                adaptor_record.h
                clean_ngs_app.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (clean_ngs ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (clean_ngs)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install clean_ngs in ${PREFIX}/bin directory
install (TARGETS clean_ngs
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/clean_ngs for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (clean_ngs)

# ----------------------------------------------------------------------------
# Setup Command Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

#Include executable pair_align in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} clean_ngs CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:clean_ngs")
  set (CPACK_PACKAGE_NAME "clean_ngs")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "clean_ngs")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Bernd Jagla <bernd.jagla@pasteur.fr>")
  set (CPACK_PACKAGE_VENDOR "Bernd Jagla <bernd.jagla@pasteur.fr>")

  seqan_configure_cpack_app (clean_ngs "clean_ngs")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:clean_ngs")

