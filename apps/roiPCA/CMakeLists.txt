# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/jagla/apps/roiPCA/CMakeLists.txt
#
# CMakeLists.txt file for roiPCA.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_jagla_apps_roiPCA)
message (STATUS "Configuring sandbox/jagla/apps/roiPCA")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (roiPCA linalg.cpp alglibmisc.cpp alglibinternal.cpp statistics.cpp specialfunctions.cpp solvers.cpp optimization.cpp interpolation.cpp integration.cpp fasttransforms.cpp diffequations.cpp dataanalysis.cpp ap.cpp linalg.h alglibmisc.h alglibinternal.h stdafx.h statistics.h specialfunctions.h solvers.h optimization.h interpolation.h integration.h fasttransforms.h diffequations.h dataanalysis.h ap.h roiPCA.cpp roiPCA.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (roiPCA ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (roiPCA)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install roiPCA in ${PREFIX}/bin directory
install (TARGETS roiPCA
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/roiPCA for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (roiPCA)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:roiPCA")
  set (CPACK_PACKAGE_NAME "roiPCA")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "roiPCA")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Bernd Jagla <bernd.jagla@pasteur.fr>")
  set (CPACK_PACKAGE_VENDOR "Bernd Jagla <bernd.jagla@pasteur.fr>")

  seqan_configure_cpack_app (roiPCA "roiPCA")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:roiPCA")

