function(DUMP_VERSION_INFO PACKAGE_NAME VERSION_VAR)
  set(VERSION_STR "const char* ${PACKAGE_NAME}_VERSION=\"${VERSION_VAR}\";\n")
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/version.h "${VERSION_STR}")

endfunction()

function(DUMP_GIT_SUMMARY PACKAGE_NAME PKG_GIT_REPO)
    
  if(NOT "${PKG_GIT_REPO}" STREQUAL "")
    set(GIT_OPTIONS --git-dir=${PKG_GIT_REPO} --work-tree=${CMAKE_SOURCE_DIR})
  else()
    set(GIT_OPTIONS "")
  endif()

  execute_process(
                  COMMAND git ${GIT_OPTIONS} log --pretty=format:'%h' -n 1
                  OUTPUT_VARIABLE GIT_REV
                  ERROR_QUIET
                  )

  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/version.h "\n/* ##### ${PACKAGE_NAME} #####\n")

  # Check whether we got any revision (which isn't
  # always the case, e.g. when someone downloaded a zip
  # file from Github instead of a checkout)
  if ("${GIT_REV}" STREQUAL "")
      file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/version.h "// No git repo found for ${PACKAGE_NAME}\n")
      set(GIT_REV "N/A")
      set(GIT_DIFF "")
      set(GIT_TAG "N/A")
      set(GIT_BRANCH "N/A")
  else()
      #execute_process(
      #    COMMAND bash -c "git ${GIT_OPTIONS} diff --quiet --exit-code || echo +"
      #    OUTPUT_VARIABLE GIT_DIFF
      #    COMMAND_ECHO STDOUT)
      execute_process(
          COMMAND git  ${GIT_OPTIONS} describe --exact-match --tags
          OUTPUT_VARIABLE GIT_TAG ERROR_QUIET)
      execute_process(
          COMMAND git ${GIT_OPTIONS} rev-parse --abbrev-ref HEAD
          OUTPUT_VARIABLE GIT_BRANCH)

      string(STRIP "${GIT_REV}" GIT_REV)
      string(SUBSTRING "${GIT_REV}" 1 7 GIT_REV)
      #string(STRIP "${GIT_DIFF}" GIT_DIFF)
      string(STRIP "${GIT_TAG}" GIT_TAG)
      string(STRIP "${GIT_BRANCH}" GIT_BRANCH)
  endif()

  set(VERSION_STR "${PACKAGE_NAME} GIT HASH=\"${GIT_REV}\";
  ${PACKAGE_NAME} GIT TAG=\"${GIT_TAG}\";
  ${PACKAGE_NAME} GIT BRANCH=\"${GIT_BRANCH}\";*/ \n\n")

  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/version.h "${VERSION_STR}")
endfunction()

function(DUMP_CPU_INFO)
    execute_process(COMMAND lscpu OUTPUT_VARIABLE CPU_INFO)
    string(REPLACE "\n" "\n// " CPU_INFO_COMMENTED "${CPU_INFO}")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/version.h "\n/* ##### CPU INFO #####\n// ${CPU_INFO_COMMENTED}*/\n")
endfunction()

function(DUMP_GPU_INFO)
    if (MaCh3_GPU_ENABLED)
        execute_process(COMMAND nvidia-smi OUTPUT_VARIABLE GPU_INFO ERROR_QUIET)
        string(REPLACE "\n" "\n// " GPU_INFO_COMMENTED "${GPU_INFO}")
        file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/version.h "\n/* ##### GPU INFO #####\n// ${GPU_INFO_COMMENTED}*/\n")
    else()
        file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/version.h "\n/* ##### GPU INFO #####\n// No GPU information available\n*/\n")
    endif()
endfunction()

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/version.h "\n// Versions of dependencies used to build MaCh3\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/version.h "\nconst char* VERSION_HEADER_LOC=\"${CMAKE_CURRENT_BINARY_DIR}/version.h\";\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/version.h "const char* COMPILER=\"${CMAKE_CXX_COMPILER_ID}\";\n")

DUMP_VERSION_INFO(COMPILER "${CMAKE_CXX_COMPILER_VERSION}" )
DUMP_VERSION_INFO(ROOT "${ROOT_VERSION}")
DUMP_VERSION_INFO(MaCh3 "${MaCh3_VERSION}")
DUMP_VERSION_INFO(YAML_CPP "${YAML_CPP_VERSION}")
DUMP_VERSION_INFO(SPDLOG "${SPDLOG_VERSION}")

DUMP_GIT_SUMMARY(MaCh3 "")
DUMP_GIT_SUMMARY(YAML_CPP "${CMAKE_CURRENT_BINARY_DIR}/_deps/yaml-cpp-src/.git")
DUMP_GIT_SUMMARY(SPDLOG "${CMAKE_CURRENT_BINARY_DIR}/_deps/spdlog-src/.git")

DUMP_CPU_INFO()
DUMP_GPU_INFO()
