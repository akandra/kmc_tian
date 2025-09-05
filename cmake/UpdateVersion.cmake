# cmake/UpdateVersion.cmake
execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty --broken
                WORKING_DIRECTORY ${SOURCE_DIR}
                OUTPUT_VARIABLE GIT_VERSION_STRING
                OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
                WORKING_DIRECTORY ${SOURCE_DIR}
                OUTPUT_VARIABLE GIT_COMMIT_HASH
                OUTPUT_STRIP_TRAILING_WHITESPACE)

configure_file(${INPUT_FILE} ${OUTPUT_FILE} @ONLY)