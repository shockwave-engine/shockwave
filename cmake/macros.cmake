macro(shockwave_set_compile_options target)
    if(MSVC)
        target_compile_options(
            ${target} PRIVATE /diagnostics:caret /nologo /W4 /WX /permissive-
        )
    else()
        target_compile_options(
            ${target}
            PRIVATE -Wall -Wextra -Wpedantic -Werror -Woverloaded-virtual
                    -Wno-unknown-pragmas -Wno-enum-compare -fsigned-char
        )
    endif()
endmacro()

# This is setup for linux only
# LCOV fails to detect unused methods defined in headers (inline)
macro(shockwave_coverage target)
    include(CodeCoverage)

    if(NOT MSVC)
        append_coverage_compiler_flags()
        # append_coverage_compiler_flags_to_target(${target})

        setup_target_for_coverage_lcov(
            NAME ${target}-coverage EXECUTABLE ${target} EXCLUDE "*/catch2/*"
            "/usr/include/*" "test/*"
        )

        add_custom_target(
            ${target}-coverage-viewer
            COMMAND xdg-open ${PROJECT_BINARY_DIR}/${target}-coverage/index.html
            DEPENDS ${target}-coverage
        )
    endif()
endmacro()

macro(shockwave_add_library target src...)
    set(src ${src...} ${ARGN})

    add_library(${target} ${src})
    shockwave_set_compile_options(${target})
    target_compile_definitions(${target} PRIVATE SHOCKWAVE_EXPORT)
endmacro()

macro(generates_compile_commands target)
    # clangd only supports a single CompilationDatabase location.
    # -> Copy the compile commands from out/[Debug | Release | ...] to out/
    add_custom_command(
        TARGET ${target}
        POST_BUILD
        COMMAND
            ${PROJECT_SOURCE_DIR}/scripts/copy-compile-commands.py
            "--from=${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json"
            "--to=${CMAKE_CURRENT_SOURCE_DIR}/out/compile_commands.json"
    )
endmacro()
