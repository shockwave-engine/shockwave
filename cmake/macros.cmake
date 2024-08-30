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
macro(shockwave_coverage target)
    include(CodeCoverage)

    if(NOT MSVC)
        append_coverage_compiler_flags_to_target(YOUR_TARGET_NAME)

        setup_target_for_coverage_lcov(
            NAME ${target}-coverage EXECUTABLE ${target} EXCLUDE "test/*"
            "*/catch2/*"
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
