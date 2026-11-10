# Shared architecture check; Python's standard library is sufficient.
if(BUILD_TESTING)
    find_program(ALGORITHM_AUDIT_PYTHON NAMES python3)
    if(ALGORITHM_AUDIT_PYTHON)
        add_test(NAME algorithm_independence
            COMMAND ${ALGORITHM_AUDIT_PYTHON} -B
                ${CMAKE_CURRENT_LIST_DIR}/check_algorithm_dependencies.py)
    endif()
endif()
