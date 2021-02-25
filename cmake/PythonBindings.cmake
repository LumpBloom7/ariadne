find_package(PythonLibs)

if(PYTHONLIBS_FOUND)
    if(NOT EXISTS ${PROJECT_SOURCE_DIR}/python/pybind11/CMakeLists.txt)
        message(STATUS "pybind11 dependency not available.")
        find_package(Git)
        if(GIT_FOUND)
            message(STATUS "Downloading pybind11 repository...")
            if (NOT EXISTS ${PROJECT_SOURCE_DIR}/.git) # Manages the case when an archive is used
                execute_process(COMMAND git init WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
                execute_process(COMMAND rm -Rf ${PROJECT_SOURCE_DIR}/python/pybind11 OUTPUT_QUIET ERROR_QUIET)
                execute_process(COMMAND git submodule add -f https://github.com/pybind/pybind11 python/pybind11 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
                execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
            else() # When using a clone
                execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
            endif()
            if(NOT EXISTS ${PROJECT_SOURCE_DIR}/python/pybind11/CMakeLists.txt)
                message(ERROR "pybind11 download failure.")
            else()
                message(STATUS "Downloaded pybind11 dependency successfully.")
            endif()
        else()
            message(ERROR "Git not found, pybind11 dependency could not be downloaded.")
        endif()
    endif()

    add_subdirectory(python/pybind11)
    add_subdirectory(python/source)

    add_library(pyariadne MODULE
        $<TARGET_OBJECTS:pyariadne-module-obj>
        $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
        $<TARGET_OBJECTS:pyariadne-algebra-submodule-obj>
        $<TARGET_OBJECTS:pyariadne-extra-submodule-obj>
    )
    set_target_properties(pyariadne PROPERTIES PREFIX "" OUTPUT_NAME "pyariadne" SUFFIX ".so")
    target_link_libraries (pyariadne PUBLIC ariadne ${GCOV_LIBRARIES} PRIVATE pybind11::module)

    add_library(pyariadne-numeric MODULE
        $<TARGET_OBJECTS:pyariadne-numeric-module-obj>
        $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
        $<TARGET_OBJECTS:ariadne-utility>
        $<TARGET_OBJECTS:ariadne-numeric>
    )
    set_target_properties(pyariadne-numeric PROPERTIES PREFIX "" OUTPUT_NAME "numeric" SUFFIX ".so")
    target_link_libraries (pyariadne-numeric PUBLIC ${GCOV_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES} PRIVATE pybind11::module)

    add_library(pyariadne-algebra MODULE
        $<TARGET_OBJECTS:pyariadne-algebra-module-obj>
        $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
        $<TARGET_OBJECTS:pyariadne-algebra-submodule-obj>
    )
    set_target_properties(pyariadne-algebra PROPERTIES PREFIX "" OUTPUT_NAME "algebra" SUFFIX ".so")
    target_link_libraries (pyariadne-algebra PUBLIC ${GCOV_LIBRARIES} ariadne PRIVATE pybind11::module)

    execute_process(COMMAND ${PYTHON_EXECUTABLE} -m site --user-site OUTPUT_VARIABLE PYARIADNE_INSTALL_DIR)
    string(STRIP ${PYARIADNE_INSTALL_DIR} PYARIADNE_INSTALL_DIR)
    install(TARGETS pyariadne DESTINATION ${PYARIADNE_INSTALL_DIR})

    add_subdirectory(python/tests)

else()
    message(STATUS "PythonLibs not found: Python bindings will not be built.")
endif()