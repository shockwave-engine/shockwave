if(EXISTS ${FROM})
    message(STATUS "Copying ${FROM} to ${TO}")
    file(COPY_FILE ${FROM} ${TO} ONLY_IF_DIFFERENT)
else()
    message(STATUS "Not copying ${FROM} since it doesn't exist")
endif()
