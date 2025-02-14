
#!/bin/bash

CLUBB_EXECUTABLE_FILE="../bin/clubb_standalone"

# Echo the case name
echo "Running $run_case"

# Run the CLUBB model
if [ -e $CLUBB_EXECUTABLE_FILE ]; then
    $CLUBB_EXECUTABLE_FILE
    CLUBB_EXIT_STATUS=$?
    if [ $CLUBB_EXIT_STATUS -eq 6 ]; then
        # The exit status of 6 is used as the success exit status in CLUBB.
        RESULT=0
    else
        RESULT=1
    fi
else
    echo "${CLUBB_EXECUTABLE_FILE} not found (did you re-compile?)"
    RESULT=1
fi