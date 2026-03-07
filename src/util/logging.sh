#!/bin/bash

# Global string to track which log files have been initialized (overwritten)
_INITIALIZED_LOG_FILES=":"

# Main logging function
log() {
    local log_file=$1
    local level=$2
    shift 2 
    local message="$@"
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    
    # Format the log message
    local log_string="[$timestamp] [$level] $message"
    
    # Check if we have initialized this specific log file yet
    if [[ "$_INITIALIZED_LOG_FILES" != *":$log_file:"* ]]; then
        # FIRST RUN: Use 'tee' without '-a' to overwrite the file
        if [[ "$level" == "ERROR" ]]; then
            echo "$log_string" | tee "$log_file" >&2
        else
            echo "$log_string" | tee "$log_file"
        fi
        
        # Add this file path to our tracking string
        _INITIALIZED_LOG_FILES="${_INITIALIZED_LOG_FILES}${log_file}:"
    else
        # SUBSEQUENT RUNS: Use 'tee -a' to append to the file
        if [[ "$level" == "ERROR" ]]; then
            echo "$log_string" | tee -a "$log_file" >&2
        else
            echo "$log_string" | tee -a "$log_file"
        fi
    fi
}

# Wrapper functions: Expecting <log_file> as $1, followed by the <message>
log_info()    { local f=$1; shift; log "$f" "INFO"    "$@"; }
log_warning() { local f=$1; shift; log "$f" "WARNING" "$@"; }
log_error()   { local f=$1; shift; log "$f" "ERROR"   "$@"; }
log_debug()   { local f=$1; shift; log "$f" "DEBUG"   "$@"; }
