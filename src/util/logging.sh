#!/bin/bash

# Main logging function — logs to stdout (INFO/DEBUG/WARNING) or stderr (ERROR)
log() {
    local level=$1
    shift
    local message="$@"
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    local log_string="[$timestamp] [$level] $message"

    if [[ "$level" == "ERROR" ]]; then
        echo "$log_string" >&2
    else
        echo "$log_string"
    fi
}

# Wrapper functions: Expecting <message> as arguments
log_info()    { log "INFO"    "$@"; }
log_warning() { log "WARNING" "$@"; }
log_error()   { log "ERROR"   "$@"; }
log_debug()   { log "DEBUG"   "$@"; }
