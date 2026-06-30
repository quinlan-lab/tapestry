#!/bin/bash

# Render trio_no_methylation_igv_session.xml from its template, substituting the
# trio sample IDs derived from the ped file. The rendered .xml is git-ignored;
# regenerate it whenever the trio (i.e. the ped file) changes.
#
# Usage:
#   Production:  ./generate_trio_igv_session.sh
#   Dev mode:    ./generate_trio_igv_session.sh --dev-dir trio_dev_data

source src/util/logging.sh
source src/util/trio_ped.sh

DEV_DIR=""

usage() {
    cat <<EOF
Usage: $0 [--dev-dir <DEV_DATA_DIR>]

Render trio_no_methylation_igv_session.xml from
trio_no_methylation_igv_session.template.xml, substituting the trio sample
IDs derived from the ped file.

Options:
  -d, --dev-dir <DEV_DATA_DIR>  Derive IDs from the dev ped instead of the
                                production ped.
  -h, --help                    Show this help message and exit.
EOF
}

# --- Argument Parsing ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dev-dir|-d)
            DEV_DIR="${2%/}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            log_error "Unknown parameter passed: $1"
            usage >&2
            exit 1
            ;;
    esac
done

# --- Configuration ---
trio_ped="trio.ped"  # production trio ped, committed to this repo
if [ -n "$DEV_DIR" ]; then
    trio_ped="${DEV_DIR}/input/trio.ped"
fi

template="trio_no_methylation_igv_session.template.xml"
rendered="trio_no_methylation_igv_session.xml"

if [ ! -f "$template" ]; then
    log_error "template not found: $template"
    exit 1
fi

# --- Derive trio sample IDs from the ped file ---
read_trio_ids "$trio_ped"

# --- Render ---
sed \
    -e "s/{{KID_ID}}/${kid_id}/g" \
    -e "s/{{DAD_ID}}/${dad_id}/g" \
    -e "s/{{MOM_ID}}/${mom_id}/g" \
    "$template" > "$rendered"

log_info "Rendered ${rendered} from ${template} (kid=${kid_id} dad=${dad_id} mom=${mom_id})"
