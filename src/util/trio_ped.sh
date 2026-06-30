#!/bin/bash

# Shared helper for deriving trio sample IDs from a ped file, so the IDs stay
# in sync with the ped rather than being hardcoded in each pipeline script.
#
# Depends on the logging functions in src/util/logging.sh (source that first).

# read_trio_ids <ped_file>
#
# Validates the ped file contains exactly one trio record and sets the global
# variables kid_id, dad_id and mom_id from it.
#
# PED columns: family_id  individual_id(kid)  paternal_id(dad)  maternal_id(mom)  sex  phenotype
# https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets?tab=readme-ov-file#accessing-controlled-samples
read_trio_ids() {
    local trio_ped="$1"

    if [ ! -f "$trio_ped" ]; then
        log_error "ped file not found: $trio_ped"
        exit 1
    fi

    local n_records
    n_records=$(awk '!/^#/ && NF >= 4' "$trio_ped" | wc -l)
    if [ "$n_records" -ne 1 ]; then
        log_error "expected exactly one trio record in ${trio_ped}, found ${n_records}"
        exit 1
    fi

    read -r kid_id dad_id mom_id < <(awk '!/^#/ && NF >= 4 { print $2, $3, $4; exit }' "$trio_ped")
    if [ -z "$kid_id" ] || [ -z "$dad_id" ] || [ -z "$mom_id" ]; then
        log_error "could not extract kid/dad/mom IDs from ped: $trio_ped"
        exit 1
    fi

    log_info "Trio extracted from ${trio_ped}: kid=${kid_id} dad=${dad_id} mom=${mom_id}"
}
