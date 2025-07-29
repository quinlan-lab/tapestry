
# TODO [by 30 July] 
# refactor dir variables in phase_meth_to_founder_haps.py into file paths and pass those as python-script params, and write a bash script (phase_meth_to_founder_haps.sh) that calls the python script
# check that Tom would be able to run phase_meth_to_founder_haps.sh using just the workflow described in the README 
# run phase_meth_to_founder_haps.sh on a pedigree, or at least a trio
#  - [ ] Update `build-iht-based-haplotype-map-and-phase-variants.sh` to:
    - Ensure `gtg-ped-map`, etc. are in the PATH
    - Run `bgzip`, etc.
# share repo with Tom N so that he can use it to phase Tom S's DNMs at TR loci to founder labels, and thereby compare with Tom S's current phasing methodology. 

# TODO [Aug: lab symposium] 
# [positive control] does my tool phase methylation levels to the correct parental haplotype at these imprinted loci: https://gemini.google.com/app/db2d1c5aa2a8ed10 ? 
# present phased methylation levels for a pedigree, 
# create a python library called "tapestry" for analyzing the phased meth data across generations, such as: given a haplotype in a child, has its methylation changed relative to the same haplotype in the parent; etc.
# and solicit feedback about other interesting questions to ask
