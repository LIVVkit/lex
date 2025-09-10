#!/bin/bash
CASE=$1
echo "LEX ON ${CASE}"

# Allow for standalone (outside of batch script) by setting WEBDIR if it's not already set
WEBDIR="${WEBDIR:-/global/cfs/projectdirs/e3sm/www/${USER}}"

# Run LIVVkit with all the configs we want for this case
# writes the output website to the user's scratch directory
mkdir -p ${SCRATCH}/lex

livv -V \
    config/${CASE}/cmb_gis.yml \
    config/${CASE}/smb_gis.yml \
    config/${CASE}/energy_e3sm_racmo_gis.yml \
    config/${CASE}/energy_e3sm_era5_gis.yml \
    config/${CASE}/energy_e3sm_merra2_merra_grid_gis.yml \
    config/${CASE}/energy_e3sm_ceres_gis.yml \
    config/${CASE}/cmb_ais.yml \
    config/${CASE}/smb_ais.yml \
    config/${CASE}/energy_e3sm_racmo_ais.yml \
    config/${CASE}/energy_e3sm_era5_ais.yml \
    config/${CASE}/energy_e3sm_merra2_merra_grid_ais.yml \
    config/${CASE}/energy_e3sm_ceres_ais.yml \
    -o $SCRATCH/lex/${CASE} &> livv_log_${CASE}.log

# Backup the existing published version of this analysis
mv ${WEBDIR}/${CASE} ${WEBDIR}/${CASE}_bkd_$(date +'%Y%m%dT%H%M')

# Move the new analysis from scratch to the published directory (prevents pre-maturely overwriting)
mv $SCRATCH/lex/${CASE} ${WEBDIR}

# Make the output directory rwxr-xr-x
chmod -R 0755 ${WEBDIR}/${CASE}

# Enter the Case output directory
pushd ${WEBDIR}/${CASE}

# Add the GROUP LINK back to the `current_runs.html` site, so all the "current" analyses are linked together
sed -i "s/\(<\!--GROUP LINK-->\)/\<a id=\"header-group\" href=\"..\/current_runs.html\">\nCurrent\ runs\n\<\/a\>/g" index.html
for htmlfile in validation/*.html
do
    sed -i "s/\(<\!-- GROUP LINK-->\)/\<a id=\"header-group\" href=\"..\/..\/current_runs.html\">\nCurrent\ runs\n\<\/a\>/g" ${htmlfile}
done
