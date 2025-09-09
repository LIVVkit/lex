#!/bin/bash
CASE_OLD=${1}
CASE_NEW=${2}
cp -R ${CASE_OLD} ${CASE_NEW}
pushd ${CASE_NEW}
sed -i "s/${CASE_OLD}/${CASE_NEW}/g" *.yml
popd
