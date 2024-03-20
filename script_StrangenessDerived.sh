#!/bin/bash 

export OPTIONS="-b --configuration json://${PWD}/jsons/json_StrangenessDerived.json --resources-monitoring 2 --aod-memory-rate-limit 1000000000 --shm-segment-size 6000000000"
#export OPTIONS="-b --configuration json://${PWD}/jsons/configuration_testj.json --resources-monitoring 2 --aod-memory-rate-limit 1000000000 --shm-segment-size 6000000000"
#o2-analysis-timestamp ${OPTIONS} \
#| o2-analysis-event-selection ${OPTIONS} \
#| o2-analysis-lf-lambdakzerobuilder ${OPTIONS} \
#| o2-analysis-lf-cascadebuilder ${OPTIONS} \
#| o2-analysis-lf-strangederivedbuilder ${OPTIONS} \
#| o2-analysis-pid-tpc-base ${OPTIONS} \
#| o2-analysis-pid-tpc ${OPTIONS} \
#| o2-analysis-pid-tof-base ${OPTIONS} \
#| o2-analysis-lf-cascadeflow ${OPTIONS}
#| o2-analysis-track-propagation ${OPTIONS} \
#| o2-analysis-track-extra ${OPTIONS} \
o2-analysis-lf-cascadespawner ${OPTIONS} \
| o2-analysis-lf-cascadeflow ${OPTIONS} --aod-file AO2D.root --aod-writer-keep dangling
#| o2-analysis-lf-derivedcascadeanalysis ${OPTIONS}









