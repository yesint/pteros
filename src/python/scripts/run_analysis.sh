#!/bin/bash
python pteros_analysis.py \
--trajectory[ \
 /media/data/semen/trajectories/grand_challenge/nowater.gro \
 /media/data/semen/trajectories/grand_challenge/nowater.xtc \
 --first_frame 0 --last_frame 10 \
] \
--task [sample1 --selection "all" --mass_weighted true]\
--task [sample1 --selection "all" --mass_weighted true]\
