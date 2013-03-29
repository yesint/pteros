#!/bin/bash
python pteros_analysis.py \
--trajectory[ \
 /media/data/semen/trajectories/grand_challenge/nowater.gro \
 /media/data/semen/trajectories/grand_challenge/nowater.xtc \
 --first_frame 0 --last_frame 10 \
] \
--task center \
--task center \
--task [sample1 --selection "all" --mass_weighted true]\
--task [sample2 --val 2]\
--task [user_script --plugin_file script1.py --val 42]
