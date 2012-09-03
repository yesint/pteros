source bigdcd.tcl

set gro_file /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_pdb2gmx.gro
set xtc_file /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_noPBC_1.xtc

proc bench {num} {
    global sel gro_file id xtc_file
    $sel frame 0
    puts [time {
    for {set i 0} {$i<100000} {incr i} {
        switch $num {
	1 {
	    $sel moveby "0.1 0.1 0.1"
        }
        2 {
	    $sel move [transaxis x 0.1]
        }
        3 {
	    measure center $sel
        }
        4 {
	    measure minmax $sel
        }
        }
    }
    }]
}


set id [mol new $gro_file waitfor all]

set sel [atomselect top all]

bench 6

exit
