source bigdcd.tcl

set xtc_file /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_noPBC_1.xtc
set gro_file /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_pdb2gmx.gro

proc bench1 {frame} {
    global all ref
    $all frame $frame
    $all move [measure fit $all $ref]
    puts "$frame: [measure rmsd $all $ref]"
}

proc bench2 {frame} {
    global sel1 sel2
    $sel1 frame $frame
    $sel2 frame $frame
    set l [measure contacts 2.5 $sel1 $sel2]
    puts "frame $frame: [llength [lindex $l 0]]"
}

proc bench3 {frame} {
    global all
    set res [$all get resid]
    foreach r $res {
	set r [atomselect top "resid $r"]
	$r delete
    }
    puts "frame: $frame"
}

proc bench_all {frame} {
    bench1 $frame
    bench2 $frame
    bench3 $frame
}


mol new $gro_file waitfor all
mol addfile $xtc_file type xtc last 1

set all [atomselect top all]
set ref [atomselect top all frame 1]

set sel1 [atomselect top "resid 1 to 100"]
set sel2 [atomselect top "resid 102 to 200"]

bigdcd bench_all $xtc_file

bigdcd_wait

exit
