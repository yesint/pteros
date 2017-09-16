mol new "/media/semen/data/semen/trajectories/2lao/after_em.gro"

set all0 [atomselect top all]
set all1 [atomselect top all]

set N 10000

animate dup top

$all0 frame 0
$all1 frame 1

set m [transaxis x 0.5 rad] 

proc p { } {
	global all0 all1 m
	$all1 move $m
	measure fit $all0 $all1 weight mass
}

time p $N