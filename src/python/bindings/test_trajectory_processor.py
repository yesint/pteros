import numpy, time
from pteros_py import *

opt = Options_tree()
opt.from_command_line("--trajectory[/media/data/semen/trajectories/grand_challenge/nowater.gro /media/data/semen/trajectories/grand_challenge/nowater.xtc --last_frame 100]")

class My_proc(Trajectory_processor):
	def pre_process(self):
		print "In Python pre_process"
	def post_process(self):
		print "In Python post_process"
	def process_frame(self,info):
		print "In Python process_frame: ", info.valid_frame
		return True

proc = My_proc(opt)
#proc = Trajectory_processor(opt)
proc.run()

#proc = Trajectory_processor1(opt)
#proc.run()
