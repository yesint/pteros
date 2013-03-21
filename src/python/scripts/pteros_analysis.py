import sys

sys.path.append("..")

from pteros import *
import pkgutil, copy

import pteros_analysis_plugins

print "+--------------------------------+"
print "+ This is pteros_analysis script +"
print "+--------------------------------+"

print "\nAvailable plugins:\n"
package = pteros_analysis_plugins
for importer, modname, ispkg in pkgutil.iter_modules(package.__path__):   
    #module = __import__(modname, fromlist="dummy")
    print "\t%s" % modname
    #print "Imported", module
    
# Now look for tasks in the command line

cmd = " ".join(sys.argv[1:])
print "\nProvided command line: '%s'" % cmd

opt = Options_tree()
opt.from_command_line(cmd)
print opt.to_json_string()

requested_tasks = set()
task_list = []
for task in opt.get_options("task"):
	task_name = task.get_value_string("")
	print "\t* Requested task '%s'" % task_name
	requested_tasks.add(task_name)
	task_list.append( task_name )
	

# Create processor class
class Processor(Trajectory_processor):
	def __init__(self,opt):
		Trajectory_processor.__init__(self,opt)
		# Init data
		self.task_list = []
		self.active_tasks = []		

	def pre_process(self):		
		for task in self.task_list:			
			# All task are active at the beginning		
			self.active_tasks.append(1)
			# We need to give each task a copy of system
			task.system = System(self.get_system())
			# Run pre_process for each task
			task.pre_process()
			
	def process_frame(self,info):
		for i in range(0,len(self.task_list)):
			if self.active_tasks[i] == 1:
				ret = self.task_list[i].process_frame(info)
				if ret == False:
					self.active_tasks[i] = 0
		if sum(self.active_tasks) > 0:
			return True
		else:
			return False
			
	def post_process(self,info):
		for task in self.task_list:				
			task.post_process(info)

# Create instance of processor
proc = Processor(opt)
# We need a container to keep all compiled tasks, otherwise they will call destructors
# and cause crash
compiled_list = []

# All pure python tasks are run sequencially on each frame	
# Try to load needed modules and connect the consumers
print "Loading needed plugins..."
for task in requested_tasks:
	module = __import__(pteros_analysis_plugins.__name__ + "." + task, fromlist="dummy")	
	if module.__file__.split('.')[-1][0:2] == "py":
		# This is pure Python plugin
		print "\t* Loaded pure Python plugin '%s'" % task
		
		# Now work with tasks		
		for tsk in task_list:
			if tsk == task:
				# create an independent instance of Task from that module
				obj = module.Task()
				# Add methods to the list of processor
				proc.task_list.append( obj )
	else:
		# This is compiled plugin
		print "\t* Loaded compiled plugin '%s'" % task

		# Now work with tasks		
		for tsk in task_list:
			if tsk == task:
				# create an independent instance of Task from that module
				obj = module.Task(proc,opt)
				compiled_list.append(obj) # Store it to prevent call of destructor!


for m in proc.task_list:
	print m


#t = Compiled_task(proc)
proc.run()
