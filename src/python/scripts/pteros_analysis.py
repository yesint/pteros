from pteros import *
import sys, pkgutil, copy

import pteros_analysis_plugins

print "+--------------------------------+"
print "+ This is pteros_analysis script +"
print "+--------------------------------+"

#---------------------------------
def available_plugins():
	print "\nAvailable plugins:\n"
	package = pteros_analysis_plugins
	for importer, modname, ispkg in pkgutil.iter_modules(package.__path__):   
		#module = __import__(modname, fromlist="dummy")
		print "\t%s" % modname
		#print "Imported", module
#---------------------------------    

# Now look for tasks in the command line

#cmd = " ".join(sys.argv[1:])
#print "\nProvided command line: '%s'" % cmd

opt = Options_tree()
opt.from_command_line(sys.argv)
print opt.to_json_string()

requested_tasks = []
task_list = []

unique_names = set()
for task in opt.get_options("task"):
	name = task.get_value_string("")
	if name not in unique_names:
		unique_names.add(name)
		requested_tasks.append(task)
		task_list.append(task)

for task in requested_tasks:
	f = task.get_value_string("plugin_file","")
	if f:
		s = "from file '%s'" % f
	else:
		s = ""
	print "\t* Requested task '%s' %s" % (task.get_value_string(""),s)
	
#--------------------------------------
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
				# We need to update frame 0 of each task with the current value
				self.task_list[i].system.setFrame_data( self.get_system().getFrame_data(0), 0)
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
#--------------------------------------

# Create instance of processor
proc = Processor(opt)
# We need a container to keep all compiled tasks, otherwise they will call destructors
# and cause a crash
compiled_list = []

# All pure python tasks are run sequencially on each frame	

# Try to load needed modules and connect the consumers
task_num = 0
print "Loading needed plugins..."
for task in requested_tasks:
	task_name = task.get_value_string("")

	# Task is loaded from the package pteros_analysis_plugins by default
	# If 'file' option is set for the task it is loaded from given file instead
	plugin_file = task.get_value_string("plugin_file","")
	
	if plugin_file:
		module = __import__(plugin_file[:-2], fromlist="dummy")	
	else:
		module = __import__(pteros_analysis_plugins.__name__ + "." + task_name, fromlist="dummy")
					
	if module.__file__.split('.')[-1][0:2] == "py":
		# This is pure Python plugin
		print "\t* Loaded pure Python plugin '%s'" % task_name
		
		# Now work with tasks		
		for tsk in task_list:
			if tsk == task:
				# create an independent instance of Task from that module
				obj = module.Task()
				# Pass options to this task. For this create an options member in instance
				obj.options = task
				# Give this task a unique textual label
				obj.label = task_name + "_id" + str(task_num);				 
				# Add methods to the list of processor
				proc.task_list.append( obj )
				task_num+=1
	else:
		# This is compiled plugin
		print "\t* Loaded compiled plugin '%s'" % task_name

		# Now work with tasks		
		for tsk in task_list:
			if tsk == task:
				# create an independent instance of Task from that module
				obj = module.Task(proc,task)
				# Give this task a unique textual label
				obj.label = task_name + "_id" + str(task_num);				 
				compiled_list.append(obj) # Store it to prevent call of destructor!
				task_num+=1


for m in proc.task_list:
	print m


#t = Compiled_task(proc)
proc.run()
