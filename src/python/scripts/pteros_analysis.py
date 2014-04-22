#!/usr/bin/python

from pteros import *
import sys, pkgutil, copy, os, imp

import pteros_analysis_plugins

#--------------------------------------
# Create processor class
class Processor(Trajectory_processor):
        def __init__(self,opt):
                Trajectory_processor.__init__(self,opt)
                # Init data
                self.task_list = []

        def pre_process(self):
                for task in self.task_list:
                        # We need to give each task a copy of system
                        task.system = System(self.get_system())
                        # Run pre_process for each task
                        task.pre_process()

        def process_frame(self,info):
                for i in range(0,len(self.task_list)):
                        # We need to update frame 0 of each task with the current value
                        self.task_list[i].system.setFrame_data( self.get_frame_ptr() , 0)
                        self.task_list[i].process_frame(info)

        def post_process(self,info):
                for task in self.task_list:
                        task.post_process(info)
#--------------------------------------
def detailed_help(module):
    if module.__file__.split('.')[-1][0:2] == "py":
        obj = module.Task()
        s = "Type:\n\tPure python plugin\n"
    else:
        obj = module.Task(proc,opt)
        s = "Type:\n\tCompiled plugin\n"

    s += "File:\n\t%s" % module.__file__

    if hasattr(obj.__class__, "help"):
        print obj.help()
    print s;
#--------------------------------------
def general_help():
    print """
    Usage:
        pteros_analysis.py --trajectory[<options>] --task[name1 <options>] --task[name2 <options>] ...

        --help
            Display usage message
        --help traj
            Help for trajectory processing options
        --help plugins
            List available analysis plugins
        --help <plugin name>
            Detailed help for particular analysis plugin
        --help all
            Detailed help for all analysis plugins and trajectory processing options
    """
#--------------------------------------

print "+--------------------------------+"
print "+ This is pteros_analysis script +"
print "+--------------------------------+"

# Check if we have at least one command-line argument
if len(sys.argv)==1:
    # Show usage message
    general_help()
    sys.exit(0)


# Parse command line
opt,task_opts = parse_command_line(sys.argv,"task")

# Display help info
if opt.has("help"):
    # Create instance of processor
    proc = Processor(opt)

    # Check if we are asked for help for single plugin
    plugin = opt("help","all").as_string()
    if plugin == "traj":
        # Show trajectory processor options
        print proc.help()
    elif plugin != "all" and plugin != "plugins":
        # Load one specific plugin
        module = __import__(pteros_analysis_plugins.__name__ + "." +plugin, fromlist="dummy")
        print "\n[ %s ]" % plugin
        detailed_help(module)
    else:
        # Load all plugins
        # Check if asked for detailed help for each of them
        if plugin=="all":
            detailed = True
            general_help()
            print proc.help()
        else:
            detailed = False

        print "----------------------"
        print "+ Available plugins: +"
        print "----------------------"

        package = pteros_analysis_plugins
        for importer, modname, ispkg in pkgutil.iter_modules(package.__path__):
            module = __import__(pteros_analysis_plugins.__name__ + "." +modname, fromlist="dummy")
            print "[ %s ]" % modname
            if detailed:
                detailed_help(module)
                print "------------------------"
        if not detailed:
            print "For detailed help for particular plugin use --help <plugin_name>"
            print "For detailed help for all available plugins use --help all"

    sys.exit(0)


requested_tasks = []
task_list = []

unique_names = set()
# Process all task options
for task in task_opts:
        name = task.get_name()
	if name not in unique_names:
		unique_names.add(name)
		requested_tasks.append(task)
	

for task in task_opts:
        f = task("plugin_file","").as_string()
	if f:
		s = "from custom plugin_file '%s'" % f
	else:
		s = ""
        print "\t* Requested task '%s' %s" % (task.get_name(),s)
	
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
        task_name = task.get_name()

	# Task is loaded from the package pteros_analysis_plugins by default
	# If 'file' option is set for the task it is loaded from given file instead
        plugin_file = task("plugin_file","").as_string()
	
	if plugin_file:
		# Get full path
		full = os.path.abspath(plugin_file)
		# extract module name
		(mod_name,ext) = os.path.splitext(os.path.basename(full))
		# Append module search path
		sys.path.append(os.path.dirname(full))
		print "\t  Using custom plugin_file '%s'" % full
		module = __import__(mod_name, fromlist="dummy")
	else:
		module = __import__(pteros_analysis_plugins.__name__ + "." + task_name, fromlist="dummy")
					
	if module.__file__.split('.')[-1][0:2] == "py":
		# This is pure Python plugin
		print "\t* Loaded pure Python plugin '%s'" % task_name
		
		# Create internal consumer for python tasks
		proc.initialize();
		
		# Now work with tasks		
                for tsk in task_opts:
                        tsk_name = tsk.get_name()
			if tsk_name == task_name:
				# create an independent instance of Task from that module
				obj = module.Task()
				# Pass options to this task. For this create an options member in instance
				obj.options = task
				# Give this task a unique textual label
				obj.label = task_name + "_id" + str(task_num);	
				print "\t\t Created plugin instance '%s'" % obj.label			 
				# Add methods to the list of processor
				proc.task_list.append( obj )
				task_num+=1
	else:
		# This is compiled plugin
		print "\t* Loaded compiled plugin '%s'" % task_name

		# Now work with tasks		
                for tsk in task_opts:
                        tsk_name = tsk.get_name()
			if tsk_name == task_name:
				# create an independent instance of Task from that module
                                obj = module.Task(proc,tsk)
				# Give this task a unique textual label
				obj.label = task_name + "_id" + str(task_num);				 
				print "\t\t Created plugin instance '%s'" % obj.label			 
				compiled_list.append(obj) # Store it to prevent call of destructor!
				task_num+=1


print "+--------------------------------+"
print "+ Trajectory processor starts... +"
print "+--------------------------------+"
proc.run()
