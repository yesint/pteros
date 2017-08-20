#!/usr/bin/env python

from pteros import *
import sys, pkgutil, copy, os, imp, signal

#--------------------------------------
def general_help():
    print """Usage:
pteros_analysis.py -f <files> <processing options>... -task name1 <task1 options> -task name2 <task2 options> ...

-help traj
    Help for trajectory processing options
-help plugins
    List available analysis plugins
-help <plugin name>
    Detailed help for particular analysis plugin
-help all
    Detailed help for all analysis plugins and trajectory processing options
"""
#--------------------------------------

class Dispatcher:
    def __init__(self):
        self.task_list = []
        self.options = None
        self.label = None

    def pre_process(self,system):
        log.info("Pre-processing...")
        for task in self.task_list:
            # We need to give each task its own system
            task.system = System(system)
            # also save system locally
            self.system = system
            # We also need to give it its own jump remover
            task.jump_remover = Jump_remover()
            # Run pre_process for each task
            task.pre_process()
        log.info("Processing frames...")

    def process_frame(self,info):
        fr = self.system.getFrame_data(0)
        for task in self.task_list:
            # We need to update frame 0 of each task with the current value
            task.system.setFrame_data(fr, 0)
            # Call jump remover
            task.jump_remover.remove_jumps(task.system)
            # Process frame
            task.process_frame(info)

    def post_process(self,info):
        log.info("Post-processing...")
        for task in self.task_list:
            task.post_process(info)

disp = Dispatcher()

#--------------------------------------

if __name__ == '__main__':
    print "+--------------------------------+"
    print "+ This is pteros_analysis script +"
    print "+--------------------------------+"

    # Create logger
    log = Logger('analysis')

    # Check if we have at least one command-line argument:
    if len(sys.argv)==1:
        # Show usage message
        general_help()
        sys.exit(0)

    # Dirty hack - set special SIGINT handler, which kills the whole program
    # instead of raising KeyboardInterrupt in Python.
    # This allows to kill the process if it runs heavy compiled extension at background
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    # Parse command line
    opt,task_opts = parse_command_line(sys.argv,"task")

    if len(task_opts)==0:
        print
        general_help()
        sys.exit(0)

    # Create trajectory reader
    reader = Trajectory_reader(opt, disp.pre_process, disp.process_frame, disp.post_process)

    # If explicitly asked for help show it
    if opt.has("help"):
        topic = opt("help","all").as_string()
        if topic == "traj":
            # Show trajectory processing options
            print reader.help()

    # Load all supplied tasks
    files_to_load = set()

    log.info("Requested tasks:")
    for task in task_opts:
        f = task.get_name()
        log.info("\t%s" % task.get_name())
        if f not in files_to_load:
            files_to_load.add(f)

    task_num = 0
    log.info("Creating task instances:")
    for f in files_to_load:
        # Get full path
        full = os.path.abspath(f)
        # extract module name
        (mod_name,ext) = os.path.splitext(os.path.basename(full))
        # Append module search path
        sys.path.append(os.path.dirname(full))
        log.info("\tLoading '%s'" % full)
        module = __import__(mod_name, fromlist="dummy")

        # Create needed number of task instances
        for task in task_opts:
            if task.get_name() == f:
                # create new instance of Task from that module
                obj = module.Task()
                # Pass options to this task. For this create an options member in instance
                obj.options = task
                # Give this task a unique textual label
                obj.label = mod_name + "_id" + str(task_num);
                task_num += 1
                # Give it a logger
                obj.log = Logger(obj.label)
                log.info("\t\t Created task instance '%s'" % obj.label)
                # Add instance to dispatcher
                disp.task_list.append( obj )

    #--------------------------------------
    # RUN!
    reader.run()




