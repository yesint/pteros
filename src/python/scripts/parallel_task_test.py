from pteros import *
import sys
from multiprocessing import Process, cpu_count, Pipe, sharedctypes, Lock
import numpy as np

class TestTask:
    def pre_process(self):
      print 'pre:',self.id
      self.sel = self.system()
      self.frames = []

    def process_frame(self,info):
      print 'proc:',self.id, info.valid_frame, self.system.getXYZ(0,0)
      for i in xrange(100):
          self.sel.rotate(0,0.1)
      self.frames.append(info.valid_frame)

    def post_process(self,info):
      print 'post:',self.id,self.frames


def run_instance(instance,id,pipe,sbox,scoord,lock):
    instance.pre_process()
    while True:
        # Signal that we are ready to get new frame
        pipe.send(1)
        
        # Get responce with data
        info, t, stop = pipe.recv()

        # Convert shared arrays back to numpy
        frame = Frame()
        
        lock.acquire()
        box = np.ctypeslib.as_array(sbox)
        coord = np.ctypeslib.as_array(scoord)
        lock.release()
        
        box.shape = (3,3)
        frame.box = box        
        
        coord.shape = (instance.system.num_atoms(),3)
        frame.set_coord_array(coord)
        
        frame.t = t
        
        # push frame to the system of instance
        instance.system.setFrame_data(frame,0)
        
        if not stop:
            instance.process_frame(info)
        else:
            instance.post_process(info)
            break
    

class Parallel_dispatcher:
    def __init__(self):
        self.processes = []
        self.pipes = []        
        self.lock = Lock()
        self.instances = []
    
    def pre_process(self,system):
        
        self.sbox = sharedctypes.RawArray('f',9)        
        self.scoord = sharedctypes.RawArray('f',system.num_atoms()*3)
        
        for i in range(cpu_count()):
            task_instance = TestTask()
            self.system = system
            task_instance.system = System(system) # Make a copy
            task_instance.id = i
            self.instances.append(task_instance)
            end1, end2 = Pipe()
            self.processes.append( Process(target=run_instance, args=(task_instance,i,end2,self.sbox,self.scoord,self.lock)) )
            self.pipes.append( end1 )
        for p in self.processes:
            p.start()
        

    def process_frame(self,info):             
        frame = self.system.getFrame_data(0)
        # Populate shared arrays
        box = frame.box
        box.shape = 9
                        
        coord = frame.get_coord_array()
        coord.shape = 3*self.system.num_atoms()
        
        self.lock.acquire()
        self.sbox[:] = box
        self.scoord[:] = coord
        self.lock.release()
        
        # Send data to first ready pipe
        done = False
        while not done:
            for pipe in self.pipes:
                if pipe.poll():
                    pipe.recv() # Remove ready tag
                    pipe.send( (info, frame.t, False) )
                    done = True
                    break
        
        
    def post_process(self,info):
        for pipe in self.pipes:            
            pipe.send( (info,0,True) )
        for p in self.processes:
            p.join()


if __name__ == '__main__':
    opt,task_opts = parse_command_line(sys.argv,"task")

    disp = Parallel_dispatcher()

    r = Trajectory_reader(opt, disp.pre_process, disp.process_frame, disp.post_process)
