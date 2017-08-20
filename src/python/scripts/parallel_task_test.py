from pteros import *
import sys
from multiprocessing import Process, cpu_count, Pipe, sharedctypes, Lock, Condition, Queue
import numpy as np
import traceback
from ctypes import Structure, c_int, c_float

class TestTask:
    def __init__(self):
        self.frames = []

    def pre_process(self):
      print 'pre:',self.id
      self.sel = self.system()
      self.result = []


    def process_frame(self,info):
      print 'proc:',self.id,'frame:', info.valid_frame, self.system.getXYZ(0,0)
      for i in xrange(2):
          self.sel.rotate(0,0.1)
      self.frames.append(info.valid_frame)
      self.result.append(info.valid_frame)

    def post_process(self,info):
      print 'post:',self.id,self.frames
      #self.result = self.id


def collector(results):
    for res in results:
        print res

#-----------------------------------------------------------------------

class _SharedInfoStruct(Structure):
    _fields_ = [('absolute_frame', c_int),
                ('absolute_time', c_float),
                ('valid_frame', c_int),
                ('first_time', c_float),
                ('last_time', c_float),
                ('first_frame', c_int),
                ('last_frame', c_int)]

    def set(self,info):
        self.absolute_frame = info.absolute_frame
        self.absolute_time = info.absolute_time
        self.valid_frame = info.valid_frame
        self.first_time = info.first_time
        self.first_frame = info.first_frame
        self.last_time = info.last_time
        self.last_frame = info.last_frame


    def get(self):
        info = Frame_info()
        info._set(
            (self.absolute_frame,
            self.absolute_time,
            self.valid_frame,
            self.first_time,
            self.last_time,
            self.first_frame,
            self.last_frame)
        )
        return info


def run_instance(instance,id,sbox,scoord,st,sinfo,cond,queue,res_queue,stop_val):

    instance.pre_process()
    while True:       
        # Tell reader that we are ready. Do not wait. If queue is full just ignore.
        try:
            queue.put_nowait(1)
        except:
            pass

        # Get the data lock
        cond.acquire()
        # Wait until notified that shared data are ready
        # (lock is released and reaquired when waiked up)
        cond.wait()

        stop = stop_val

        if not stop:
            # Recover shared data
            box = np.ctypeslib.as_array(sbox)
            coord = np.ctypeslib.as_array(scoord)
            t = st.value
            info = sinfo.get()

            # Convert shared arrays back to numpy
            frame = Frame()

            box.shape = (3,3)
            frame.box = box

            coord.shape = (instance.system.num_atoms(),3)
            frame.set_coord_array(coord)

            frame.t = t

        # Release lock
        cond.release()
        
        if not stop:
            instance.system.setFrame_data(frame,0)
            instance.process_frame(info)
        else:
            instance.post_process(info)
            # Send magic result variable back to caller
            res_queue.put(instance.result)
            break
    

class Parallel_dispatcher:
    def __init__(self):
        self.processes = []
        self.pipes = []        
        self.lock = Lock()
        self.cond = Condition(self.lock)
        self.queue = Queue(1)
        self.Ncpu = cpu_count()
        self.res_queue = Queue(self.Ncpu)
    
    def pre_process(self,system):
        try:
            print('\tRunning pure-python parallel task using {} processes'.format(self.Ncpu))

            # For Frame
            self.sbox = sharedctypes.RawArray('f',9)
            self.scoord = sharedctypes.RawArray('f',system.num_atoms()*3)
            self.st = sharedctypes.RawValue('f')
            # For Frame_info
            self.sinfo = sharedctypes.RawValue(_SharedInfoStruct)
            # Stop value
            self.stop = sharedctypes.RawValue('b',False)

            for i in range(self.Ncpu):
                task_instance = TestTask()
                self.system = system
                task_instance.system = System(system) # Make a copy
                task_instance.id = i
                task_instance.result = None

                self.processes.append( Process(target=run_instance, args=(
                        task_instance,
                        i,
                        self.sbox,
                        self.scoord,
                        self.st,
                        self.sinfo,
                        self.cond, self.queue, self.res_queue,
                        self.stop)) )

            for p in self.processes:
                p.start()

        except Exception as e:
            traceback.print_exc()
            raise
        

    def process_frame(self,info):             
        try:
            # Wait until some worker says it's ready
            self.queue.get()

            frame = self.system.getFrame_data(0)

            # Populate shared arrays
            box = frame.box
            box.shape = 9

            coord = frame.get_coord_array()
            coord.shape = 3*self.system.num_atoms()

            # Update shared variables
            self.cond.acquire()

            self.sbox[:] = box
            self.scoord[:] = coord
            self.st = frame.t
            self.sinfo.set(info)

            # Notify one waiting worker that shared data are updated
            self.cond.notify()
            self.cond.release()

        except Exception as e:
            traceback.print_exc()
            raise
               
        
    def post_process(self,info):
        try:
            self.cond.acquire()
            self.stop.value = True
            self.cond.release()

            for p in self.processes:
                # While process is alive waik it up until it accepts stop
                while p.is_alive():
                    self.cond.acquire()
                    self.cond.notify()
                    self.cond.release()

            results = []
            while not self.res_queue.empty():
                results.append( self.res_queue.get() )

            for p in self.processes:
                p.join()

            print('\tRunning aggregator handler...')
            collector(results)

        except Exception as e:
            traceback.print_exc()
            raise

if __name__ == '__main__':
    opt,task_opts = parse_command_line(sys.argv,"task")

    disp = Parallel_dispatcher()

    r = Trajectory_reader(opt, disp.pre_process, disp.process_frame, disp.post_process)
    r.run()
