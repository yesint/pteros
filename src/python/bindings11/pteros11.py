from _pteros11 import *


class Task_base(Task_base_):
    def __init__(self,opt):
        Task_base_.__init__(self,opt)
        # Set class name for logger initialization on C++ side
        self._class_name = type(self).__name__

