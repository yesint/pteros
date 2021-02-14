from pteros import *

class Task(Task_base):
    def pre_process(self):
        self.log.info("{} pre process".format(self.label))

    def process_frame(self,info):
        self.log.info("{} process_frame {}".format(self.label,info.valid_frame))

    def post_process(self,info):
        self.log.info("{} post_process {}".format(self.label,info.valid_frame))
