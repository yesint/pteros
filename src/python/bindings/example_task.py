from pteros import *

class MyTask(Task_base):
    def pre_process(self):
        self.log.info("{} pre process".format(self.id))

    def process_frame(self,info):
        self.log.info("{} process_frame {}".format(self.id,info.valid_frame))

    def post_process(self,info):
        self.log.info("{} post_process {}".format(self.id,info.valid_frame))
