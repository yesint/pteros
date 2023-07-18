from pteros import *

class Sample1(Task_base):
	def pre_process(self):
                sel_text = self.options("selection").as_string()
                self.sel = self.system(sel_text)
                self.use_mass = self.options("mass_weighted","false").as_bool()
                self.log.info("Working on selection '{}'".format(self.sel.get_text()))
                self.log.info("There are {} atoms in selection".format(self.sel.size()))
		
	def post_process(self,info):
                self.log.info("Finished!")
		
	def process_frame(self,info):
                self.log.info("Frame {} time {}".format(info.absolute_frame,info.absolute_time))
                self.log.info("Selection center: {}".format(self.sel.center(self.use_mass)))
