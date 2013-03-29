from pteros import *

class Task:		
	def pre_process(self):
		sel_text = self.options.get_value_string("selection")
		self.sel = Selection(self.system, sel_text )
		self.use_mass = self.options.get_value_bool("mass_weighted",False)
		print "Working on selection '%s'" % self.sel.get_text()
		print "There are %i atoms in selection" % self.sel.size()
		
	def post_process(self,info):
		print "Finished!"
		
	def process_frame(self,info):
		print "Frame ", info.absolute_frame, " time ", info.absolute_time
		print "Selection center: ", self.sel.center(self.use_mass)
		return True # Returning True means allow continuing with the next frame		
