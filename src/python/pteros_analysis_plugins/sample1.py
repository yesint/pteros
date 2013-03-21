from pteros import *

class Task:
	def pre_process(self):
		print "pre_process in sample1"
	def post_process(self,info):
		print "post_process in sample1"
	def process_frame(self,info):
		print "process_frame in sample1. Frame: ", info.valid_frame
