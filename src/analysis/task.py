sel = Selection(system)

def pre_process():
  sel.modify("name CA")
  print "In python pre_process:"
  print "There are",sel.size()," CA atoms"

def process_frame(info):
  print "In python process_frame at frame",info.valid_frame," center:",sel.center()

def post_process(info):
  print "In python post_process",info.valid_frame
