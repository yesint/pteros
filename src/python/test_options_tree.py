import numpy
from pteros_py import *

o = Options_tree()
o.from_command_line("--trajectory [a.xtc b.pdb c.xtc --last_frame 100] --dump --kill")
print o.to_command_line()
s =  o.to_json_string()
print s

o1 = Options_tree()
o1.from_json_string(s)
print o1.to_json_string()

opt = o.get_option("trajectory")
print opt.get_values_string("")

print o.get_value_int("trajectory/last_frame")
print o.get_value_double("trajectory/last_frame_nonexist",3.14)

