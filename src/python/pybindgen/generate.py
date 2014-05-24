#!/usr/bin/python

from pybindgen import *
from pybindgen.typehandlers.base import ForwardWrapperBase

import sys

IN = Parameter.DIRECTION_OUT
OUT = Parameter.DIRECTION_OUT

class EigenParam(Parameter):

    DIRECTIONS = []
    CTYPES = []
    PYARRAY_DTYPE = ''
    CPP_DTYPE = ''
    NDIM = 1
    DIM1 = -1
    DIM2 = -1
    
    def __init__(self, ctype, name, direction=1):
        super(EigenParam,self).__init__(ctype, name, direction)
      
    def convert_python_to_c(self, wrapper):
        assert isinstance(wrapper, ForwardWrapperBase)
        
        dim1=''
        dim2=''
        
        if self.direction & Parameter.DIRECTION_IN:
            # Declare pyobjects            
            arg = wrapper.declarations.declare_variable("PyObject*", self.name+'_arg','NULL')
            arr = wrapper.declarations.declare_variable("PyObject*", self.name+'_arr','NULL')
            # Parse pyobject from python params
            wrapper.parse_params.add_parameter('O', ['&'+arg], self.value)
            # Now we need to create an array from obtained object        
            wrapper.before_call.write_code('%s = PyArray_FROM_OTF(%s, %s, NPY_IN_ARRAY);' % (arr,arg,self.PYARRAY_DTYPE))
            # Add cleanup for this temporary array
            wrapper.before_call.add_cleanup_code("Py_XDECREF(%s);" % arr)
            # Check if conversion is OK
            wrapper.before_call.write_error_check('%s==NULL' % arr, 'PyErr_SetString(PyExc_TypeError,"Convertion to PyArray failed!");')
            # Now check number of dimensions
            wrapper.before_call.write_error_check('PyArray_NDIM(%s)!=%i' % (arr,self.NDIM), 
                                'PyErr_SetString(PyExc_TypeError,"%iD array expected!");' % self.NDIM)
            # Now check sizes for static arrays
            if self.DIM1!=-1 and self.DIM2!=-1:
                if self.NDIM==1:
                    wrapper.before_call.write_error_check('PyArray_DIM(%s,0)!=%i' % (arr,self.DIM1), 
                                    'PyErr_SetString(PyExc_RuntimeError,"Need only %i elements in array!");' % self.DIM1)
                    dim1 = "%i" % self.DIM1
                else:  
                    wrapper.before_call.write_error_check('PyArray_DIM(%s,0)!=%i' % (arr,self.DIM1), 
                                    'PyErr_SetString(PyExc_RuntimeError,"Need only %i elements in dimension 1!");' % self.DIM1)
                    wrapper.before_call.write_error_check('PyArray_DIM(%s,1)!=%i' % (arr,self.DIM2), 
                                    'PyErr_SetString(PyExc_RuntimeError,"Need only %i elements in dimension 2!");' % self.DIM2)
                    dim1 = "%i" % self.DIM1
                    dim2 = "%i" % self.DIM2
            else:
                # We have dynamic-size array. Extract its dimensions form pyarray
                if self.NDIM==1:
                    dim1 = 'PyArray_DIM(%s,0)' % arr
                else:
                    dim1 = 'PyArray_DIM(%s,0)' % arr
                    dim2 = 'PyArray_DIM(%s,1)' % arr
           
    
            # Everything is Ok, do a mapping
            # 1D mapping
            if self.NDIM==1:
                wrapper.before_call.write_code('Eigen::Map<%s> %s_map((%s*)PyArray_DATA(%s),%s,1);' %
                                                    (self.CPP_DTYPE,self.name,self.CAST_DTYPE,arr,dim1))
            else: # 2D mapping
                wrapper.before_call.write_code('Eigen::Map<%s> %s_map((%s*)PyArray_DATA(%s),%s,%s);' %
                                                    (self.CPP_DTYPE,self.name,self.CAST_DTYPE,arr,dim2,dim1))
                
            wrapper.call_params.append(self.name+'_map')    

        if self.direction & Parameter.DIRECTION_OUT:
            # We don't need to parse this parameter at all ince its out and will be converted to return tuple element

            assert self.DIM1!=-1 and self.DIM2!=-1  # Dynamic size is not allowed for out params!

            # This is ordinary static case
            dim1 = "%i" % self.DIM1
            dim2 = "%i" % self.DIM2
                    
            # Create mapped array
            if self.NDIM==1:
                wrapper.declarations.declare_variable('npy_intp','%s_sz[1]={%s}' % (self.name,dim1))
                wrapper.declarations.declare_variable('PyObject*','%s_to_map' % self.name, 'PyArray_SimpleNew(1, %s_sz, %s)' %
                                                            (self.name,self.PYARRAY_DTYPE))
                wrapper.declarations.declare_variable('Eigen::Map<%s>' % self.CPP_DTYPE,'%s_call((%s*)PyArray_DATA(%s_to_map),%s,1)' %
                                                            (self.name,self.CAST_DTYPE,self.name,dim1))
            else: # 2D case
                wrapper.declarations.declare_variable('npy_intp','%s_sz[2]={%s,%s}' % (self.name,dim1,dim2))
                wrapper.declarations.declare_variable('PyObject*','%s_to_map' % self.name, 'PyArray_SimpleNew(2, %s_sz, %s)' %
                                                            (self.name,self.PYARRAY_DTYPE))
                wrapper.declarations.declare_variable('Eigen::Map<%s>' % self.CPP_DTYPE,'%s_call((%s*)PyArray_DATA(%s_to_map),%s,%s)' %
                                                            (self.name,self.CAST_DTYPE,self.name,dim2,dim1))                
                                                                            
            # Add call param
            wrapper.call_params.append(self.name+'_call')
            # And return the param in output tuple
            wrapper.build_params.add_parameter("N", [self.name+'_to_map'])  


class Vector3fParam(EigenParam):
    DIRECTIONS = [Parameter.DIRECTION_IN, Parameter.DIRECTION_OUT]
    CTYPES = ['Eigen::Vector3f&','Vector3f_const_ref','Vector3f_ref']
    PYARRAY_DTYPE = 'NPY_FLOAT'
    CPP_DTYPE = 'Eigen::Vector3f'
    CAST_DTYPE = 'float'
    NDIM = 1
    DIM1 = 3
    DIM2 = 1

class Vector3iParam(EigenParam):
    DIRECTIONS = [Parameter.DIRECTION_IN, Parameter.DIRECTION_OUT]
    CTYPES = ['Eigen::Vector3i&','Vector3i_const_ref','Vector3i_ref']
    PYARRAY_DTYPE = 'NPY_INT'
    CPP_DTYPE = 'Eigen::Vector3i'
    CAST_DTYPE = 'int'
    NDIM = 1
    DIM1 = 3
    DIM1 = 1
    
class Matrix3fParam(EigenParam):
    DIRECTIONS = [Parameter.DIRECTION_IN, Parameter.DIRECTION_OUT]
    CTYPES = ['Eigen::Matrix3f&','Matrix3f_const_ref','Matrix3f_ref']
    PYARRAY_DTYPE = 'NPY_FLOAT'
    CPP_DTYPE = 'Eigen::Matrix3f'
    CAST_DTYPE = 'float'
    NDIM = 2
    DIM1 = 3
    DIM2 = 3
    
class MatrixXfParam(EigenParam):
    DIRECTIONS = [Parameter.DIRECTION_IN]
    CTYPES = ['Eigen::MatrixXf&','MatrixXf_const_ref','MatrixXf_ref']
    PYARRAY_DTYPE = 'NPY_FLOAT'
    CPP_DTYPE = 'Eigen::MatrixXf'
    CAST_DTYPE = 'float'
    NDIM = 2
    DIM1 = -1
    DIM2 = -1
              
#----------------------------------------------------------------------------------------------------------

class Vector3fReturn(ReturnValue):
    CTYPES = ['Eigen::Vector3f','Vector3f']
    
    def __init__(self, ctype):
        saved = ctype.ctype
        super(Vector3fReturn, self).__init__("#define _a_hack_to_disable_retval_generation")
        self.ctype = saved

    def convert_c_to_python(self, wrapper):
        wrapper.declarations.declare_variable('npy_intp','retval_sz[1] = {3}')
        wrapper.declarations.declare_variable('PyObject*','retval_to_map = PyArray_SimpleNew(1, retval_sz, NPY_FLOAT)')
        wrapper.declarations.declare_variable('Eigen::Map<Eigen::Vector3f>','retval((float*)PyArray_DATA(retval_to_map),3,1)')
        wrapper.build_params.add_parameter("N", ['retval_to_map'],prepend=True)

class Matrix3fReturn(ReturnValue):
    CTYPES = ['Eigen::Matrix3f','Matrix3f']
    
    def __init__(self, ctype):
        saved = ctype.ctype
        super(Matrix3fReturn, self).__init__("#define _a_hack_to_disable_retval_generation")
        self.ctype = saved

    def convert_c_to_python(self, wrapper):
        wrapper.declarations.declare_variable('npy_intp','retval_sz[2] = {3,3}')
        wrapper.declarations.declare_variable('PyObject*','retval_to_map = PyArray_SimpleNew(2, retval_sz, NPY_FLOAT)')
        wrapper.declarations.declare_variable('Eigen::Map<Eigen::Matrix3f>','retval((float*)PyArray_DATA(retval_to_map),3,3)')
        wrapper.build_params.add_parameter("N", ['retval_to_map'],prepend=True)


class MatrixXfReturn(ReturnValue):
    CTYPES = ['Eigen::MatrixXf','MatrixXf',]
    
    def __init__(self, ctype, dim1=None, dim2=None):
        saved = ctype.ctype
        super(MatrixXfReturn, self).__init__("#define _a_hack_to_disable_retval_generation")
        self.ctype = saved
        self.dim1_code = dim1
        self.dim2_code = dim2

    def convert_c_to_python(self, wrapper):
        assert dim1!=None and dim2!=None
        wrapper.declarations.declare_variable('npy_intp','retval_sz[2] = {%s,%s}' %
                                                 (self.dim1_code,self.dim2_code))
        wrapper.declarations.declare_variable('PyObject*','retval_to_map = PyArray_SimpleNew(2, retval_sz, NPY_FLOAT)')
        wrapper.declarations.declare_variable('Eigen::Map<Eigen::Matrix3f>','retval((float*)PyArray_DATA(retval_to_map),retval_sz[1],retval_sz[0])')
        wrapper.build_params.add_parameter("N", ['retval_to_map'],prepend=True)

#-----------------------------------------------------

# Callback parameter for System::load
class load_callback_Param(Parameter):
    DIRECTIONS = [Parameter.DIRECTION_IN]
    CTYPES = ['load_callback']

    def convert_python_to_c(self,wrapper):
        assert isinstance(wrapper, ForwardWrapperBase)
        py_cb = wrapper.declarations.declare_variable("PyObject*", self.name,"NULL")
        wrapper.parse_params.add_parameter('O', ['&'+py_cb], self.name, optional=True)        
        
        wrapper.before_call.write_error_check("%s && !PyCallable_Check(%s)" % (py_cb,py_cb),
                                              """PyErr_SetString(PyExc_TypeError, "callback must be callable");""")
        wrapper.before_call.write_code("Py_XINCREF(%s);" % py_cb)
        # Make a lambda, which will excute callback
        wrapper.before_call.write_code("""
std::function<bool(pteros::System*,int)> %s_lam=nullptr;
if(%s){
    %s_lam = [&](pteros::System* sys, int fr)->bool {        
        PyObject* ret = PyObject_CallFunction(%s, (char*) "Oi", self, fr);
        bool res = (ret == Py_True) ? true : false;
        Py_DECREF(ret);
        return res;
    };
}
        """ % (self.name,py_cb,self.name,py_cb))        
        wrapper.before_call.add_cleanup_code("Py_XDECREF(%s);" % py_cb)
        wrapper.call_params.append(self.name+'_lam')                                              

#==================================================================================

# Support for l-value returns

def lvalue_return(mod,cl,method_name,retval,params):
    class_name = mod.cpp_namespace+'::'+cl.name
    # Make decl string
    params_decl_str=''
    for p in params:
        params_decl_str += p[0]+' '+p[1]+','
    params_decl_str=params_decl_str[:-1]
    # Make call string
    params_call_str=''
    for p in params:
        params_call_str += p[1]+','
    params_call_str=params_call_str[:-1]
    # Getter
    mod.header.writeln("// Wrappers for lvalue-returning method %s::%s" % (class_name,method_name))
    mod.header.writeln("namespace %s {" % mod.cpp_namespace)
    mod.header.writeln("%s get_%s(%s& obj, %s){" % (retval[0],method_name,class_name,params_decl_str))
    mod.header.writeln("  return obj.%s(%s);" % (method_name,params_call_str))
    mod.header.writeln("}")
    # Setter
    if retval[0]!='int' and retval[0]!='float':
        mod.header.writeln("void set_%s(%s& obj, %s, const %s& data){" % 
                        (method_name,class_name,params_decl_str,retval[0]))
    else:
        mod.header.writeln("void set_%s(%s& obj, %s, %s data){" % 
                        (method_name,class_name,params_decl_str,retval[0]))
    mod.header.writeln("  obj.%s(%s) = data;" % (method_name,params_call_str))
    mod.header.writeln("}")
    mod.header.writeln("} // namespace")
    # Add wrappers
    getter_params = [param(class_name,'obj')]+params
    if retval[0]!='int' and retval[0]!='float':
        setter_params = [param(class_name,'obj')]+params+[param('const %s&' % retval[0],'data')]
    else:
        setter_params = [param(class_name,'obj')]+params+[param('%s' % retval[0],'data')]
    cl.add_function_as_method('get_%s' % method_name,retval,getter_params,custom_name='get%s' % method_name)
    cl.add_function_as_method('set_%s' % method_name,None,setter_params,custom_name='set%s' % method_name)
    
#==================================================================================

mod = Module('pteros',cpp_namespace='pteros')
mod.add_include('<numpy/noprefix.h>')
mod.add_include('"pteros/core/system.h"')
mod.add_include('"pteros/core/selection.h"')

mod.after_init.write_code('import_array();')

# Pre-register all classes:
Energy_components = mod.add_struct('Energy_components')
Frame = mod.add_class('Frame')
System = mod.add_class('System')
Selection = mod.add_class('Selection')
Periodic_box = mod.add_class('Periodic_box')
Periodic_box = mod.add_class('Atom')

# Pre-register stl containers:
mod.add_container('std::vector<int>', 'int', 'vector')
mod.add_container('std::vector<char>', 'char', 'vector')
mod.add_container('std::vector<std::string>', 'std::string', 'vector')
mod.add_container('std::vector<float>', 'float', 'vector')
mod.add_container('std::vector<pteros::Selection>', 'pteros::Selection', 'vector')

#~~~~~~~~~~~~~~~~~~~~~~~~~~
# Energy_components
#~~~~~~~~~~~~~~~~~~~~~~~~~~
Energy_components.add_constructor([])
Energy_components.add_instance_attribute('total','float')
Energy_components.add_instance_attribute('lj_14','float')
Energy_components.add_instance_attribute('q_14','float')
Energy_components.add_instance_attribute('lj_sr','float')
Energy_components.add_instance_attribute('q_sr','float')
Energy_components.add_method('to_str',retval('std::string'),[])

#~~~~~~~~~~~~~~~~~~~~~~~~~~
# Frame
#~~~~~~~~~~~~~~~~~~~~~~~~~~
Frame.add_constructor([])

#~~~~~~~~~~~~~~~~~~~~~~~~~~
# System
#~~~~~~~~~~~~~~~~~~~~~~~~~~
System.add_constructor([])
System.add_constructor([param('std::string','fname')])
System.add_constructor([param('const System&','other')])
System.add_method('append',None,[param('const System&','sys')])
System.add_method('append',None,[param('const Selection&','sel')])
System.add_method('num_frames', retval('int'),[])
System.add_method('num_atoms', retval('int'),[])
System.add_method('load', None ,[param('std::string','fname'),                                 
                                 param('int','b',default_value='0'),
                                 param('int','e',default_value='-1'),
                                 param('int','skip',default_value='0'),
                                 param('load_callback','on_frame',default_value='0')
                                ])
System.add_method('frame_dup',retval('int'),[param('int','fr')])
System.add_method('frame_append',None,[param('const Frame&','fr')])
System.add_method('frame_copy',None,[param('int','fr1'),param('int','fr2')])
System.add_method('frame_delete',None,[param('int','b',default_value='0'),
                                       param('int','e',default_value='-1')])
lvalue_return(mod,System,'Box',retval('pteros::Periodic_box'),[param('int','fr')])
lvalue_return(mod,System,'Time',retval('float'),[param('int','fr')])
lvalue_return(mod,System,'XYZ',retval('Eigen::Vector3f'),[param('int','ind'),param('int','fr')])
lvalue_return(mod,System,'Atom_data',retval('pteros::Atom'),[param('int','fr')])
lvalue_return(mod,System,'Frame_data',retval('pteros::Frame'),[param('int','fr')])
System.add_method('dssp',None,[param('std::string','fname')])
System.add_method('dssp','std::string',[])
System.add_method('atoms_dup',None,[param('const std::vector<int>&','ind'),
                                    param('Selection*','res_sel',default_value='NULL',
                                                                 transfer_ownership=False)])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Selection
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Selection.add_constructor([])
Selection.add_constructor([param('const System&','sys'),param('std::string','str')])
Selection.add_constructor([param('const System&','sys')])
Selection.add_constructor([param('const System&','sys'),param('int','ind1'),param('int','ind2')])
Selection.add_constructor([param('const Selection&','sel')])
Selection.add_method('size',retval('int'),[])
Selection.add_method('append',None,[param('const Selection&','sel')])
Selection.add_method('append',None,[param('int','ind')])


#==========================================
mod.generate(FileCodeSink(open('bindings.cpp','w')))

