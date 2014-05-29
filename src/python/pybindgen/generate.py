#!/usr/bin/python

from pybindgen import *
from pybindgen.typehandlers.base import ForwardWrapperBase

import sys

IN = Parameter.DIRECTION_OUT
OUT = Parameter.DIRECTION_OUT

class EigenParam(Parameter):

    DIRECTIONS = []
    CTYPES = []
    PYARRAY_DTYPE = None
    CPP_DTYPE = None
    CAST_DTYPE = None
    NDIM = 1
    DIM1 = -1
    DIM2 = -1
   
    #def __init__(self, ctype, name, direction=1):
    #    super(EigenParam,self).__init__(ctype, name, direction)
      
    def convert_python_to_c(self, wrapper):
        assert isinstance(wrapper, ForwardWrapperBase)
        
        dim1=''
        dim2=''
        
        if self.direction & Parameter.DIRECTION_IN:
            # Declare pyobjects            
            arg = wrapper.declarations.declare_variable("PyObject*", self.name+'_arg','NULL')
            arr = wrapper.declarations.declare_variable("PyObject*", self.name+'_arr','NULL')

            # Parse pyobject from python params
            if self.default_value==None:
                wrapper.parse_params.add_parameter('O', ['&'+arg], self.value)
                # Now we need to create an array from obtained object        
                wrapper.before_call.write_code('%s = PyArray_FROM_OTF(%s, %s, NPY_IN_ARRAY);' % 
                                                (arr,arg,self.PYARRAY_DTYPE))
            else:                            
                wrapper.parse_params.add_parameter('O', ['&'+arg], self.value, optional=True)
                wrapper.before_call.write_code('if(!%s){' % arg)
                wrapper.before_call.indent()
                # Just create a pyarray with default value
                if self.NDIM==1:
                    wrapper.before_call.write_code('npy_intp %s_sz[1]={%s};' 
                                                            % (self.name,self.DIM1))
                    wrapper.before_call.write_code('%s = PyArray_SimpleNew(1, %s_sz, %s);' 
                            % (arr,self.name,self.PYARRAY_DTYPE))
                    wrapper.before_call.write_code('Eigen::Map<%s> %s_def((%s*)PyArray_DATA(%s),%s,%s);' 
                            % (self.CPP_DTYPE,self.name,self.CAST_DTYPE,arr,self.PYARRAY_DTYPE,self.DIM1))
                    wrapper.before_call.write_code('%s_def = %s;' % (self.name,self.default_value))
                else: # 2D case
                    wrapper.before_call.write_code('npy_intp %s_sz[2]={%s,%s};' 
                                                            % (self.name,self.DIM1,self.DIM2))
                    wrapper.before_call.write_code('%s = PyArray_SimpleNew(2, %s_sz, %s);' 
                            % (arr,self.name,self.PYARRAY_DTYPE))
                    wrapper.before_call.write_code('Eigen::Map<%s> %s_def((%s*)PyArray_DATA(%s),%s,%s,%s);' 
                            % (self.CPP_DTYPE,self.name,self.CAST_DTYPE,arr,self.PYARRAY_DTYPE,
                               self.DIM2,self.DIM1))
                    wrapper.before_call.write_code('%s_def = %s;' % (self.name,self.default_value))
                        
                wrapper.before_call.unindent()
                wrapper.before_call.write_code('} else {')
                wrapper.before_call.indent()
                wrapper.before_call.write_code('%s = PyArray_FROM_OTF(%s, %s, NPY_IN_ARRAY);' % 
                                                (arr,arg,self.PYARRAY_DTYPE))
                wrapper.before_call.unindent()
                wrapper.before_call.write_code('}')

            # Add cleanup for this temporary array
            wrapper.before_call.add_cleanup_code("Py_XDECREF(%s);" % arr)
            # Check if conversion is OK
            wrapper.before_call.write_error_check('%s==NULL' % arr, 
                            'PyErr_SetString(PyExc_TypeError,"Convertion to PyArray failed!");')
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
                
            # Do call
            # For Affine3f we need special treatment, so the call string will be passed
            if self.ctype == 'Eigen::Affine3f&':
                wrapper.before_call.write_code('Eigen::Affine3f %s_call(%s_map);' % (self.name,self.name))
                wrapper.call_params.append(self.name+'_call')
            else:
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
                wrapper.declarations.declare_variable('PyObject*','%s_to_map' % self.name, 
                    'PyArray_SimpleNew(1, %s_sz, %s)' % (self.name,self.PYARRAY_DTYPE))
                wrapper.declarations.declare_variable('Eigen::Map<%s>' % self.CPP_DTYPE,
                    '%s_call((%s*)PyArray_DATA(%s_to_map),%s,1)' % (self.name,self.CAST_DTYPE,self.name,dim1))
            else: # 2D case
                wrapper.declarations.declare_variable('npy_intp','%s_sz[2]={%s,%s}' % (self.name,dim1,dim2))
                wrapper.declarations.declare_variable('PyObject*','%s_to_map' % self.name, 
                    'PyArray_SimpleNew(2, %s_sz, %s)' % (self.name,self.PYARRAY_DTYPE))
                wrapper.declarations.declare_variable('Eigen::Map<%s>' % self.CPP_DTYPE,
                    '%s_call((%s*)PyArray_DATA(%s_to_map),%s,%s)' %
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

class Affine3fParam(EigenParam):
    DIRECTIONS = [Parameter.DIRECTION_IN, Parameter.DIRECTION_OUT]
    CTYPES = ['Eigen::Affine3f&']
    PYARRAY_DTYPE = 'NPY_FLOAT'
    CPP_DTYPE = 'Eigen::Matrix4f'
    CAST_DTYPE = 'float'
    NDIM = 2
    DIM1 = 4
    DIM2 = 4   
    # There is a special treatment for Affine3f in base class!

#----------------------------------------------------------------
# Support for std vectors
class std_vector_wrapper(Parameter):
    DIRECTIONS = []
    CTYPES = []    
    
    def pylist_element_to_cpp(self, wrapper):
        # MUST define cpp_el in the code!
        pass

    def cppvector_element_to_py(self, wrapper):
        # MUST define py_el in the code!
        pass 
        
    def extra_cleanup(self,wrapper):
        pass           
    
    def convert_python_to_c(self, wrapper):        
        assert isinstance(wrapper, ForwardWrapperBase)
        # declare a pyobject which holds a sequence
        arg = wrapper.declarations.declare_variable("PyObject*", self.name+'_arg','NULL')
        # declare a c++ vector
        cpp_vec = wrapper.declarations.declare_variable('std::vector<%s>' % self.CPP_DTYPE, self.name+'_cpp_vec')
        
        if self.direction & Parameter.DIRECTION_IN:
            # Parse pyobject from python params
            wrapper.parse_params.add_parameter('O!', ['&PyList_Type','&'+arg], self.value)
            # var for sequence length
            sz = wrapper.declarations.declare_variable('int',self.name+'_sz')            
            # Get size
            wrapper.before_call.write_code('%s = PyList_Size(%s);' % (sz,arg))
            # Cycle over elements
            wrapper.before_call.write_code('for(Py_ssize_t i=0;i<%s;++i){' % sz)
            wrapper.before_call.indent()
            # Get i-th element
            wrapper.before_call.write_code('PyObject* el = PyList_GetItem(%s,i);' % arg)
            
            self.pylist_element_to_cpp(wrapper) # This could be overloaded for different types
            
            # Copy this element to c++ vector
            wrapper.before_call.write_code('%s.push_back(cpp_el);' % cpp_vec)
            # decref i-th element
            wrapper.before_call.write_code('Py_DECREF(el);')
            
            self.extra_cleanup(wrapper)

            wrapper.before_call.unindent()
            wrapper.before_call.write_code('}')
            # Now add call param
            wrapper.call_params.append(cpp_vec)

        if self.direction & Parameter.DIRECTION_OUT:
            # Add call param
            wrapper.call_params.append(cpp_vec)
            # After call cycle over the vector and create a sequence from it            
            wrapper.after_call.write_code('PyObject* outlist_%s = PyList_New(%s.size());' % (self.name,cpp_vec))
            wrapper.after_call.write_code('for(size_t i=0;i<%s.size();++i){' % cpp_vec)
            wrapper.after_call.indent()
            
            self.cppvector_element_to_py(wrapper)
            
            wrapper.after_call.write_code('PyList_SET_ITEM(outlist_%s,i,py_el);' % (self.name))
            wrapper.after_call.unindent()
            wrapper.after_call.write_code('}')
            # And return the param in output tuple
            wrapper.build_params.add_parameter("N", ['outlist_'+self.name])  


class vector_of_EigenVectors(std_vector_wrapper):    
    PYARRAY_DTYPE = ''
    CPP_DTYPE = ''
    CAST_DTYPE = ''
    DIM = -1

    def pylist_element_to_cpp(self, wrapper):
        # MUST define cpp_el in the code!
        # Create a PyArray from i-th element
        wrapper.before_call.write_code('PyObject* arr = PyArray_FROM_OTF(el, %s, NPY_IN_ARRAY);' % 
                                                (self.PYARRAY_DTYPE))
        # Check if convertion Ok
        wrapper.before_call.write_code('if(arr==NULL) PyErr_SetString(PyExc_TypeError,"Convertion to PyArray failed!");')
        # Map i-th element to Eigen
        wrapper.before_call.write_code('Eigen::Map<%s> cpp_el((%s*)PyArray_DATA(arr),%s,1);' %
                                        (self.CPP_DTYPE,self.CAST_DTYPE,self.DIM))
    def extra_cleanup(self,wrapper):
        wrapper.before_call.write_code('Py_XDECREF(arr);')
    
    def cppvector_element_to_py(self, wrapper):
        # MUST define py_el in the code!
        wrapper.after_call.write_code('npy_intp arr_sz[1] = {3};')
        wrapper.after_call.write_code('PyObject* py_el = PyArray_SimpleNew(1, arr_sz, %s);' % self.PYARRAY_DTYPE)
        wrapper.after_call.write_code('Eigen::Map<%s> mp((%s*)PyArray_DATA(arr),3,1);' % 
                                            (self.CPP_DTYPE,self.CAST_DTYPE))
        wrapper.after_call.write_code('mp = %s[i];' % cpp_vec)

            
class vector_of_Vector3f(vector_of_EigenVectors):
    DIRECTIONS = [Parameter.DIRECTION_IN, Parameter.DIRECTION_OUT]
    CTYPES = ['std::vector<Eigen::Vector3f>&']
    PYARRAY_DTYPE = 'NPY_FLOAT'
    CPP_DTYPE = 'Eigen::Vector3f'
    CAST_DTYPE = 'float'
    DIM = 3
    
class vector_of_Vector2i(vector_of_EigenVectors):
    DIRECTIONS = [Parameter.DIRECTION_IN, Parameter.DIRECTION_OUT]
    CTYPES = ['std::vector<Eigen::Vector2i>&']
    PYARRAY_DTYPE = 'NPY_INT'
    CPP_DTYPE = 'Eigen::Vector2i'
    CAST_DTYPE = 'int'
    DIM = 2

class vector_of_int(std_vector_wrapper):    
    DIRECTIONS = [Parameter.DIRECTION_IN, Parameter.DIRECTION_OUT]
    CTYPES = ['std::vector<int>&']
    CPP_DTYPE = 'int'
    
    def pylist_element_to_cpp(self, wrapper):
        # MUST define cpp_el in the code!
        # Create an int from i-th element
        wrapper.before_call.write_code('int cpp_el = PyInt_AsLong(el);')
        # Check if convertion Ok
        wrapper.before_call.write_error_check('PyErr_Occurred()',
            'PyErr_SetString(PyExc_TypeError,"Convertion to int failed!");')        
    
    def cppvector_element_to_py(self, wrapper):
        # MUST define py_el in the code!
        wrapper.after_call.write_code('PyObject* py_el = PyInt_FromLong(%s[i]);' % cpp_vec)


class vector_of_float(std_vector_wrapper):
    DIRECTIONS = [Parameter.DIRECTION_IN, Parameter.DIRECTION_OUT]
    CTYPES = ['std::vector<float>&']
    CPP_DTYPE = 'float'
    
    def pylist_element_to_cpp(self, wrapper):
        # MUST define cpp_el in the code!
        # Create an int from i-th element
        wrapper.before_call.write_code('float cpp_el = (float)PyFloat_AsDouble(el);')
        # Check if convertion Ok
        wrapper.before_call.write_error_check('PyErr_Occurred()',
            'PyErr_SetString(PyExc_TypeError,"Convertion to float failed!");')        
    
    def cppvector_element_to_py(self, wrapper):
        # MUST define py_el in the code!
        wrapper.after_call.write_code('PyObject* py_el = PyFloat_FromDouble((double)%s[i]);' % cpp_vec)


class vector_of_string(std_vector_wrapper):
    DIRECTIONS = [Parameter.DIRECTION_IN, Parameter.DIRECTION_OUT]
    CTYPES = ['std::vector<std::string>&']
    CPP_DTYPE = 'std::string'
    
    def pylist_element_to_cpp(self, wrapper):
        # MUST define cpp_el in the code!
        # Create an int from i-th element
        wrapper.before_call.write_code('char* buf = PyString_AsString(el);')
        # Check if convertion Ok
        wrapper.before_call.write_error_check('buf==NULL',
            'PyErr_SetString(PyExc_TypeError,"Convertion to string failed!");')        
        wrapper.before_call.write_code('std::string cpp_el(buf);')
        
    def cppvector_element_to_py(self, wrapper):
        # MUST define py_el in the code!
        wrapper.after_call.write_code('PyObject* py_el = PyString_FromString(%s[i].c_str());' % cpp_vec)

class vector_of_char(std_vector_wrapper):
    DIRECTIONS = [Parameter.DIRECTION_IN, Parameter.DIRECTION_OUT]
    CTYPES = ['std::vector<char>&']
    CPP_DTYPE = 'char'
    
    def pylist_element_to_cpp(self, wrapper):
        # MUST define cpp_el in the code!
        # Create an int from i-th element
        wrapper.before_call.write_code('char* buf = PyString_AsString(el);')
        # Check if convertion Ok
        wrapper.before_call.write_error_check('buf==NULL',
            'PyErr_SetString(PyExc_TypeError,"Convertion to char failed!");')
        wrapper.before_call.write_error_check('strlen(buf)!=1',
            'PyErr_SetString(PyExc_TypeError,"Single character expected!");')        
        wrapper.before_call.write_code('char cpp_el(buf[0]);')
        
    def cppvector_element_to_py(self, wrapper):
        # MUST define py_el in the code!
        wrapper.after_call.write_code('PyObject* py_el = PyString_FromFormat("\%c",%s[i]);' % cpp_vec)
    
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
        wrapper.declarations.declare_variable('Eigen::Map<Eigen::Vector3f>',
            'retval((float*)PyArray_DATA(retval_to_map),3,1)')
        wrapper.build_params.add_parameter("N", ['retval_to_map'],prepend=True)


class Affine3fReturn(ReturnValue):
    CTYPES = ['Eigen::Affine3f']
    
    def __init__(self, ctype):
        saved = ctype.ctype
        super(Affine3fReturn, self).__init__("#define _a_hack_to_disable_retval_generation")
        self.ctype = saved

    def convert_c_to_python(self, wrapper):
        wrapper.declarations.declare_variable('npy_intp','retval_sz[2] = {4,4}')
        wrapper.declarations.declare_variable('PyObject*','retval_to_map = PyArray_SimpleNew(2, retval_sz, NPY_FLOAT)')
        wrapper.declarations.declare_variable('Eigen::Affine3f','retval')
        wrapper.declarations.declare_variable('Eigen::Map<Eigen::Matrix4f>',
            'retval_map((float*)PyArray_DATA(retval_to_map),4,4)')
        wrapper.after_call.write_code('retval_map = retval.matrix();')
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
        wrapper.declarations.declare_variable('Eigen::Map<Eigen::Matrix3f>',
            'retval((float*)PyArray_DATA(retval_to_map),3,3)')
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
        assert self.dim1_code!=None and self.dim2_code!=None
        wrapper.declarations.declare_variable('npy_intp','retval_sz[2] = {%s,%s}' %
                                                 (self.dim2_code,self.dim1_code))
        wrapper.declarations.declare_variable('PyObject*','retval_to_map = PyArray_SimpleNew(2, retval_sz, NPY_FLOAT)')
        wrapper.declarations.declare_variable('Eigen::Map<Eigen::MatrixXf>','retval((float*)PyArray_DATA(retval_to_map),retval_sz[1],retval_sz[0])')
        wrapper.build_params.add_parameter("N", ['retval_to_map'],prepend=True)


class std_vector_return(ReturnValue):
    CTYPES = []
    def cppvector_element_to_py(self,wrapper):
        pass
        
    def convert_c_to_python(self, wrapper):
        # After call cycle over the vector and create a sequence from it            
        wrapper.after_call.write_code('PyObject* retlist = PyList_New(retval.size());')
        wrapper.after_call.write_code('for(size_t i=0;i<retval.size();++i){')
        wrapper.after_call.indent()
        
        self.cppvector_element_to_py(wrapper)
        
        wrapper.after_call.write_code('PyList_SET_ITEM(retlist,i,py_el);')
        wrapper.after_call.unindent()
        wrapper.after_call.write_code('}')
        # And return the param in output tuple
        wrapper.build_params.add_parameter("N", ['retlist'],prepend=True)

class std_vector_int_return(std_vector_return):
    CTYPES = ['std::vector<int>']
    def cppvector_element_to_py(self,wrapper):
        wrapper.after_call.write_code('PyObject* py_el = PyInt_FromLong(retval[i]);')

class std_vector_float_return(std_vector_return):
    CTYPES = ['std::vector<float>']
    def cppvector_element_to_py(self,wrapper):
        wrapper.after_call.write_code('PyObject* py_el = PyFloat_FromDouble((double)retval[i]);')

class std_vector_string_return(std_vector_return):
    CTYPES = ['std::vector<std::string>']
    def cppvector_element_to_py(self,wrapper):
        wrapper.after_call.write_code('PyObject* py_el = PyString_FromString(retval[i].c_str());')

class std_vector_char_return(std_vector_return):
    CTYPES = ['std::vector<char>']
    def cppvector_element_to_py(self,wrapper):
        wrapper.after_call.write_code('PyObject* py_el = PyString_FromStringAndSize(&retval[i],1);')


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
mod.header.writeln('using namespace pteros;')
mod.header.writeln('using namespace Eigen;')

# Pre-register all classes:
Energy_components = mod.add_struct('Energy_components')
Frame = mod.add_class('Frame')
System = mod.add_class('System')
Selection = mod.add_class('Selection')
Periodic_box = mod.add_class('Periodic_box')
Periodic_box = mod.add_class('Atom')

# Pre-register stl containers:
mod.add_container('std::vector<pteros::Selection>', 'pteros::Selection', 'vector')
mod.add_container('std::vector<pteros::Atom>', 'pteros::Atom', 'vector')

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
System.add_method('atoms_add',None,[param('const std::vector<pteros::Atom>&','atm'),
                                    param('const std::vector<Eigen::Vector3f>&','crd'),
                                    param('Selection*','res_sel',default_value='NULL',
                                                                 transfer_ownership=False)])
System.add_method('atoms_delete',None,[param('const std::vector<int>&','ind')])
System.add_method('distance',retval('float'),[param('int','i'),
                                              param('int','j'),
                                              param('int','fr'),
                                              param('bool','is_periodic',default_value='true'),
                                              param('Vector3i_const_ref','dims',default_value='Eigen::Vector3i::Ones()')
                                             ])
System.add_method('wrap_all',None,[param('int','fr'),
                                   param('Vector3i_const_ref','dims_to_wrap',default_value='Eigen::Vector3i::Ones()')
                                  ])
System.add_method('non_bond_energy',retval('pteros::Energy_components'),
                    [param('const std::vector<Eigen::Vector2i>&','nlist'),
                     param('int','fr'),
                     param('bool','is_periodic',default_value='true')
                    ])
System.add_method('clear',None,[])
System.add_method('force_field_ready',retval('bool'),[])
System.add_method('assign_resindex',None,[])
System.add_method('sort_by_resindex',None,[])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Selection
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Selection.add_constructor([])
Selection.add_constructor([param('const System&','sys'),param('std::string','str')])
Selection.add_constructor([param('const System&','sys')])
Selection.add_constructor([param('const System&','sys'),param('int','ind1'),param('int','ind2')])
Selection.add_constructor([param('const Selection&','sel')])
Selection.add_method('append',None,[param('const Selection&','sel')])
Selection.add_method('append',None,[param('int','ind')])
Selection.add_method('modify',None,[param('const System&','sys'),param('std::string','str')])
Selection.add_method('modify',None,[param('std::string','str')])
Selection.add_method('modify',None,[param('const System&','sys')])
Selection.add_method('modify',None,[param('const System&','sys'),param('int','ind1'),param('int','ind2')])
Selection.add_method('modify',None,[param('const std::vector<int>&','ind')])
Selection.add_method('apply',None,[])
Selection.add_method('update',None,[])
Selection.add_method('clear',None,[])
Selection.add_method('get_frame',retval('int'),[])
Selection.add_method('set_frame',None,[param('int','fr')])
Selection.add_method('get_system',retval('System*',reference_existing_object=True),[])
Selection.add_method('get_text',retval('std::string'),[])
Selection.add_method('get_index',retval('std::vector<int>'),[])
Selection.add_method('get_chain',retval('std::vector<char>'),[])
Selection.add_method('set_chain',None,[param('const std::vector<char>&','data')])
Selection.add_method('set_chain',None,[param('char','data')])
Selection.add_method('get_unique_chain',retval('std::vector<char>'),[])
Selection.add_method('get_resid',retval('std::vector<int>'),[])
Selection.add_method('set_resid',None,[param('const std::vector<int>&','data')])
Selection.add_method('set_resid',None,[param('int','data')])
Selection.add_method('get_unique_resid',retval('std::vector<int>'),[])
Selection.add_method('get_resindex',retval('std::vector<int>'),[])
Selection.add_method('get_unique_resindex',retval('std::vector<int>'),[])
Selection.add_method('get_name',retval('std::vector<std::string>'),[])
Selection.add_method('set_name',None,[param('const std::vector<std::string>&','data')])
Selection.add_method('set_name',None,[param('std::string&','data')])
Selection.add_method('get_resname',retval('std::vector<std::string>'),[])
Selection.add_method('set_resname',None,[param('const std::vector<std::string>&','data')])
Selection.add_method('set_resname',None,[param('std::string&','data')])
Selection.add_method('get_xyz',retval('Eigen::MatrixXf',dim1='3',dim2='self->obj->size()'),[])

Selection.add_method('size',retval('int'),[])
#==========================================
mod.generate(FileCodeSink(open('bindings.cpp','w')))

