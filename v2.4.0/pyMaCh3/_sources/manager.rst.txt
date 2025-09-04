Manager
=======

This module handles the high level stuff like config options and YAML stuff. The main class is the :py:class:`pyMaCh3._pyMaCh3.manager.Manager` class. \
You can read more about the manager and config files on `the wiki page <https://github.com/mach3-software/MaCh3/wiki/01.-Manager-and-config-handling>`_. \
YAML stuff works essentially the same as in the c++ version but with some caveats. 
The main difference is that in the python version, the way that you access the actual data of a yaml node is different due to the way the python binding of the c++ code works. 

Lets take the following YAML snippet as an example ::
   
   Example: 
      Int: 1234
      Str: 'TestString'
      IntArray: ['0', '1', '2'] 
      StrArray: ['val1', 'val2', 'val3'] 

In the c++ version you would access these by doing ::

   // for the integer 
   int testInt = node['Int'].as<int>() 
   // for the string 
   std::string testStr = node['Str'].as<std::string>() 
   // for the integer array 
   ... = node['IntArray'].as<std::vector<int>>() 
   // and for the string array 
   ... = node['StrArray'].as<std::vector<std::string>>() 

In the python version however things are different. Arrays are interpreted as arrays of YAML nodes and all of the underlying data are stored as strings. \
So to access the data as above you would need to do ::

   test_int = int(node['Int'].data()) 
   test_str = node['Str'].data() 
   ## then for the arrays 
   int1 = int(node['IntArray'][0].data()) 
   int2 = int(node['IntArray'][1].data()) 
   int3 = int(node['IntArray'][2].data()) 

   str1 = node['StrArray'][0].data() 
   str2 = node['StrArray'][1].data() 
   str3 = node['StrArray'][2].data() 

Parsing arrays can be made a bit less painful using list comprehension ::

   int_list = [int(i.data()) for i in node['IntArray']] 
   str_list = [i.data() for i in node['StrArray']] 


.. automodapi:: pyMaCh3._pyMaCh3.manager
   :members:
   :undoc-members:
   :show-inheritance:
