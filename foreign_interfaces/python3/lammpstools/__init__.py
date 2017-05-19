__all__ = [ 'block_data.py', 'test_me.py', 'dump_reader.py' ]

import sys, os
sys.path.append( os.path.dirname(__file__) )

print("Path in __init__.py is: ", sys.path)

import block_data
import dump_reader
import test
