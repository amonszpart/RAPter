#!/usr/bin/python

import os
import sys

if len(sys.argv) != 3:
    print 'Usage: divide.py a b, returns a/b'

print int(float(sys.argv[1]) / float(sys.argv[2]))