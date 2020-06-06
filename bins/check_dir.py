#!/usr/bin/python3
import os
import sys

for i in sys.argv:
    if not os.path.exists(str(i)):
        os.makedirs(str(i))
