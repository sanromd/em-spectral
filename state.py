#!/usr/bin/env python
# encoding: utf-8
#
class State(object):

    _initialized = 1 # sets the initialized flag
    t = 0 #initializes at time t = 0 by default
    parameters  = None
    q           = None
    dimensions  = None
    grid        = None
    solver      = None
    aux         = None
    _q_old      = None
    def __init__(self, **kargs):
        cls_dict = self.__class__.__dict__
        for key in kargs:
            assert key in cls_dict
        self.__dict__.update(kargs)

    def etar(self,t=0):
        pass