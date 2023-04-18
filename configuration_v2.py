#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import configparser

# I used the PCR-GLOBWB methodology as guide to learn how to use the code
# Should I refer to this, or tell of it in any way??

class Configuration(object):
    
    def __init__(self, iniFile):
        
        object.__init__(self)
        
        if iniFile is None:
            raise Exception('Error: No configuration file found')
        
        self.workDir    = os.getcwd()
        
        self.iniFilePath = os.path.join(self.workDir, iniFile)
        
        
        
    def parseConfigFile(self, fileName):
        
        config = configparser.ConfigParser()