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
        
        # Check if ini-file is given
        if iniFile is None:
            raise Exception('Error: No configuration file found')
        
        self.workDir    = os.getcwd()
        
        self.iniFilePath = os.path.join(self.workDir, iniFile)
        
        
        
    def parseConfigFile(self, fileName):
        
        # Activate parser and read ini-file
        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(fileName)
        
        # Read the different available setting groups
        self.groups = config.sections()
        
        # Read all settings for every group
        for group in self.groups:
            vars(self)[group] = {}              # to instantiate self.GENERAL
            options = config.options(group)
            for opt in options:
                val = config.get(group, opt)
                self.__getattribute__(group)[opt] = val