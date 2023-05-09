#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import configparser

# I used the PCR-GLOBWB methodology as guide to learn how to use the code
# Should I refer to this, or tell of it in any way??

class Configuration(object):
    
    def __init__(self, ini_file):
        
        object.__init__(self)
        
        # Check if ini-file is given
        if ini_file is None:
            raise Exception('Error: No configuration file found')
        
        self.work_dir    = os.getcwd()
        
        self.ini_file_path = os.path.join(self.work_dir, ini_file)
        
        self.parseConfigFile(self.ini_file_path)
        
        
        
    def parseConfigFile(self, file_name):
        
        # Activate parser and read ini-file
        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(file_name)
        
        # Read the different available setting groups
        self.groups = config.sections()
        
        # Read all settings for every group
        for group in self.groups:
            vars(self)[group] = {}              # to instantiate self.GENERAL
            options = config.options(group)
            for opt in options:
                val = config.get(group, opt)
                self.__getattribute__(group)[opt] = val

# For in-file testing
if __name__ == '__main__':
    configuration = Configuration("F:/Projecten intern (2023)/Stage Steven Hosper/Model/v1/config.ini")
    