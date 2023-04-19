# -*- coding: utf-8 -*-
"""
Created on April 12 2023

Create a moving average plot of a csv data file for x amount of timesteps

@author: steven.hosper
"""

# Importing Libraries
import configuration as config
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Import time-series data
discharge = pd.read_csv(config.path + f'/output/{config.scenario}/maximumDischarge.csv',
                        sep=';',
                        names=["Time", "Discharge"])

discharge.head()

# Creating exponential running mean for x amount of seconds (now 60)
x = 3600
discharge[f'EWMA{x}'] = discharge['Discharge'].ewm(span=x).mean()

discharge.head()


# Plotting
plt.style.use('default')
discharge[['Discharge', f'EWMA{x}']].plot(label='Discharge',
                                          figsize=(16,8))
plt.axvline(900, color = 'r')
plt.axvline(2700, color = 'r')
#plt.axvline(2700, color = 'g')
#plt.axvline(3600, color = 'g')


plt.show()