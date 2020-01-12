# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import pandas as pd

fire = pd.read_csv('LFB Incident data from January 2017.csv', encoding = 'latin1')
dwelling2018 = fire.loc[(fire['CalYear'] == 2018) & (fire['IncidentGroup'] != 'False Alarm') & (fire['PropertyCategory'] == 'Dwelling')]
vehicle2018 = fire.loc[(fire['CalYear'] == 2018) & (fire['IncidentGroup'] != 'False Alarm') & (fire['PropertyCategory'] == 'Road Vehicle')]

dwelling2018.to_csv('2018 dwelling fire.csv')
