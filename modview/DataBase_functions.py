#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sqlite3
conn = sqlite3.connect('LoaderDB3.db') # connection
conn.row_factory = sqlite3.Row
c = conn.cursor()


# In[2]:


def Prepare_id(**kwargs):
    # This function runs through items in **kwargs and uses them to build an sql statement to execute.
    # At the end, it returns all database rows whose column names and entries match with **kwargs 
    # Example usage so far has only been **kwargs = {t0='2018 Sep 15'}. Namely, return all observational objects
    # with the specific t0.
    Data = "SELECT * FROM OBJ_ID"
    i = 0
    for key, value in kwargs.items():
        if i == 0:
            Data += " WHERE "
        else:
            Data += " AND "
        Data += "{}='{}'".format(key, value)
        i += 1
    Data += ";" 
    return c.execute(Data)


def format_change(DICT):
    # This function turns 
    limits_keys = ['t0','t1','lon','lat']
    for n in range(len(DICT)):
        DICT[n]['limits'] = {key: DICT[n][key] for key in DICT[n] if key in limits_keys}
        for key in limits_keys: del DICT[n][key]
            
def id_maker(**kwargs):
    Prepare_id(**kwargs)
    A = [dict(row) for row in c.fetchall()]
    format_change(A)
    return A

def Data_extraction(obj_id, n, **kwargs):
    ''' obj_id is a list of dictionaries (returned by Prepare_id) each of which includes
    all fields that identify observational objects matching some request. 
    n is the index of the item of interest within obj_id. '''
    
    combiner_id = obj_id[n]['COMBINER_ID'] # use combiner table to pull datasources that make up obj_id[n]
    common_str = 'SELECT * FROM COMBINER AS c INNER JOIN DATASOURCES as d ON c.Data_obj = d.Data_obj WHERE c.COMBINER_ID = %d'
    i = 1
    for key, value in kwargs.items(): # examples so far have not used **kwargs at all. 
        common_str += " AND d." # add new clauses to database statement to execute
        common_str += "{}='{}'".format(key, value)
        i += 1
    common_str += ";"
    c.execute(common_str %combiner_id) # execute statement (returns all entries in datasources table that corresponds to obj_id[0]
    return [dict(row) for row in c.fetchall()] 

def path_maker(DataSources, data_dir):
    # Take output from data_extraction and return a dictionary that can be run through loader.assemble
    Dir = {}
    for kk in range(len(DataSources)):
        if DataSources[kk]['Source'] not in Dir:
            if DataSources[kk]['Source'] != 'ADCP':
                Dir[DataSources[kk]['Source']] = [data_dir + DataSources[kk]['Filepath']]
            else:
                Dir[DataSources[kk]['Source']] = [data_dir + 'ADCP/' + DataSources[kk]['Filepath']]
        elif DataSources[kk]['Source'] in Dir:
            if DataSources[kk]['Source'] != 'ADCP':
                Dir[DataSources[kk]['Source']].append(data_dir + DataSources[kk]['Filepath'])
            else:
                Dir[DataSources[kk]['Source']].append(data_dir + 'ADCP/' + DataSources[kk]['Filepath'])
    return Dir

