#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sqlite3
conn = sqlite3.connect('LoaderDB3.db') # connection
conn.row_factory = sqlite3.Row
c = conn.cursor()


# In[2]:


def Prepare_id(**kwargs):
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
    combiner_id = obj_id[n]['COMBINER_ID']
    common_str = 'SELECT * FROM COMBINER AS c INNER JOIN DATASOURCES as d ON c.Data_obj = d.Data_obj WHERE c.COMBINER_ID = %d'
    i = 1
    for key, value in kwargs.items():
        common_str += " AND d."
        common_str += "{}='{}'".format(key, value)
        i += 1
    common_str += ";"
    c.execute(common_str %combiner_id)
    return [dict(row) for row in c.fetchall()] 

def path_maker(DataSources, data_dir):
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

