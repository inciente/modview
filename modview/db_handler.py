import sqlite3
''' Functions to search pathfiles in datasets.db and format them to
act as input for loader.assemble objects '''

conn = sqlite3.connect('datasets.db') # connect to database
conn.row_factory = sqlite3.Row
c = conn.cursor();

def Prepare_id(**kwargs):
    ''' Run through items in **kwargs and translate to sql statement.
    Return all rows in table OBJ_ID whose column names and entries 
    match with **kwargs '''
    statement = "SELECT * FROM OBJ_ID"
    i = 0
    for key, value in kwargs.items():
        if i == 0:
            statement += " WHERE "
        else:
            statement += " AND "
        statement += "{}='{}'".format(key, value)
        i += 1
    statement += ";"
    return c.execute(statement)

# Comment from Noel: is this function finished?
def format_change(DICT):
    # This function turns
    limits_keys = ['t0','t1','lon','lat']
    for n in range(len(DICT)):
        DICT[n]['limits'] = {key: DICT[n][key] for key in DICT[n] \
                if key in limits_keys}
        for key in limits_keys: del DICT[n][key]

def id_maker(**kwargs):
    ''' Execute search statement with **kwargs and change format of output
    into format understood by loader.assemble '''
    Prepare_id(**kwargs)
    A = [dict(row) for row in c.fetchall()]
    format_change(A)
    return A

''' seems like Data_extraction returns dictionaries with the filepath for 
datasources within observational object obj_id[n].'''
def Data_extraction(obj_id, n, **kwargs):
    ''' obj_id is a list of dictionaries (returned by Prepare_id) 
    each of which includes all fields that identify observational 
    objects matching some search statement. 
    n is the index of the item of interest within obj_id.'''
    # use combiner table to pull datasources that make up obj_id[n]
    combiner_id = obj_id[n]['COMBINER_ID'] 
    common_str = 'SELECT * FROM COMBINER AS c INNER JOIN DATASOURCES' \
            ' as d ON c.Data_obj = d.Data_obj WHERE c.COMBINER_ID = %d'
    i = 1
    for key, value in kwargs.items(): # examples so far have not used **kwargs 
        common_str += " AND d." # add new clauses to database statement to execute
        common_str += "{}='{}'".format(key, value)
        i += 1
    common_str += ";"
    c.execute(common_str %combiner_id) # returns datasources linked to obj_id[0]
    return [dict(row) for row in c.fetchall()] 

# Comment from Noel: is the idea here that every separate assemble request must
# come with the directory data_dir of corresponding data? 
# Directories need to be within datasources table, but I see that they're not. 
# Many of the if else statements in this function will be eliminated once that's done.
def path_maker(DataSources, data_dir):
    # Take output from data_extraction and return a dictionary that can be run through loader.assemble
    Dir = {}
    for kk in range(len(DataSources)):
        if DataSources[kk]['Source'] not in Dir:
            if DataSources[kk]['Source'] != 'ADCP':
                Dir[DataSources[kk]['Source']] = \
                        [data_dir + DataSources[kk]['Filepath']]
            else:
                Dir[DataSources[kk]['Source']] = \
                        [data_dir + 'ADCP/' + DataSources[kk]['Filepath']]
        elif DataSources[kk]['Source'] in Dir:
            if DataSources[kk]['Source'] != 'ADCP':
                Dir[DataSources[kk]['Source']].append( \
                        data_dir + DataSources[kk]['Filepath'])
            else:
                Dir[DataSources[kk]['Source']].append(data_dir + 'ADCP/' + DataSources[kk]['Filepath'])
    return Dir

