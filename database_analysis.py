import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import json
from datetime import datetime

import Bio
import Bio.PDB
import Bio.SeqRecord

cwd = os.getcwd()
#datadir = cwd + '/pdb'
datadir = '/vault/privateer_database/pdb'
def file_paths(root_directory):
    filepathlist = []
    for root, dirs, files in os.walk(root_directory):
        for f in files:
            filepathlist.append(os.path.join(root,f))
    return filepathlist

def get_year_pdb(jsonfilepath):
    jsonfilename = os.path.basename(jsonfilepath)
    pdbcode = jsonfilename.rpartition('.')[0]
    print('Looking for pdb file corresponding to ' + jsonfilepath)
    pdbfilepath = '/vault/pdb/pdb' + pdbcode + '.ent'
    if os.path.isfile(pdbfilepath):
        pdbheader = Bio.PDB.parse_pdb_header(pdbfilepath)
        datestring = pdbheader['deposition_date']
        dt = datetime.strptime(datestring, '%Y-%m-%d')
        return int(dt.year)
    else:
        return 123456789

def get_res_pdb(jsonfilepath):
    jsonfilename = os.path.basename(jsonfilepath)
    pdbcode = jsonfilename.rpartition('.')[0]
    pdbfilepath = '/vault/pdb/pdb' + pdbcode + '.ent'
    pdbheader = Bio.PDB.parse_pdb_header(pdbfilepath)
    res = pdbheader['resolution']
    return res

def glycans_per_year(databasedir):
    print('Analysing database at ' + databasedir)
    filepathlist = file_paths(databasedir)
    glycans = np.zeros(len(filepathlist))
    nglycans = np.zeros(len(filepathlist))
    oglycans = np.zeros(len(filepathlist))
    years = np.zeros(len(filepathlist))
    for i in range(len(filepathlist)):
        jsonfile = filepathlist[i]
        years[i] = get_year_pdb(jsonfile)
        with open(jsonfile, 'r') as f:
            data = json.load(f)
            glycandata = data['glycans']
            nglycan = glycandata['n-glycan']
            oglycan = glycandata['o-glycan']
            sglycan = glycandata['s-glycan']
            cglycan = glycandata['c-glycan']
            ligand = glycandata['ligand']
            if len(nglycan)>0:
                nglycans[i] = 1
            if len(oglycan)>0:
                oglycans[i] = 1
            if len(nglycan)>0 or len(oglycan)>0 or len(sglycan)>0 or len(cglycan)>0 or len(ligand)>0:
                glycans[i] = 1
    years = np.delete(years, np.where(years == 123456789)[0])

    start = np.min(years)
    end = np.max(years)
    year_range = np.arange(start,end)
    nglycansperyear = np.zeros(len(year_range))
    oglycansperyear = np.zeros(len(year_range))
    glycansperyear = np.zeros(len(year_range))
    for i in range(len(year_range)):
        for j in range(len(years)):
            if years[j] == year_range[i]:
                nglycansperyear[i] += nglycans[j]
                oglycansperyear[i] += oglycans[j]
                glycansperyear[i] += glycans[j]
    plt.plot(year_range, glycansperyear, label='Carbohydrate entries in PDB')
    plt.plot(year_range, nglycansperyear, label='N-glycosylated')
    plt.plot(year_range, oglycansperyear, label='O-glycosylated')
    plt.xlabel('Release Year')
    plt.ylabel('Carbohydrate Entries in PDB')
    plt.legend()
    plt.savefig(cwd+'/glycosylated_per_year.png')
    plt.close()                
    return

glycans_per_year(datadir)
