import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import json
from datetime import datetime
import requests
from bs4 import BeautifulSoup
import gzip

import Bio
import Bio.PDB
import Bio.SeqRecord

cwd = os.getcwd()
#datadir = cwd + '/pdb'
datadir = '/vault/privateer_database/pdb/'
pdbdir = '/vault/pdb_mirror/data/structures/all/pdb/' # Up to date
# pdbdir = '/vault/pdb/' #Old
def file_paths(root_directory):
    """
    Function to find all the files in the database
    """
    filepathlist = []
    for root, dirs, files in os.walk(root_directory):
        for f in files:
            filepathlist.append(os.path.join(root,f))
    return filepathlist

def get_pdb_path(pdb_code: str) -> str:
    """Get PDB path, if the PDB is zipped, then unzip and return the new path 

    Args:
        pdb_code (str): 4 letter PDB code

    Returns:
        str: path to the unzipped PDB file
    """
    gzipped_pdb = f"/vault/pdb_mirror/data/structures/all/pdb/pdb{pdb_code}.ent.gz"
    unzipped_pdb = f"/vault/tmp_extracted_pdbs/pdb{pdb_code}.ent"

    if os.path.exists(unzipped_pdb):
        return unzipped_pdb

    with gzip.open(gzipped_pdb, 'rb') as f_in:
        with open(unzipped_pdb, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return unzipped_pdb

def get_year_pdb_old_mirror(jsonfilepath, pdbdir):
    """
    Function to find the year the structure corresponding to a particular database entry was deposited
    Version of function to run with up to date mirror /vault/pdb
    Have to check if files are there as missing most recent ones.
    """
    jsonfilename = os.path.basename(jsonfilepath)
    pdbcode = jsonfilename.rpartition('.')[0]
    # print('Looking for pdb file corresponding to ' + jsonfilepath)
    pdbfilepath = pdbdir + 'pdb' + pdbcode + '.ent'
    if os.path.isfile(pdbfilepath): #Checking if file exists
        pdbheader = Bio.PDB.parse_pdb_header(pdbfilepath)
        datestring = pdbheader['deposition_date']
        dt = datetime.strptime(datestring, '%Y-%m-%d')
        return int(dt.year)
    else:
        print('Failed to find corresponding pdb at ' + pdbfilepath)
        return 123456789

def get_year_pdb(jsonfilepath, pdbdir):
    """
    Function to find the year the structure corresponding to a particular database entry was deposited
    Version of function to run with up to date mirror /vault/pdb_mirror/data/structures/all/pdb
    Can no longer check file existence as os doesn't like .gz files.
    But doesn't currently work as Bio.PDB.parse_pdb_header also won't open so using the old mirror for now.
    """
    jsonfilename = os.path.basename(jsonfilepath)
    pdbcode = jsonfilename.rpartition('.')[0]
    # print('Looking for pdb file corresponding to ' + jsonfilepath)
    pdbfilepath = get_pdb_path(pdbcode)
    pdbheader = Bio.PDB.parse_pdb_header(pdbfilepath)
    datestring = pdbheader['deposition_date']
    dt = datetime.strptime(datestring, '%Y-%m-%d')
    return int(dt.year)


def get_res_pdb(jsonfilepath, pdbdir):
    """
    Function to find resolution of structure corresponding to particular database entry
    """
    jsonfilename = os.path.basename(jsonfilepath)
    pdbcode = jsonfilename.rpartition('.')[0]
    print('Looking for pdb file corresponding to ' + jsonfilepath)
    pdbfilepath = pdbdir + 'pdb' + pdbcode + '.ent.gz'
    pdbheader = Bio.PDB.parse_pdb_header(pdbfilepath)
    res = pdbheader['resolution']
    return res

def save_csv(year_range, depositionsperyear, glycansperyear, nglycansperyear, oglycansperyear, cglycansperyear, sglycansperyear, ligandsperyear):
    assert(len(year_range) == len(depositionsperyear))
    assert(len(year_range) == len(glycansperyear))
    assert(len(year_range) == len(nglycansperyear))
    assert(len(year_range) == len(oglycansperyear))
    assert(len(year_range) == len(sglycansperyear))
    assert(len(year_range) == len(cglycansperyear))
    assert(len(year_range) == len(ligandsperyear))

    output_file = "glycosylation_per_year.csv"
    with open(output_file, "w") as output_file: 
        output_file.write("Year,TotalDepositions,TotalGlycosylation,NGlycosylation,OGlycosyation,CGlycosyation,SGlycosyation,Ligands\n")
        for year, depo, total, n, o, c, s, l in zip(year_range, depositionsperyear, glycansperyear,nglycansperyear, oglycansperyear, cglycansperyear, sglycansperyear, ligandsperyear):
            output_file.write(f"{year},{depo},{total},{n},{o},{c},{s},{l}\n")

def save_json(year_range, depositionsperyear, glycansperyear, nglycansperyear, oglycansperyear, cglycansperyear, sglycansperyear, ligandsperyear):
    assert(len(year_range) == len(depositionsperyear))
    assert(len(year_range) == len(glycansperyear))
    assert(len(year_range) == len(nglycansperyear))
    assert(len(year_range) == len(oglycansperyear))
    assert(len(year_range) == len(sglycansperyear))
    assert(len(year_range) == len(cglycansperyear))
    assert(len(year_range) == len(ligandsperyear))

    output_file = "glycosylation_per_year.json"
    data = {}
    for year, depo, total, n, o, c, s, l in zip(year_range, depositionsperyear, glycansperyear,nglycansperyear, oglycansperyear, cglycansperyear, sglycansperyear, ligandsperyear):
        data[str(year)]={"totalDepositions": int(depo), "totalGlycans": int(total) ,"nGlycans": int(n), "oGlycans": int(o), "cGlycans": int(c), "sGlycans": int(s), "ligands": int(l)}
    with open(output_file, "w") as output_file: 
        json.dump(data, output_file)

def glycans_per_year(databasedir, pdbdir):
    print('Analysing database at ' + databasedir)
    filepathlist = file_paths(databasedir)
    glycans = np.zeros(len(filepathlist))
    nglycans = np.zeros(len(filepathlist))
    oglycans = np.zeros(len(filepathlist))
    sglycans = np.zeros(len(filepathlist))
    cglycans = np.zeros(len(filepathlist))
    ligands = np.zeros(len(filepathlist))
    years = np.zeros(len(filepathlist))
    for i in range(len(filepathlist)):
        jsonfile = filepathlist[i]
        years[i] = get_year_pdb(jsonfile,pdbdir) 
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
            if len(sglycan)>0:
                sglycans[i] = 1
            if len(cglycan)>0:
                cglycans[i] = 1
            if len(ligand)>0:
                ligands[i] = 1
            if len(nglycan)>0 or len(oglycan)>0 or len(sglycan)>0 or len(cglycan)>0 or len(ligand)>0:
                glycans[i] = 1
    # years = np.delete(years, np.where(years == 123456789)[0])
    start = np.min(years)
    end = np.max(years)
    year_range = np.arange(start,end, dtype=np.intc)
    nglycansperyear = np.zeros(len(year_range), dtype=np.intc)
    oglycansperyear = np.zeros(len(year_range), dtype=np.intc)
    sglycansperyear = np.zeros(len(year_range), dtype=np.intc)
    cglycansperyear = np.zeros(len(year_range), dtype=np.intc)
    ligandsperyear = np.zeros(len(year_range), dtype=np.intc)
    glycansperyear = np.zeros(len(year_range), dtype=np.intc)
    for i in range(len(year_range)):
        for j in range(len(years)):
            if years[j] == year_range[i]:
                nglycansperyear[i] += nglycans[j]
                oglycansperyear[i] += oglycans[j]
                sglycansperyear[i] += sglycans[j]
                cglycansperyear[i] += cglycans[j]
                ligandsperyear[i] += ligands[j]
                glycansperyear[i] += glycans[j]         
    return year_range, glycansperyear, nglycansperyear, oglycansperyear, sglycansperyear, cglycansperyear, ligandsperyear

def depositions_per_year():
    years = []
    depositions = []
    URL = 'https://www.wwpdb.org/stats/deposition'
    page = requests.get(URL)
    soup = BeautifulSoup(page.content, 'html.parser')
    table = soup.find('h3', string='Number of Structures Released per year').find_next('table') #Finds table under the relevant header
    #tables = soup.find_all('table', attrs={'class':'table table-striped table-bordered text-right'}) #Alternate method for just choosing the second table
    #table = tables[1]
    rows = table.find_all('tr')
    for row in rows:
        data = row.find_all('td')
        data = [ele.text.strip() for ele in data]
        if len(data) > 1: # Get rid of empty row from header
            years.append(int(data[0]))
            depositions.append(int(data[1]))    
    return years, depositions

def plot_and_save_per_year_summary(databasedir,pdbdir):
    year_range, glycansperyear, nglycansperyear, oglycansperyear, sglycansperyear, cglycansperyear, ligandsperyear = glycans_per_year(databasedir, pdbdir)
    years, depositions = depositions_per_year()
    deposperyear = np.zeros(len(year_range))
    for i in range(len(year_range)):
        for j in range(len(years)):
            if years[j] == year_range[i]:
                deposperyear[i] = depositions[j]
    # save_csv(year_range, deposperyear, glycansperyear, nglycansperyear, oglycansperyear, cglycansperyear, sglycansperyear, ligandsperyear)
    save_json(year_range, deposperyear, glycansperyear, nglycansperyear, oglycansperyear, cglycansperyear, sglycansperyear, ligandsperyear)     
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    colour = 'tab:blue'
    ax2.bar(years, depositions, label='Structures Deposited in PDB', color=colour)
    ax2.set_ylabel('Structured Deposited in PDB', color=colour)
    ax1.plot(year_range, glycansperyear, label='Carbohydrate entries in PDB')
    ax1.plot(year_range, nglycansperyear, label='N-glycosylated')
    ax1.plot(year_range, oglycansperyear, label='O-glycosylated')
    ax1.plot(year_range, sglycansperyear, label='S-glycosylated')
    ax1.plot(year_range, cglycansperyear, label='C-glycosylated')
    ax1.plot(year_range, ligandsperyear, label='Ligands')
    ax1.set_xlabel('Release Year')
    ax1.set_ylabel('Carbohydrate Entries in PDB')
    ax1.set_ylim(0,np.max(glycansperyear)+0.01*np.max(glycansperyear))
    ax1.legend()
    plt.tight_layout()
    plt.savefig(cwd+'/glycosylated_per_year.png')
    plt.close()  
    return

plot_and_save_per_year_summary(datadir,pdbdir)
