import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import json
from datetime import datetime
import requests
from bs4 import BeautifulSoup
import gzip

import gemmi

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


def get_year_pdb(jsonfilepath):
    jsonfilename = os.path.basename(jsonfilepath)
    pdbcode = jsonfilename.rpartition('.')[0]
    try:
        pdbfilepath = get_pdb_path(pdbcode)
        #print(f'Looking for pdb at {pdbfilepath}')
        st = gemmi.read_structure(pdbfilepath)
        metadata = st.info
        year = metadata['_pdbx_database_status.recvd_initial_deposition_date'].split('-')[0]
    except:
        try:
            ciffilepath = f'/vault/pdb_mirror/data/structures/all/mmCIF/{pdbcode}.cif.gz'
            print(f'Failed to find corresponding pdb to {pdbcode} looking for mmCIF at {ciffilepath}')
            st = gemmi.read_structure(ciffilepath)
            metadata = st.info
            year = metadata['_pdbx_database_status.recvd_initial_deposition_date'].split('-')[0]
        except:
            print(f'Failed to find corresponding mmCIF to {pdbcode}')
            year = 123456789
    return int(year)

def save_csv(year_range, depositionsperyear, glycansperyear, nglycansperyear, oglycansperyear, cglycansperyear, sglycansperyear, ligandsperyear, redo):
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

def save_json(year_range, depositionsperyear, glycansperyear, nglycansperyear, oglycansperyear, cglycansperyear, sglycansperyear, ligandsperyear, redo):
    assert(len(year_range) == len(depositionsperyear))
    assert(len(year_range) == len(glycansperyear))
    assert(len(year_range) == len(nglycansperyear))
    assert(len(year_range) == len(oglycansperyear))
    assert(len(year_range) == len(sglycansperyear))
    assert(len(year_range) == len(cglycansperyear))
    assert(len(year_range) == len(ligandsperyear))

    if redo:
        output_file = "glycosylation_per_year_redo.json"
    else:
        output_file = "glycosylation_per_year.json"
    data = {}
    for year, depo, total, n, o, c, s, l in zip(year_range, depositionsperyear, glycansperyear,nglycansperyear, oglycansperyear, cglycansperyear, sglycansperyear, ligandsperyear):
        data[str(year)]={"totalDepositions": int(depo), "totalGlycans": int(total) ,"nGlycans": int(n), "oGlycans": int(o), "cGlycans": int(c), "sGlycans": int(s), "ligands": int(l)}
    with open(output_file, "w") as output_file: 
        json.dump(data, output_file)

def glycans_per_year(databasedir):
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
            if len(sglycan)>0:
                sglycans[i] = 1
            if len(cglycan)>0:
                cglycans[i] = 1
            if len(ligand)>0:
                ligands[i] = 1
            if len(nglycan)>0 or len(oglycan)>0 or len(sglycan)>0 or len(cglycan)>0 or len(ligand)>0:
                glycans[i] = 1
    years = np.delete(years, np.where(years == 123456789)[0])
    start = np.min(years)
    end = np.max(years)
    year_range = np.arange(start,end+1, dtype=np.intc)
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

def plot_and_save_per_year_summary(databasedir,redo):
    year_range, glycansperyear, nglycansperyear, oglycansperyear, sglycansperyear, cglycansperyear, ligandsperyear = glycans_per_year(databasedir)
    years, depositions = depositions_per_year()
    deposperyear = np.zeros(len(year_range))
    for i in range(len(year_range)):
        for j in range(len(years)):
            if years[j] == year_range[i]:
                deposperyear[i] = depositions[j]
    # save_csv(year_range, deposperyear, glycansperyear, nglycansperyear, oglycansperyear, cglycansperyear, sglycansperyear, ligandsperyear, redo)
    save_json(year_range, deposperyear, glycansperyear, nglycansperyear, oglycansperyear, cglycansperyear, sglycansperyear, ligandsperyear, redo)     
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    colour = 'tab:blue'
    ax2.bar(year_range, deposperyear, label='Structures Deposited in PDB', color=colour, alpha=0.5)
    ax2.set_ylabel('Structured Deposited in PDB', color=colour)
    ax2.set_ylim(0,np.max(depositions)+0.01*np.max(depositions))
    ax2.set_xlim(np.min(year_range),np.max(year_range))
    ax1.plot(year_range, glycansperyear, color='tab:orange', label='Carbohydrate entries in PDB')
    ax1.plot(year_range, nglycansperyear, color='tab:green', label='N-glycosylated')
    ax1.plot(year_range, oglycansperyear, color='tab:red', label='O-glycosylated')
    ax1.plot(year_range, sglycansperyear, color='tab:purple', label='S-glycosylated')
    ax1.plot(year_range, cglycansperyear, color='tab:brown', label='C-glycosylated')
    ax1.plot(year_range, ligandsperyear, color='tab:pink', label='Ligands')
    ax1.set_xlabel('Release Year')
    ax1.set_ylabel('Carbohydrate Entries in PDB')
    ax1.set_ylim(0,np.max(glycansperyear)+0.01*np.max(glycansperyear))
    ax1.legend()
    plt.tight_layout()
    if redo:
        plt.savefig('glycosylated_per_year_redo.png')
    else:
        plt.savefig('glycosylated_per_year.png')
    plt.close()  
    return

if __name__ == "__main__":
    redo = True
    if redo:
        datadir = '/vault/privateer_database/pdbredo/'
    else:
        datadir = '/vault/privateer_database/pdb/'
    plot_and_save_per_year_summary(datadir, redo)