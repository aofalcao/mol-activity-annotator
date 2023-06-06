#!/usr/bin/python
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.inchi import InchiToInchiKey

import json
import requests

def get_chem_data(smiles):
    try:
        m = Chem.MolFromSmiles(smiles)
    except:
        m = Chem.MolFromSmiles(smiles.strip(), sanitize=False)
        m.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(m,Chem.SANITIZE_SYMMRINGS|Chem.SANITIZE_SETCONJUGATION|Chem.SANITIZE_SETHYBRIDIZATION)

    inchi=Chem.inchi.MolToInchi(m)
    inchikey=InchiToInchiKey(inchi)
    return inchikey

def make_mol(smiles):
    mols={}
    ick=get_chem_data(smiles)
    mols["User1"]=(smiles, ick, "NA")
    return mols
                   


def read_molecules(fname):
    #This must be the structure or it will output an error: (mol_id, activity, SMILES), separated by tabs
    #the activity need not be a float. If char, then it will be a classification model
    fil=open(fname, "rt")
    lins=fil.readlines()
    fil.close()
    
    mols={}
    for lin in lins:   #### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        molid, act, smiles = lin.strip().split("\t")
        try:
            act=float(act)
        except:
            pass
        ick=get_chem_data(smiles)
        mols[molid]=(smiles, ick, act)
    return mols

#this is the heart of the matter
# for one molecule
#first we get to know whether the molecule exists
#then we get the chemblid of the molecule
#then we get the assays
#finally we get the targets

def read_url(url):
    try:
        response = requests.get(url)
        data = json.loads(response.text)
        return data
    except:
        return None

def get_chembl_id(inchikey):
    base_url="https://www.ebi.ac.uk/chembl/api/data/molecule?molecule_structures__standard_inchi_key="
    suffix="&format=json"
    url=base_url + inchikey + suffix
    data=read_url(url)
    if data is None: return False
    if len(data["molecules"])==1:
        return data["molecules"][0]['molecule_chembl_id']
    else:
        return False

def get_activities_chembl(mol_id):
    """
    For a given molecule (mol_id in chembl) will return the possible activities
    """
    base_url="https://www.ebi.ac.uk/chembl/api/data/activity?molecule_chembl_id="
    suffix="&format=json&limit=0"
    url=base_url + mol_id + suffix
    data=read_url(url)
    if data is None: return False
    
    activities=[]
    try:
        for act in data["activities"]:
            D={"assay_id": act["assay_chembl_id"], "target_id": act['target_chembl_id'],
               "assay_type": act['standard_type'], "units":  act['standard_units'], 
               "value": act['standard_value'], "rel":act['standard_relation'], "pvalue": act['pchembl_value']}

            activities.append(D)
    except:
        print("************************************")
        print(data)
        
    return activities
    
def get_targets_chembl(tid):
    base_url = "https://www.ebi.ac.uk/chembl/api/data/target?target_chembl_id="
    suffix = "&format=json"

    url = base_url + tid + suffix
    data=read_url(url)
    if data is None: return False
    
    #print(len(data["targets"]))
    targ=data["targets"][0]
    organism=targ["organism"]
    pname=targ["pref_name"]
    tcomps=[]
    for comp in targ["target_components"]:
        acc="" if "accession" not in comp else comp["accession"]
        gene=""
        if len(comp['target_component_synonyms'])>0:
            for syn in comp["target_component_synonyms"]:
                if syn["syn_type"]=="GENE_SYMBOL":
                    gene=syn["component_synonym"].upper()
        D={"organism":organism, "pname": pname, "uniprot": acc, "gene": gene}
        tcomps.append(D)
    return tcomps
    
def get_target_info_chembl(inchikey):
    mol_id=get_chembl_id(inchikey)
    acts=get_activities_chembl(mol_id)
    targets={}
    for act in acts:
        targs=get_targets_chembl(act["target_id"])
        if len(targs)>= 1: targets[act["target_id"]]=targs
    return targets
    


# a better approach might be:
#1. get all the assays of all molecules
#2. get all the targets 
#3 assemble everything for a nice output

def get_targets_info_chembl(mol_list, data_type="inchikeys", silent=False):
    """
    for each molecule within mol list gets the possible active targets where it has been
    identified as such
    Returns: 
       mol_ids: a dictionnary with inchikey as key and the id as value (only when datatype is "inchikey")
       mol_acts: a dictionary with mol_id as key and a list of targets as value
       btargets: extra information about each target (organism, gene and uniprot_id) 
    """
    mol_ids={}
    mol_acts={}
    if not silent: print("Processing molecules...")
    if data_type=="inchikeys":
        i=0
        for ick in mol_list:
            mol_id=get_chembl_id(ick)
            if mol_id: 
                mol_ids[ick]=mol_id
                acts=get_activities_chembl(mol_id)
                mol_acts[mol_id]=acts
                i+=1
                if i % 100 == 0 and not silent: print("\tProcessed %4d molecules" % i )
    elif data_type=="chemblids":
        i=0
        for mol_id in mol_list:
            acts=get_activities_chembl(mol_id)
            mol_acts[mol_id]=acts
            i+=1
            if i % 100 == 0 and not silent: print("\tProcessed %4d molecules" % i )
    else:
        print("Wrong data type. Exiting")
        return
    #now assemble all the chembl targets found
    ctargets=set()
    for mid in mol_acts:
        acts=mol_acts[mid]
        for act in acts:
            try:
                ctargets.add(act["target_id"])
            except:
                print("----------error---------------->", act)
                pass
    if not silent: print("Getting the Chembl target data...")
    #finally process all targets at once and get the biological targets
    btargets={}
    i=0
    #print("====>", len(ctargets))
    #print(ctargets)

    #this part below is perhaps overkill as we are getting biological information for each target
    #even if such target will never appear on the report. This should be done perhaps ONLY AFTER 
    # the target has passed the filter!
    for ctarg in ctargets: 
        btargs=get_targets_chembl(ctarg)
        for bt in btargs:
            #print("BT:", bt, ctarg)
            btargets.setdefault(ctarg, []).append((bt["organism"], bt["gene"], bt["uniprot"]))
        i+=1
        if i % 100 == 0 and not silent: print("\tProcessed %4d targets" % i )
    if not silent: print("Done!")
    #mol_acts - dic whith key chembl_molid and activities as value
    #btargets - dic with chembl_molid as key and tuple with org gene and uniprot code
    return mol_ids, mol_acts, btargets

def check_activity(act):
    mult_units={"nM": 1, "uM": 1000, "mM":1000000}
    #print(act["assay_type"], "...",act['standard_type'])
    if act["assay_type"] in ["Ki", "IC50", "Kd", "EC50", "AC50", "Potency"]:
        rel=act["rel"]
        unit=act["units"]
        if unit in mult_units: 
            mult=mult_units[unit]
        else:
            #print("Unit not found:", unit, act["assay_type"])
            return "NA"
        try:
            v = float(act["value"])*mult
        except:
            return "N"
        if rel=="=": 
            return "A" if v < 10000  else "N"

        elif rel in [">", ">="]:
            return "NA" if v < 10000  else "N"

        elif rel in ["<", "<="]:
            return "A" if v < 10000  else "NA"
            
    elif act["assay_type"] in ["INH", "Inhibition"]:
        rel=act["rel"]
        try:
            v = float(act["value"])
        except:
            return "N"
        
        return "A" if v>0  else "N"
    
    elif act["assay_type"] in ["pKi", "pIC50"]:
        try:
            v = float(act["value"])
        except:
            return "N"
        return "A" if v < 5.0 else "N"     # 10,000 nM
    elif act["assay_type"] in ["Activity"]:
        try:
            v = float(act["value"])
        except:
            return "N"
        unit=act["units"]
        if unit=="&":
            return "A" if v > 0 else "N"
        else: return "NA"
    else:
        #print("Assay Type Unknown",act['assay_type'], act['value'], act['units'])
        return "NA"
    
    

def gene_counting(mols, mol_active_genes):
    #This is where the report is concluded
    #this is where we count genes
    genes=[]
    uniqs=set()
    for mid in mols:
        if mid in mol_active_genes:
            for g in mol_active_genes[mid]: 
                genes.append(g)
            uniqs=uniqs | mol_active_genes[mid]
    gene_counts=[]
    for g in uniqs: gene_counts.append((genes.count(g), g))
    gene_counts.sort()
    gene_counts.reverse()
    return gene_counts


def write_report_aggregated(mols, mol_active_genes, fname=None):

    if fname is not None: fil=open(fname, "wt")
    gene_counts=gene_counting(mols,mol_active_genes)
    N=len(mols)
    for c, g in gene_counts:
        s="%s\t%5d\t%6.4f" %(g,c,c/N)
        if fname is None: 
            print(s)
        else:
            fil.write(s + "\n")
    if fname is not None: fil.close()
    
def write_report_detail(mols, mol_active_genes, fname=None, append=False):
    #this is where we detail for each molecule the potential genes
    fmode="wt"
    suf=""
    if append==True: 
        fmode="at"
        suf="_"
    if fname is not None: fil=open(fname, fmode)
    for mid in mols:
        s="%s:" % (suf+mid)
        if mid in mol_active_genes:
            for g in mol_active_genes[mid]:
                if len(g.strip())>0: s+= "\t%s" % (g)
        if fname is None: 
            print(s)
        else:
            fil.write(s+"\n")
        
    if fname is not None: fil.close()




def make_reports(mols, mol_ids, mol_acts, targets, search="A", ofield="G"):
    """
     mols is  is a dictionary with the original as keys and with values including fp, inchikey and activity
     mol_ids is a dictionary with the ick as key and value as the chembl_id
     mol_acts is a dictionary with chemble mol_id and an array with all the activities for each chembl target
             the activities are a dictionary with assay_id, target_id, assay_type, units, value, rel
     targets are biological targets discovered: for each chembl target we have an array of biological targets. 
     ofield is a string that can have [O]rganism or [G]ene or [U]niprot or all
     This function returns for each molecule its active genes. Also a set of all targets used (chembl_ids)
     For the moment I'm just going to process the elements with one Biological target
        report specifies which report we will produce (A=Aggregated; D=Detailed)
    """ 
    good_targs=set()

    act_filter=[]        
    if "A" in search: act_filter.append("A")
    if "U" in search: act_filter.append("NA")
    if "N" in search: act_filter.append("N")
    
    mol_active_genes={}
    int_mol_ids=list(mols.keys())
    for mid in int_mol_ids:
        ick = mols[mid][1]
        if ick in mol_ids: 
            db_id =mol_ids[ick]
            for act in mol_acts[db_id]:
                is_active=check_activity(act)
                #print(mid, act["target_id"], is_active)
                if is_active in act_filter:
                    tid=act["target_id"]
                    if tid in targets: 
                        #this is the constraint - we may change it in the future
                        #here we are just processing the targets with ONE gene
                        #if there is more than one, we will not consider them as the essay 
                        #is not specific
                        if len(targets[tid])==1:
                            fs=[]
                            if "O" in ofield: fs.append(targets[tid][0][0].upper()) #organism
                            if "G" in ofield: fs.append(targets[tid][0][1].upper()) #gene
                            if "U" in ofield: fs.append(targets[tid][0][2].upper()) #uniprot
                            output="_".join(fs)
                            #if len(ofield)>1: output="'"+output+"'" #add ticks to aid parsing
                            mol_active_genes.setdefault(mid, set()).add(output)
                            good_targs.add((tid, output))

    return mol_active_genes, good_targs


def make_reports_x(mol_acts, targets, search="A", ofield="G"):
    """
     this is basically the same as make_reports, but simplified when only chembl molecules are being tested
     this is the case when queries from similarities are being processed
     mol_acts is a dictionary with chembl mol_id and an array with all the activities for each chembl target
             the activities are a dictionary with assay_id, target_id, assay_type, units, value, rel
     targets are biological targets discovered: for each chembl target we have an array of biological targets. 
     This function returns for each molecule its active genes. Also a set of all targets used (chembl_ids)
     For the moment I'm just going to process the elements with one Biological target
    """
    good_targs=set()
    act_filter=[]        
    if "A" in search: act_filter.append("A")
    if "U" in search: act_filter.append("NA")
    if "N" in search: act_filter.append("N")
    
    mol_active_genes={}
    for mid in mol_acts:
        for act in mol_acts[mid]:
            is_active = check_activity(act)
            #print(mid, act["target_id"], is_active)
            if is_active in act_filter:
                tid = act["target_id"]
                if tid in targets: 
                    #this is the constraint - we may change it in the future
                    #here we are just processing the targets with ONE gene
                    #if there is more than one, we will not consider them as the essay 
                    #is not specific
                    if len(targets[tid]) == 1:
                        fs=[]
                        if "O" in ofield: fs.append(targets[tid][0][0].upper()) #organism
                        if "G" in ofield: fs.append(targets[tid][0][1].upper()) #gene
                        if "U" in ofield: fs.append(targets[tid][0][2].upper()) #uniprot
                        output="_".join(fs)
                        #if len(ofield)>1: output="'"+output+"'"
                        mol_active_genes.setdefault(mid, set()).add(output)
                        good_targs.add((tid, output))

                        #gene=targets[tid][0][1].upper()
                        #mol_active_genes.setdefault(mid, set()).add(gene)
                        #good_targs.add((tid, gene))
                        ##print(mid, tid, "===>", targets[tid][0][1].upper())
    
    return mol_active_genes, good_targs

def get_sims(smiles, sim_thr=0.9, chembl_server="https://www.ebi.ac.uk"):
    """
    Gets the similar molecules to a given compound by querying Chembl
    """
    sim_mols=[]
    ssim=str(int(sim_thr*100))
    smi1=smiles.replace("%", "%25")
    smi2=smi1.replace("#", "%23")
    smi3=smi2.replace("\\", "%5C")
    smi4=smi3.replace("+", "%2B")
    #https://www.ebi.ac.uk/chembl/api/data/similarity/c1ccccc1CCN/70.xml

    #base_url = "https://www.ebi.ac.uk/chembl/api/data/similarity/"
    #suffix = ".json"
    #url=base_url + smi3 + "/" + ssim + suffix
    url=chembl_server + "/chembl/api/data/similarity.json?smiles="+smi4+"&limit=0&similarity=" + ssim
    
    #print(url)
    data=read_url(url)
    if data is None: return False
    the_mols=data["molecules"]
    while data["page_meta"]["next"] is not None:
        url = chembl_server+data["page_meta"]["next"]
        data=read_url(url)
        the_mols+=data["molecules"]
    if len(the_mols)>0:
        for mol in the_mols:
            sim_mols.append(mol["molecule_chembl_id"])
    return sim_mols

def get_all_sims(mols, thres, silent):
    sim_mols=[]
    i=0
    sim_mols_detail={}
    for mid in mols:
        smi=mols[mid][0]
        #print("****-->", smi)
        sims=get_sims(smi, thres)
        sim_mols += sims
        sim_mols_detail[mid]=sims
        i+=1
        if i % 100 == 0 and not silent: print("\tProcessed %4d similarities" % i )
    #returns a clean set with ALL the distinct molecules from chembl in the neighborhood of the data set
    # and a Dictionary with, for each molecule in the data is present get all its similarities
    #for m in sim_mols_detail:
    #    print("deb:", m, sim_mols_detail[m])
    #print("**************** SIMS ******************")
    #for s in sim_mols: print(s)
    #print("SIZE:", len(sim_mols))
    #print("****************************************")
    return set(sim_mols), sim_mols_detail

def merge_T0T1(mols, mols_t0, mol_genes_T0, mols_t1, mol_genes_T1):
    #mols is the original dic with keys with the original ID
    #mols_T0 is a list with ick as key and the chembl_ids aqs value
    #mol_genes_t0 is the gene data for each molecule of T0 found on Chem
    #mols_t1 is a Dictionary with keys the chemblids and keys the similar ids
    # Deeply fixed (06/06/23) - probably a nasty bug/oversight
    mol_genes_T01={}
    for mid in mols:
        genes=list(mol_genes_T0[mid])
        #if this molecule has chembl neighbours
        if mid in mols_t1:
            adjs=mols_t1[mid]
            for a in adjs:
                if a in mol_genes_T1:
                    genes=genes+list(mol_genes_T1[a])
        mol_genes_T01[mid]=set(genes)
    #returns the genes referenced in T0 and T1, if any
    return mol_genes_T01

def get_activity_profile(fname=None, do_sims=False, thres=0.7, report="A", write_files=True, to_screen=False,
                         join_aggregates=False, silent=False, search="A", single_mol=None, ofield="G"):
    """This is the main entry point.
    fname - file name in .SAR format whith molecules to query
    do_sims - flag that specifies whther or not we will do similarity search
    thres - threshold for similarity search (0.7 - 1.0)
    report - type of report (A - Aggregate; D Detail)
    write_files - flag that indicates whether files are written
    silent - flag that indicates whether the function verbosely explains what is doing
    search - type of search (A- actives; N - non actives; U - unknowns)
    single_mol - if only one molecule is being processed (SMILES or None)
    """
    #this is the main function
    #see if we have to create file names
    if write_files==True:
        #requires a .sar extension for this to work properly
        root=fname.replace(".sar", "")
        detail_0=root+"_annotations_T0.txt"
        detail_1=root+"_annotations_T1.txt"
        aggregate_0=root+"_anno_count_T0.txt"
        aggregate_1=root+"_anno_count_T1.txt"
        if join_aggregates:
            full_aggregate=root+"_anno_count_full.txt"
            #if there is aggregation of data we must also join the details file
            detail_0=root+"_annotations_TA.txt"
            detail_1=root+"_annotations_TA.txt"
                
    else:
        detail_0=None
        detail_1=None
        aggregate_0=None
        aggregate_1=None
        full_aggregate=None

    if single_mol is not None:
        mols=make_mol(single_mol)
    else:
        if fname is not None:
            mols = read_molecules(fname)
        else:
            print("No input file! Nothing to do. Exiting")
            exit()
       
    
    mids=list(mols.keys())
    icks=[]
    smis=[]
    for mid in mids: 
        icks.append(mols[mid][1])
        smis.append(mols[mid][0])

    #
    # this part will get the real molecules out of chembl
    mol_ids, mol_acts, targets = get_targets_info_chembl(icks, silent=silent)

    mol_active_genes, target_ids=make_reports(mols, mol_ids, mol_acts, targets, search=search, ofield=ofield)
    
    if "A" in report and (write_files==True or to_screen==True): 
        write_report_aggregated(mols, mol_active_genes, fname=aggregate_0)

    #this part will extract the similar ones
    #this is the eXtendeded search, thus the 'x'
    #First get the similars for each molecule
    if do_sims==True:
        if not silent: print("Get the similars...", end=" ")
        sim_mols, sim_mols_detail = get_all_sims(mols, thres, silent)
        
        
        # extract the selfs (these were handed down before)
        sim_mols_x = sim_mols - set(mol_ids.values())
        if not silent: print("N. of similar molecules:", len(sim_mols_x))
        if not silent: print("\tNOTE: If the app stops responding for a long time and there is\n\t\tno network activity, please press ENTER")

        sim_mols_x=list(sim_mols_x)
        mol_ids_x, mol_acts_x, targets_x = get_targets_info_chembl(sim_mols_x, "chemblids", silent=silent)
        #mol_ids_x should be empty
        #mol_acts_x should have the activities for each molecule
        # targets, the info about the respective targets
        #at this phase we know everything
        
        #here is philosophical issue that might be solved with more options. But what might be interesting is that 
        #we may always want to have Actives in the similars. Even if the search may include non-actives or unknowns
        #we may look for only situations for which we know the similar molecule is active in. [eventually include unknowns]
        #so the following line will be replaced by the next one, where the search must select only ACTIVE similars
        #mol_active_genes_x, target_ids_x = make_reports_x(mol_acts_x, targets_x, search=search)
        mol_active_genes_x, target_ids_x = make_reports_x(mol_acts_x, targets_x, search="A", ofield=ofield)
        
        #target_ids contains a set of all valid targets (chembl_ids) found
        #target_ids.update(target_ids_x)
        
        if "A" in report  and (write_files==True or to_screen==True):
            write_report_aggregated(sim_mols_x, mol_active_genes_x, fname=aggregate_1)

        genes_T01 = merge_T0T1(mols, mol_ids, mol_active_genes, sim_mols_detail, mol_active_genes_x)
        
        if "D" in report and (write_files==True or to_screen==True): 
            write_report_detail(mols, mol_active_genes, fname=detail_0)
            write_report_detail(mols, genes_T01,        fname=detail_1)
    
        return {"mols": mols, "active_genes": mol_active_genes, 
                "mols_x": sim_mols_x, "active_genes_x": mol_active_genes_x, "target_ids": target_ids,
                "target_ids_x": target_ids_x}
    else:
        if "D" in report  and (write_files==True or to_screen==True): 
            write_report_detail(mols, mol_active_genes, fname=detail_0)
        
        return {"mols": mols, "active_genes": mol_active_genes,
                "mols_x": None, "active_genes_x": None, "target_ids": target_ids,
                "target_ids_x": None}
    

    
#def write_report_aggregated(mols, mol_active_genes, fname=None):    
#def write_report_detail(mols, mol_active_genes, fname=None, append=False):


def print_help():
    s= """
Usage: This is a Python3 tool and requires an enviroment where RDkit and requests are installed.
       To run type python csanno.py  -in [input .sar file] [options] 

Control parameters:
    -in file_name - the data set to annotate (required) (.sar format)
    -sim thr - similarity threshold [if absent (default) no similarity search will be performed]
    -mol [molecule in SMILES format] - molecule for which analysis is going to be performed.  
         Incompatible with -in option
    -report [AD]: A - aggregate report;
                  D - Detailed report (one line for each molecule)
                  Both options simultaneusly are permitted and both reports are produced
    -search [ANU] A - searches for actives in the database (default);
                  N - searches non-actives;
                  U - searches unknowns (compounds for which activity is not determined)
                  All options can be included simultaneously 
    -ofield [GUO] G - outputs the unique gene (default);
                  U - outputs the uniprot id;
                  O - outputs the organism
                  All options can be included simultaneously 

Output Control Options:
    -joint_aggs - produces a report with conjoined similarities and exact matches
    -silent - no intermediate output at all
    -help - this screen
    -nofiles - no files are created
    -toscreen - writes output to screen
    -pickle - writes a Python pickle for easy posterior analysis (see internal documentation)
"""
    print(version)
    print(s.strip())
    

version="CSANNO - (C) 2019/2023 - Andre O. Falcao DI/FCUL version 0.3.20230605"
if __name__=="__main__":
    import sys
    import time
    start_t=time.time()
    if len(sys.argv)< 2:
        print("Invalid option. Exiting")
        exit()

    #defaults
    report = "A"  #default is only the aggregates
    search = "A"  #default only actives 
    ofield="G"    #default, only outputs the Gene
    do_sims = False
    joint_aggs = False
    sim_thr = 0.9
    writefiles = True
    silent = False
    writepickle = False
    to_screen=False
    mol=None
    input_file=None
    

    
    bin_args=["-in", "-sim", "-search", "-report", "-mol", "-ofield"]
    all_args=bin_args[:]+["-silent", "-help", "-nofiles", "-pickle", "-joint_aggs", "-toscreen"]

    inp_fs=sys.argv
    method = inp_fs[1]

    if '-help' in inp_fs or '--help' in inp_fs:
        print_help()
        exit()

    #error checker
    for arg in inp_fs:
        if arg[0]=="-" and arg not in all_args:
            print("'%s' Unknown Option. Exiting" % arg)
            exit()
	    
    try:
        #arg_screener
        input_args={}
        for arg in bin_args:
            if arg in inp_fs:
                pos_val = inp_fs.index(arg)+1
                input_args[arg]=inp_fs[pos_val]


        if '-in'  in input_args:
            input_file  = input_args["-in"]
        if '-sim' in input_args:
            sim_thr  = float(input_args["-sim"])
            do_sims = True
        if '-mol' in input_args:
            mol  = input_args["-mol"]
            silent = True
            writefiles = False
            to_screen = True
            input_file = None  #if -mol we will ignore any input file
        if '-search' in input_args:  search  = input_args["-search"]
        if '-report' in input_args:  report  = input_args["-report"]
        if '-ofield' in input_args:  ofield  = input_args["-ofield"]
        
        if '-nofiles' in inp_fs:  writefiles = False
        if '-toscreen' in inp_fs:  to_screen = True
        if '-silent' in inp_fs: silent = True
        if '-pickle' in inp_fs: writepickle = True
        if '-joint_aggs' in inp_fs: joint_aggs = True
        
        
    except:
        print("Illegal Option or Value. Exiting")
        exit()

    if input_file is None and mol is None: 
        print("No input file! Nothing to do. Exiting")
        exit()

    if len(set(search) - set("AUN") )>0 or len(search)==0 or len(search)>3:
        print("Invalid search option: %s  - Exiting" % search)
        exit()

    if len(set(ofield) - set("OGU") )>0 or len(ofield)==0 or len(ofield)>3:
        print("Invalid output field option: %s  - Exiting" % ofield)
        exit()


    if report not in ["A", "D", "AD"]:
        print("Invalid report option: %s - Exiting" % report)
        exit()

    if writefiles==False and to_screen==False: 
        print("No files and no screen output. Nothing to see! Exiting")
        exit()

    #print(input_file, do_sims, sim_thr, report, search, joint_aggs)
    if not silent: print(version)
    if not silent: print("\tNOTE: If the app stops responding for a long time and there is\n\t\tno network activity, please press ENTER")
    rep=get_activity_profile(fname=input_file, do_sims=do_sims, thres=sim_thr, report=report,
                             search=search, join_aggregates=joint_aggs, write_files=writefiles, to_screen=to_screen,
                             silent=silent, single_mol=mol, ofield=ofield)

    if writepickle==True:
        import pickle
        fname=input_file.replace(".sar", "")+"_annotations.pickle"
        pickle.dump(rep, open(fname, "wb"))
    #print("%9.4f secs" % (time.time() - start_t))
