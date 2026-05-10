#!/usr/bin/python

#Chemical Space Annotation - csanno

#Copyright © 2022-2026 Andre O. Falcao - BioISI - DI/FCUL

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software
#and associated documentation files (the "Software"), to deal in the Software without restriction, 
#including without limitation the rights to use, copy, modify, merge, publish, distribute, 
#publicense, and/or sell copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or 
#substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
#PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
#FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
#OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
#DEALINGS IN THE SOFTWARE.

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.inchi import InchiToInchiKey
import pickle
import os

import json
import requests

DEFAULT_CHEMBL_SERVER = "https://www.ebi.ac.uk"
REQUEST_TIMEOUT = 30
TARGET_CACHE_JSON = "./btargets_cache.json"
TARGET_CACHE_PICKLE = "./btargets.pickle"


def warn(msg, silent=False):
    if not silent:
        print("WARNING: %s" % msg)


def normalise_chembl_server(chembl_server):
    if chembl_server is None or str(chembl_server).strip() == "":
        return DEFAULT_CHEMBL_SERVER
    return str(chembl_server).rstrip("/")


def safe_float(x):
    try:
        if x is None:
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def get_chem_data(smiles):
    if smiles is None:
        raise ValueError("Missing SMILES")

    smiles = str(smiles).strip()
    if smiles == "":
        raise ValueError("Empty SMILES")

    m = Chem.MolFromSmiles(smiles)
    if m is None:
        try:
            m = Chem.MolFromSmiles(smiles, sanitize=False)
            if m is not None:
                m.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(m, Chem.SANITIZE_SYMMRINGS |
                                    Chem.SANITIZE_SETCONJUGATION |
                                    Chem.SANITIZE_SETHYBRIDIZATION)
        except Exception as exc:
            raise ValueError("Invalid SMILES '%s': %s" % (smiles, exc))

    if m is None:
        raise ValueError("Invalid SMILES '%s'" % smiles)

    inchi = Chem.inchi.MolToInchi(m)
    inchikey = InchiToInchiKey(inchi)
    if inchikey is None or inchikey == "":
        raise ValueError("Could not generate InChIKey for SMILES '%s'" % smiles)
    return inchikey


def make_mol(smiles):
    mols = {}
    ick = get_chem_data(smiles)
    mols["User1"] = (smiles, ick, "NA")
    return mols
                   


def read_molecules(fname, silent=False):
    #This must be the structure: (mol_id, activity, SMILES), separated by tabs
    #the activity need not be a float. If char, then it will be a classification model
    #blank lines and lines starting with # are ignored
    if fname is None or str(fname).strip() == "":
        raise ValueError("No input file provided")

    if not os.path.isfile(fname):
        raise ValueError("Input file not found: %s" % fname)

    mols = {}
    with open(fname, "rt") as fil:
        for lno, lin in enumerate(fil, 1):
            lin = lin.rstrip("\n")
            if lin.strip() == "" or lin.lstrip().startswith("#"):
                continue

            fields = lin.split("\t")
            if len(fields) != 3:
                warn("Skipping line %d: expected 3 tab-separated fields, found %d" % (lno, len(fields)), silent)
                continue

            molid, act, smiles = fields
            molid = molid.strip()
            smiles = smiles.strip()

            if molid == "":
                warn("Skipping line %d: empty molecule id" % lno, silent)
                continue

            try:
                act = float(act)
            except (TypeError, ValueError):
                act = act.strip()

            try:
                ick = get_chem_data(smiles)
            except ValueError as exc:
                warn("Skipping line %d (%s): %s" % (lno, molid, exc), silent)
                continue

            if molid in mols:
                warn("Duplicate molecule id '%s' on line %d; replacing previous entry" % (molid, lno), silent)

            mols[molid] = (smiles, ick, act)

    if len(mols) == 0:
        raise ValueError("No valid molecules were read from %s" % fname)

    return mols


#this is the heart of the matter
# for one molecule
#first we get to know whether the molecule exists
#then we get the chemblid of the molecule
#then we get the assays
#finally we get the targets


def read_url(url, params=None, timeout=REQUEST_TIMEOUT, silent=True):
    try:
        response = requests.get(url, params=params, timeout=timeout)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.Timeout:
        warn("Request timed out: %s" % url, silent)
        return None
    except requests.exceptions.RequestException as exc:
        warn("Request failed for %s: %s" % (url, exc), silent)
        return None
    except ValueError:
        warn("Could not decode JSON response from %s" % url, silent)
        return None


def load_target_cache(silent=False):
    if os.path.isfile(TARGET_CACHE_JSON):
        try:
            with open(TARGET_CACHE_JSON, "rt") as fil:
                data = json.load(fil)
            if isinstance(data, dict):
                return data
        except (IOError, ValueError) as exc:
            warn("Could not read target JSON cache: %s" % exc, silent)

    # Pickle is intentionally not loaded by default because loading an untrusted pickle
    # can execute arbitrary code. Set CSANNO_TRUST_PICKLE_CACHE=1 if this old cache file
    # is trusted and you want to reuse it.
    if os.path.isfile(TARGET_CACHE_PICKLE):
        if os.environ.get("CSANNO_TRUST_PICKLE_CACHE") == "1":
            try:
                with open(TARGET_CACHE_PICKLE, "rb") as fil:
                    data = pickle.load(fil)
                if isinstance(data, dict):
                    warn("Loaded legacy pickle cache because CSANNO_TRUST_PICKLE_CACHE=1", silent)
                    return data
            except Exception as exc:
                warn("Could not read legacy pickle cache: %s" % exc, silent)
        else:
            warn("Ignoring legacy btargets.pickle cache. Set CSANNO_TRUST_PICKLE_CACHE=1 to trust it.", silent)

    return {}


def save_target_cache(btargets, silent=False):
    try:
        with open(TARGET_CACHE_JSON, "wt") as fil:
            json.dump(btargets, fil)
    except IOError as exc:
        warn("Could not write target JSON cache: %s" % exc, silent)


def get_chembl_id(inchikey, chembl_server=DEFAULT_CHEMBL_SERVER, silent=False):
    chembl_server = normalise_chembl_server(chembl_server)
    url = chembl_server + "/chembl/api/data/molecule.json"
    params = {"molecule_structures__standard_inchi_key": inchikey, "format": "json"}
    data = read_url(url, params=params, silent=silent)
    if data is None:
        return False

    molecules = data.get("molecules", [])
    if len(molecules) == 1:
        return molecules[0].get('molecule_chembl_id', False)
    else:
        return False


def get_activities_chembl(mol_id, chembl_server=DEFAULT_CHEMBL_SERVER, silent=False):
    """
    For a given molecule (mol_id in chembl) will return the possible activities
    """
    if mol_id is None or mol_id is False:
        return []

    chembl_server = normalise_chembl_server(chembl_server)
    url = chembl_server + "/chembl/api/data/activity.json"
    params = {"molecule_chembl_id": mol_id, "format": "json", "limit": 0}
    data = read_url(url, params=params, silent=silent)
    if data is None:
        return []
    
    activities = []
    for act in data.get("activities", []):
        try:
            D = {"assay_id": act.get("assay_chembl_id"),
                 "target_id": act.get('target_chembl_id'),
                 "assay_type": act.get('standard_type'),
                 "units": act.get('standard_units'), 
                 "value": act.get('standard_value'),
                 "rel": act.get('standard_relation'),
                 "pvalue": act.get('pchembl_value')}
            if D["target_id"] is not None:
                activities.append(D)
        except (AttributeError, KeyError) as exc:
            warn("Skipping malformed activity for %s: %s" % (mol_id, exc), silent)

    return activities
    

def get_targets_chembl(tid, chembl_server=DEFAULT_CHEMBL_SERVER, silent=False):
    #gets the target information
    # this could perhaps be sped up if we keep a cache of already queried tids
    if tid is None or tid is False:
        return []

    chembl_server = normalise_chembl_server(chembl_server)
    url = chembl_server + "/chembl/api/data/target.json"
    params = {"target_chembl_id": tid, "format": "json"}
    data = read_url(url, params=params, silent=silent)
    if data is None:
        return []
    
    targets = data.get("targets", [])
    if len(targets) == 0:
        return []

    targ = targets[0]
    organism = targ.get("organism", "")
    pname = targ.get("pref_name", "")
    tcomps = []
    for comp in targ.get("target_components", []):
        acc = comp.get("accession", "")
        gene = ""
        for syn in comp.get("target_component_synonyms", []):
            if syn.get("syn_type") == "GENE_SYMBOL":
                gene = syn.get("component_synonym", "").upper()
        D = {"organism": organism, "pname": pname, "uniprot": acc, "gene": gene}
        tcomps.append(D)
    return tcomps
    

def get_target_info_chembl(inchikey, chembl_server=DEFAULT_CHEMBL_SERVER, silent=False):
    mol_id = get_chembl_id(inchikey, chembl_server=chembl_server, silent=silent)
    if not mol_id:
        return {}
    acts = get_activities_chembl(mol_id, chembl_server=chembl_server, silent=silent)
    targets = {}
    for act in acts:
        targs = get_targets_chembl(act.get("target_id"), chembl_server=chembl_server, silent=silent)
        if len(targs) >= 1:
            targets[act["target_id"]] = targs
    return targets
    


# a better approach might be:
#1. get all the assays of all molecules
#2. get all the targets 
#3 assemble everything for a nice output


def get_targets_info_chembl(mol_list, data_type="inchikeys", silent=False, chembl_server=DEFAULT_CHEMBL_SERVER):
    """
    for each molecule within mol list gets the possible active targets where it has been
    identified as such
    Returns: 
       mol_ids: a dictionary with inchikey as key and the id as value (only when datatype is "inchikey")
       mol_acts: a dictionary with mol_id as key and a list of targets as value
       btargets: extra information about each target (organism, gene and uniprot_id) 
    """
    mol_ids = {}
    mol_acts = {}
    if not silent:
        print("Processing molecules...")

    if data_type == "inchikeys":
        i = 0
        for ick in mol_list:
            mol_id = get_chembl_id(ick, chembl_server=chembl_server, silent=silent)
            if mol_id: 
                mol_ids[ick] = mol_id
                acts = get_activities_chembl(mol_id, chembl_server=chembl_server, silent=silent)
                mol_acts[mol_id] = acts
                i += 1
                if i % 100 == 0 and not silent:
                    print("\tProcessed %4d molecules" % i)
    elif data_type == "chemblids":
        i = 0
        for mol_id in mol_list:
            acts = get_activities_chembl(mol_id, chembl_server=chembl_server, silent=silent)
            mol_acts[mol_id] = acts
            i += 1
            if i % 100 == 0 and not silent:
                print("\tProcessed %4d molecules" % i)
    else:
        print("Wrong data type. Exiting")
        return {}, {}, {}

    #now assemble all the chembl targets found
    ctargets = set()
    for mid in mol_acts:
        acts = mol_acts[mid]
        for act in acts:
            tid = act.get("target_id")
            if tid is not None:
                ctargets.add(tid)

    if not silent:
        print("Getting the Chembl target data...")

    #finally process all targets at once and get the biological targets
    i = 0

    #this part below is perhaps overkill as we are getting biological information for each target
    #even if such target will never appear on the report. This should be done perhaps ONLY AFTER 
    # the target has passed the filter!

    btargets = load_target_cache(silent=silent)
    update_cache = False
    for ctarg in ctargets: 
        if ctarg in btargets: 
            continue

        update_cache = True
        if not silent:
            print("\tthis one is new! ->", ctarg)
        btargs = get_targets_chembl(ctarg, chembl_server=chembl_server, silent=silent)
        btargets[ctarg] = []
        for bt in btargs:
            btargets[ctarg].append((bt.get("organism", ""), bt.get("gene", ""), bt.get("uniprot", "")))
        #just going to count those targets for which we need to retrieve them from chembl
        i += 1
        if i % 100 == 0 and not silent:
            print("\tProcessed %4d targets" % i)

    #if any change in btargets we will write them
    if update_cache == True:
        save_target_cache(btargets, silent=silent)
    if not silent:
        print("Done!")
        
    #mol_acts - dic whith key chembl_molid and activities as value
    #btargets - dic with chembl target id as key and tuple with org gene and uniprot code
    return mol_ids, mol_acts, btargets


def check_activity(act):
    mult_units = {"nM": 1, "uM": 1000, "µM": 1000, "mM": 1000000}
    assay_type = act.get("assay_type")
    rel = act.get("rel")
    unit = act.get("units")

    if rel is None or rel == "":
        rel = "="

    # Standard concentration readouts. Values are converted to nM.
    if assay_type in ["Ki", "IC50", "Kd", "EC50", "AC50", "Potency"]:
        if unit in mult_units: 
            mult = mult_units[unit]
        else:
            return "NA"

        v = safe_float(act.get("value"))
        if v is None:
            return "N"
        v = v * mult

        if rel == "=": 
            return "A" if v < 10000 else "N"

        elif rel in [">", ">="]:
            return "NA" if v < 10000 else "N"

        elif rel in ["<", "<="]:
            return "A" if v < 10000 else "NA"

        else:
            return "NA"
            
    elif assay_type in ["INH", "Inhibition"]:
        v = safe_float(act.get("value"))
        if v is None:
            return "N"

        if rel in ["=", ">", ">="]:
            return "A" if v > 0 else "N"
        elif rel in ["<", "<="]:
            return "NA" if v > 0 else "N"
        else:
            return "NA"
    
    elif assay_type in ["pKi", "pIC50"]:
        v = safe_float(act.get("value"))
        if v is None:
            return "N"

        # pKi/pIC50 are -log10(M). A value of 5 corresponds to 10,000 nM.
        # Higher p-values mean stronger activity, so this is the opposite direction
        # from raw Ki/IC50 concentration values.
        if rel == "=":
            return "A" if v > 5.0 else "N"
        elif rel in [">", ">="]:
            return "A" if v > 5.0 else "NA"
        elif rel in ["<", "<="]:
            return "NA" if v > 5.0 else "N"
        else:
            return "NA"

    elif assay_type in ["Activity"]:
        v = safe_float(act.get("value"))
        if v is None:
            return "N"
        if unit == "&":
            return "A" if v > 0 else "N"
        else:
            return "NA"
    else:
        return "NA"
    
    

def gene_counting(mols, mol_active_genes):
    #This is where the report is concluded
    #this is where we count genes
    genes = []
    uniqs = set()
    for mid in mols:
        if mid in mol_active_genes:
            for g in mol_active_genes[mid]: 
                genes.append(g)
            uniqs = uniqs | mol_active_genes[mid]
    gene_counts = []
    for g in uniqs:
        gene_counts.append((genes.count(g), g))
    gene_counts.sort(key=lambda x: (-x[0], x[1]))
    return gene_counts



def write_report_aggregated(mols, mol_active_genes, fname=None):
    fil = None
    if fname is not None:
        fil = open(fname, "wt")

    gene_counts = gene_counting(mols, mol_active_genes)
    N = len(mols)
    if N == 0:
        if fil is not None:
            fil.close()
        return

    for c, g in gene_counts:
        if len(str(g).strip()) == 0:
            continue
        s = "%s\t%5d\t%6.4f" % (g, c, c/N)
        if fname is None: 
            print(s)
        else:
            fil.write(s + "\n")
    if fil is not None:
        fil.close()
    

def write_report_detail(mols, mol_active_genes, fname=None, append=False):
    #this is where we detail for each molecule the potential genes
    fmode = "wt"
    suf = ""
    if append == True: 
        fmode = "at"
        suf = "_"

    fil = None
    if fname is not None:
        fil = open(fname, fmode)

    for mid in mols:
        s = "%s:" % (suf + mid)
        if mid in mol_active_genes:
            for g in sorted(mol_active_genes[mid]):
                if len(str(g).strip()) > 0:
                    s += "\t%s" % (g)
        if fname is None: 
            print(s)
        else:
            fil.write(s + "\n")
        
    if fil is not None:
        fil.close()




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
    good_targs = set()

    act_filter = []        
    if "A" in search:
        act_filter.append("A")
    if "U" in search:
        act_filter.append("NA")
    if "N" in search:
        act_filter.append("N")
    
    mol_active_genes = {}
    int_mol_ids = list(mols.keys())
    for mid in int_mol_ids:
        ick = mols[mid][1]
        if ick in mol_ids: 
            db_id = mol_ids[ick]
            for act in mol_acts.get(db_id, []):
                is_active = check_activity(act)
                if is_active in act_filter:
                    tid = act.get("target_id")
                    if tid in targets: 
                        #this is the constraint - we may change it in the future
                        #here we are just processing the targets with ONE gene
                        #if there is more than one, we will not consider them as the essay 
                        #is not specific
                        if len(targets[tid]) == 1:
                            fs = []
                            if "O" in ofield:
                                fs.append(str(targets[tid][0][0]).upper()) #organism
                            if "G" in ofield:
                                fs.append(str(targets[tid][0][1]).upper()) #gene
                            if "U" in ofield:
                                fs.append(str(targets[tid][0][2]).upper()) #uniprot
                            output = "_".join(fs)
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
    good_targs = set()
    act_filter = []        
    if "A" in search:
        act_filter.append("A")
    if "U" in search:
        act_filter.append("NA")
    if "N" in search:
        act_filter.append("N")
    
    mol_active_genes = {}
    for mid in mol_acts:
        for act in mol_acts[mid]:
            is_active = check_activity(act)
            if is_active in act_filter:
                tid = act.get("target_id")
                if tid in targets: 
                    #this is the constraint - we may change it in the future
                    #here we are just processing the targets with ONE gene
                    #if there is more than one, we will not consider them as the essay 
                    #is not specific
                    if len(targets[tid]) == 1:
                        fs = []
                        if "O" in ofield:
                            fs.append(str(targets[tid][0][0]).upper()) #organism
                        if "G" in ofield:
                            fs.append(str(targets[tid][0][1]).upper()) #gene
                        if "U" in ofield:
                            fs.append(str(targets[tid][0][2]).upper()) #uniprot
                        output = "_".join(fs)
                        mol_active_genes.setdefault(mid, set()).add(output)
                        good_targs.add((tid, output))
    
    return mol_active_genes, good_targs


def get_sims(smiles, sim_thr=0.9, chembl_server=DEFAULT_CHEMBL_SERVER, silent=False):
    """
    Gets the similar molecules to a given compound by querying Chembl
    """
    sim_mols = []
    ssim = str(int(sim_thr * 100))
    chembl_server = normalise_chembl_server(chembl_server)

    url = chembl_server + "/chembl/api/data/similarity.json"
    params = {"smiles": smiles, "limit": 0, "similarity": ssim}
    data = read_url(url, params=params, silent=silent)
    if data is None:
        return []

    the_mols = data.get("molecules", [])
    page_meta = data.get("page_meta", {})
    while page_meta.get("next") is not None:
        next_url = chembl_server + page_meta.get("next")
        data = read_url(next_url, silent=silent)
        if data is None:
            break
        the_mols += data.get("molecules", [])
        page_meta = data.get("page_meta", {})

    if len(the_mols) > 0:
        for mol in the_mols:
            mid = mol.get("molecule_chembl_id")
            if mid is not None:
                sim_mols.append(mid)
    return sim_mols


def get_all_sims(mols, thres, silent, chembl_server=DEFAULT_CHEMBL_SERVER):
    sim_mols = []
    i = 0
    sim_mols_detail = {}
    for mid in mols:
        smi = mols[mid][0]
        sims = get_sims(smi, thres, chembl_server=chembl_server, silent=silent)
        sim_mols += sims
        sim_mols_detail[mid] = sims
        i += 1
        if i % 100 == 0 and not silent:
            print("\tProcessed %4d similarities" % i)
    #returns a clean set with ALL the distinct molecules from chembl in the neighborhood of the data set
    # and a Dictionary with, for each molecule in the data is present get all its similarities
    return set(sim_mols), sim_mols_detail


def merge_T0T1(mols, mols_t0, mol_genes_T0, mols_t1, mol_genes_T1):
    #mols is the original dic with keys with the original ID
    #mols_T0 is a list with ick as key and the chembl_ids aqs value
    #mol_genes_t0 is the gene data for each molecule of T0 found on Chem
    #mols_t1 is a Dictionary with keys the input ids and values the similar chembl ids
    # Deeply fixed (06/06/23) - probably a nasty bug/oversight
    mol_genes_T01 = {}
    for mid in mols:
        genes = []
        if mid in mol_genes_T0:
            genes = list(mol_genes_T0[mid])
        #if this molecule has chembl neighbours
        if mid in mols_t1:
            adjs = mols_t1[mid]
            for a in adjs:
                if a in mol_genes_T1:
                    genes = genes + list(mol_genes_T1[a])
        mol_genes_T01[mid] = set(genes)
    #returns the genes referenced in T0 and T1, if any
    return mol_genes_T01


def merge_T1(mols, mols_t1, mol_genes_T1):
    #returns only the genes referenced in T1, grouped by the original molecule id
    mol_genes_T1_by_input = {}
    for mid in mols:
        genes = []
        if mid in mols_t1:
            adjs = mols_t1[mid]
            for a in adjs:
                if a in mol_genes_T1:
                    genes = genes + list(mol_genes_T1[a])
        mol_genes_T1_by_input[mid] = set(genes)
    return mol_genes_T1_by_input


def write_reports_for_channel(mols, mol_active_genes, report, fname_aggregate, fname_detail, write_files, to_screen):
    if "A" in report:
        if write_files == True:
            write_report_aggregated(mols, mol_active_genes, fname=fname_aggregate)
        if to_screen == True:
            write_report_aggregated(mols, mol_active_genes, fname=None)

    if "D" in report:
        if write_files == True:
            write_report_detail(mols, mol_active_genes, fname=fname_detail)
        if to_screen == True:
            write_report_detail(mols, mol_active_genes, fname=None)


def get_activity_profile(fname=None, do_sims=False, thres=0.7, report="A", write_files=True, to_screen=False,
                         join_aggregates=False, silent=False, search="A", single_mol=None, ofield="G",
                         chembl_server=DEFAULT_CHEMBL_SERVER):
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
    full_aggregate = None
    detail_A = None
    if write_files == True:
        #requires a .sar extension for this to work properly
        if fname is not None:
            root = fname.replace(".sar", "")
        else:
            root = "csanno_single_mol"
        detail_0 = root + "_annotations_T0.txt"
        detail_1 = root + "_annotations_T1.txt"
        detail_A = root + "_annotations_TA.txt"
        aggregate_0 = root + "_anno_count_T0.txt"
        aggregate_1 = root + "_anno_count_T1.txt"
        full_aggregate = root + "_anno_count_full.txt"
    else:
        detail_0 = None
        detail_1 = None
        detail_A = None
        aggregate_0 = None
        aggregate_1 = None
        full_aggregate = None

    if single_mol is not None:
        mols = make_mol(single_mol)
    else:
        if fname is not None:
            mols = read_molecules(fname, silent=silent)
        else:
            print("No input file! Nothing to do. Exiting")
            exit()
       
    mids = list(mols.keys())
    icks = []
    smis = []
    for mid in mids: 
        icks.append(mols[mid][1])
        smis.append(mols[mid][0])

    #
    # this part will get the real molecules out of chembl

    mol_ids, mol_acts, targets = get_targets_info_chembl(icks, silent=silent, chembl_server=chembl_server)

    mol_active_genes, target_ids = make_reports(mols, mol_ids, mol_acts, targets, search=search, ofield=ofield)
    
    if do_sims == True:
        if "A" in report:
            if write_files == True:
                write_report_aggregated(mols, mol_active_genes, fname=aggregate_0)
            if to_screen == True:
                write_report_aggregated(mols, mol_active_genes, fname=None)
    else:
        write_reports_for_channel(mols, mol_active_genes, report, aggregate_0, detail_0, write_files, to_screen)

    #this part will extract the similar ones
    #this is the eXtendeded search, thus the 'x'
    #First get the similars for each molecule
    if do_sims == True:
        if not silent:
            print("Get the similars...", end=" ")
        sim_mols, sim_mols_detail = get_all_sims(mols, thres, silent, chembl_server=chembl_server)
        
        # extract the selfs (these were handed down before)
        sim_mols_x = sim_mols - set(mol_ids.values())
        if not silent:
            print("N. of similar molecules:", len(sim_mols_x))

        sim_mols_x = list(sim_mols_x)
        mol_ids_x, mol_acts_x, targets_x = get_targets_info_chembl(sim_mols_x, "chemblids", silent=silent,
                                                                    chembl_server=chembl_server)
        #mol_ids_x should be empty
        #mol_acts_x should have the activities for each molecule
        # targets, the info about the respective targets
        #at this phase we know everything
        
        # Similarity reports are intended as evidence for molecules known to be active.
        # Exact-match reports honour the user's -search option; the similarity layer remains active-only.
        mol_active_genes_x, target_ids_x = make_reports_x(mol_acts_x, targets_x, search="A", ofield=ofield)
        
        if "A" in report and (write_files == True or to_screen == True):
            if write_files == True:
                write_report_aggregated(sim_mols_x, mol_active_genes_x, fname=aggregate_1)
            if to_screen == True:
                write_report_aggregated(sim_mols_x, mol_active_genes_x, fname=None)

        genes_T1_by_input = merge_T1(mols, sim_mols_detail, mol_active_genes_x)
        genes_T01 = merge_T0T1(mols, mol_ids, mol_active_genes, sim_mols_detail, mol_active_genes_x)
        
        if "D" in report and (write_files == True or to_screen == True):
            if write_files == True:
                write_report_detail(mols, mol_active_genes, fname=detail_0)
                write_report_detail(mols, genes_T1_by_input, fname=detail_1)
            if to_screen == True:
                write_report_detail(mols, mol_active_genes, fname=None)
                write_report_detail(mols, genes_T1_by_input, fname=None)

        if join_aggregates == True:
            if "A" in report and (write_files == True or to_screen == True):
                if write_files == True:
                    write_report_aggregated(mols, genes_T01, fname=full_aggregate)
                if to_screen == True:
                    write_report_aggregated(mols, genes_T01, fname=None)
            if "D" in report and (write_files == True or to_screen == True):
                if write_files == True:
                    write_report_detail(mols, genes_T01, fname=detail_A)
                if to_screen == True:
                    write_report_detail(mols, genes_T01, fname=None)
    
        return {"mols": mols, "active_genes": mol_active_genes, 
                "mols_x": sim_mols_x, "active_genes_x": mol_active_genes_x, "target_ids": target_ids,
                "target_ids_x": target_ids_x}
    else:
        if join_aggregates == True:
            warn("-joint_aggs was requested but -sim was not used; no joined similarity report created", silent)
        
        return {"mols": mols, "active_genes": mol_active_genes,
                "mols_x": None, "active_genes_x": None, "target_ids": target_ids,
                "target_ids_x": None}
    

    
#def write_report_aggregated(mols, mol_active_genes, fname=None):    
#def write_report_detail(mols, mol_active_genes, fname=None, append=False):



def print_help():
    s = """
Usage: This is a Python tool and requires an enviroment where RDkit and requests are installed.
       To run type: python csanno.py  -in [input .sar file] [options] OR
                    python csanno.py  -mol [molecule SMILES] [options]
Control parameters:
    -in file_name - the data set to annotate (required) (.sar format)
    -sim thr - similarity threshold [if absent (default) no similarity search will be performed]
    -mol [molecule in SMILES format] - molecule for which analysis is going to be performed.  
         Incompatible with -in option
    -report [AD]: A - aggregate report;
                  D - Detailed report (one line for each molecule)
                  Both options simultaneously are permitted and both reports are produced
    -search [ANU] A - searches for actives in the database (default);
                  N - searches non-actives in exact ChEMBL matches; 
                  U - searches unknowns in exact ChEMBL matches;
                  Similarity reports are still active-only evidence
                  All options can be included simultaneously 
    -ofield [GUO] G - outputs the unique gene (default);
                  U - outputs the uniprot id;
                  O - outputs the organism;
                  All options can be included simultaneously
    -chembl_server URL - ChEMBL-compatible server [default: https://www.ebi.ac.uk]

Output Control Options:
    -joint_aggs - with -sim, produces joined T0+T1 aggregate and detail reports
    -silent - no intermediate output at all
    -help - this screen
    -nofiles - no files are created
    -toscreen - writes output to screen
    -pickle - writes a Python pickle for easy posterior analysis (see internal documentation)
"""
    print(version)
    print(s.strip())
    

version = "CSANNO - (C) 2019/2026 - Andre O. Falcao DI/FCUL version 0.4.20260510"
if __name__ == "__main__":
    import sys
    import time
    start_t = time.time()
    if len(sys.argv) < 2:
        print("Invalid option. Exiting")
        exit()

    #defaults
    report = "A"  #default is only the aggregates
    search = "A"  #default only actives 
    ofield = "G"    #default, only outputs the Gene
    do_sims = False
    joint_aggs = False
    sim_thr = 0.9
    writefiles = True
    silent = False
    writepickle = False
    to_screen = False
    mol = None
    input_file = None
    chembl_server = DEFAULT_CHEMBL_SERVER
    
    bin_args = ["-in", "-sim", "-search", "-report", "-mol", "-ofield", "-chembl_server"]
    all_args = bin_args[:] + ["-silent", "-help", "--help", "-nofiles", "-pickle", "-joint_aggs", "-toscreen"]

    inp_fs = sys.argv

    if '-help' in inp_fs or '--help' in inp_fs:
        print_help()
        exit()

    #error checker
    for arg in inp_fs[1:]:
        if len(arg) > 0 and arg[0] == "-" and arg not in all_args:
            print("'%s' Unknown Option. Exiting" % arg)
            print("Run 'csanno.py --help'")
            exit()
	    
    try:
        #arg_screener
        input_args = {}
        for arg in bin_args:
            if arg in inp_fs:
                pos_val = inp_fs.index(arg) + 1
                if pos_val >= len(inp_fs) or inp_fs[pos_val] in all_args:
                    raise ValueError("Missing value for %s" % arg)
                input_args[arg] = inp_fs[pos_val]

        if '-in' in input_args:
            input_file = input_args["-in"]
        if '-sim' in input_args:
            sim_thr = float(input_args["-sim"])
            do_sims = True
        if '-mol' in input_args:
            mol = input_args["-mol"]
            writefiles = False
            to_screen = True
            input_file = None  #if -mol we will ignore any input file
        if '-search' in input_args:
            search = input_args["-search"]
        if '-report' in input_args:
            report = input_args["-report"]
        if '-ofield' in input_args:
            ofield = input_args["-ofield"]
        if '-chembl_server' in input_args:
            chembl_server = normalise_chembl_server(input_args["-chembl_server"])
        
        if '-nofiles' in inp_fs:
            writefiles = False
        if '-toscreen' in inp_fs:
            to_screen = True
        if '-silent' in inp_fs:
            silent = True
        if '-pickle' in inp_fs:
            writepickle = True
        if '-joint_aggs' in inp_fs:
            joint_aggs = True
        
    except Exception as exc:
        print("Illegal Option or Value. Exiting")
        print(exc)
        print("Run 'csanno.py --help'")
        exit()

    if input_file is None and mol is None: 
        print("No input file! Nothing to do. Exiting")
        exit()

    if len(set(search) - set("AUN")) > 0 or len(search) == 0 or len(search) > 3:
        print("Invalid search option: %s  - Exiting" % search)
        exit()

    if len(set(ofield) - set("OGU")) > 0 or len(ofield) == 0 or len(ofield) > 3:
        print("Invalid output field option: %s  - Exiting" % ofield)
        exit()

    if report not in ["A", "D", "AD"]:
        print("Invalid report option: %s - Exiting" % report)
        exit()

    if sim_thr <= 0.0 or sim_thr > 1.0:
        print("Invalid similarity threshold: %s - Exiting" % sim_thr)
        exit()

    if writefiles == False and to_screen == False: 
        print("No files and no screen output. Nothing to see! Exiting")
        exit()

    try:
        if not silent:
            print(version)
        rep = get_activity_profile(fname=input_file, do_sims=do_sims, thres=sim_thr, report=report,
                                   search=search, join_aggregates=joint_aggs, write_files=writefiles,
                                   to_screen=to_screen, silent=silent, single_mol=mol, ofield=ofield,
                                   chembl_server=chembl_server)

        if writepickle == True:
            if input_file is not None:
                fname = input_file.replace(".sar", "") + "_annotations.pickle"
            else:
                fname = "csanno_single_mol_annotations.pickle"
            pickle.dump(rep, open(fname, "wb"))
    except ValueError as exc:
        print("Input error: %s" % exc)
        exit()
    except KeyboardInterrupt:
        print("Interrupted by user")
        exit()
    #print("%9.4f secs" % (time.time() - start_t))
