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
import time

import json
import copy
import requests

try:
    import yaml
except ImportError:
    yaml = None

from csanno_reports import (
    build_run_metadata, normalise_report_format, write_combined_report,
    make_reports, make_reports_x, merge_T0T1, merge_T1,
    write_report_aggregated, write_report_detail, write_reports_for_channel, write_run_info,
)

DEFAULT_CHEMBL_SERVER = "https://www.ebi.ac.uk"
REQUEST_TIMEOUT = 30
REQUEST_RETRIES = 2
REQUEST_RETRY_DELAY_SECONDS = 0.75
RETRIABLE_HTTP_STATUS_CODES = set([429, 500, 502, 503, 504])
TARGET_CACHE_JSON = "./btargets_cache.json"
TARGET_CACHE_PICKLE = "./btargets.pickle"
DEFAULT_ACTIVITY_RULES_FILE = "csanno_default_rules.yaml"

# Collected during one CSANNO run. These are reported once, instead of printing
# one warning line for every transient ChEMBL failure.
REQUEST_FAILURES = []

DEFAULT_ACTIVITY_RULES = {
    "activity_rules_version": "csanno-default-0.5",
    "profile_name": "CSANNO default activity rules",
    "profile_description": "Historical CSANNO thresholds with corrected pKi/pIC50 direction.",
    "default_result": "NA",
    "missing_relation_default": "=",
    "allowed_results": ["A", "N", "NA"],
    "unit_multipliers": {"nM": 1, "uM": 1000, "µM": 1000, "mM": 1000000},
    "rules": [
        {
            "name": "low_concentration_activity",
            "standard_types": ["Ki", "IC50", "Kd", "EC50", "AC50", "Potency"],
            "value_field": "value",
            "units_field": "units",
            "relation_field": "rel",
            "unit_conversion": True,
            "converted_unit": "nM",
            "missing_value_result": "N",
            "missing_units_result": "NA",
            "unsupported_relation_result": "NA",
            "relations": {
                "=":  {"operator": "<", "threshold": 10000, "true_result": "A",  "false_result": "N"},
                ">":  {"operator": "<", "threshold": 10000, "true_result": "NA", "false_result": "N"},
                ">=": {"operator": "<", "threshold": 10000, "true_result": "NA", "false_result": "N"},
                "<":  {"operator": "<", "threshold": 10000, "true_result": "A",  "false_result": "NA"},
                "<=": {"operator": "<", "threshold": 10000, "true_result": "A",  "false_result": "NA"}
            }
        },
        {
            "name": "inhibition_positive_response",
            "standard_types": ["INH", "Inhibition"],
            "value_field": "value",
            "relation_field": "rel",
            "missing_value_result": "N",
            "unsupported_relation_result": "NA",
            "relations": {
                "=":  {"operator": ">", "threshold": 0, "true_result": "A",  "false_result": "N"},
                ">":  {"operator": ">", "threshold": 0, "true_result": "A",  "false_result": "N"},
                ">=": {"operator": ">", "threshold": 0, "true_result": "A",  "false_result": "N"},
                "<":  {"operator": ">", "threshold": 0, "true_result": "NA", "false_result": "N"},
                "<=": {"operator": ">", "threshold": 0, "true_result": "NA", "false_result": "N"}
            }
        },
        {
            "name": "p_activity_high_value_active",
            "standard_types": ["pKi", "pIC50"],
            "value_field": "value",
            "relation_field": "rel",
            "missing_value_result": "N",
            "unsupported_relation_result": "NA",
            "relations": {
                "=":  {"operator": ">", "threshold": 5.0, "true_result": "A",  "false_result": "N"},
                ">":  {"operator": ">", "threshold": 5.0, "true_result": "A",  "false_result": "NA"},
                ">=": {"operator": ">", "threshold": 5.0, "true_result": "A",  "false_result": "NA"},
                "<":  {"operator": ">", "threshold": 5.0, "true_result": "NA", "false_result": "N"},
                "<=": {"operator": ">", "threshold": 5.0, "true_result": "NA", "false_result": "N"}
            }
        },
        {
            "name": "generic_activity_positive_response",
            "standard_types": ["Activity"],
            "value_field": "value",
            "units_field": "units",
            "relation_field": "rel",
            "units_allowed": ["&"],
            "result_if_unit_mismatch": "NA",
            "missing_value_result": "N",
            "unsupported_relation_result": "NA",
            "relations": {
                "=": {"operator": ">", "threshold": 0, "true_result": "A", "false_result": "N"}
            }
        }
    ]
}


def warn(msg, silent=False):
    if not silent:
        print("WARNING: %s" % msg)


def reset_request_failures():
    del REQUEST_FAILURES[:]


def _clean_request_params(params):
    if params is None:
        return {}
    try:
        return {str(key): str(value) for key, value in params.items()}
    except AttributeError:
        return {"params": str(params)}


def record_request_failure(url, params=None, status_code=None, error=None):
    REQUEST_FAILURES.append({
        "url": str(url),
        "params": _clean_request_params(params),
        "status_code": status_code,
        "error": str(error) if error is not None else "",
    })


def get_request_failure_summary(example_limit=8):
    by_status = {}
    by_endpoint = {}
    for item in REQUEST_FAILURES:
        status = str(item.get("status_code") or "network/non-http")
        by_status[status] = by_status.get(status, 0) + 1
        endpoint = item.get("url", "")
        by_endpoint[endpoint] = by_endpoint.get(endpoint, 0) + 1
    examples = REQUEST_FAILURES[:example_limit]
    return {
        "count": len(REQUEST_FAILURES),
        "by_status": by_status,
        "by_endpoint": by_endpoint,
        "examples": examples,
    }


def print_request_failure_summary(silent=False):
    summary = get_request_failure_summary()
    if summary.get("count", 0) == 0 or silent:
        return
    pieces = []
    for status, count in sorted(summary.get("by_status", {}).items()):
        pieces.append("%s: %d" % (status, count))
    warn("Some ChEMBL/HTTP requests failed after retry (%d total; %s). "
         "The run continued, but affected molecules were treated as having no retrievable data. "
         "See the JSON metadata request_warnings block for examples."
         % (summary.get("count", 0), ", ".join(pieces)), silent)


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




def normalise_relation(rel, activity_rules=None):
    if activity_rules is None:
        activity_rules = DEFAULT_ACTIVITY_RULES
    if rel is None or str(rel).strip() == "":
        return activity_rules.get("missing_relation_default", "=")
    rel = str(rel).strip()
    if rel in ["'='", '"="']:
        return "="
    return rel


def copy_default_activity_rules():
    return copy.deepcopy(DEFAULT_ACTIVITY_RULES)


def load_activity_rules(rules_file=None, silent=False):
    """
    Loads activity classification rules from YAML.
    If no rule file is supplied, CSANNO first looks for csanno_default_rules.yaml
    beside this script and in the current directory. If no file is found, the
    internal default rules are used.
    """
    candidate_files = []
    if rules_file is not None and str(rules_file).strip() != "":
        candidate_files.append(str(rules_file).strip())
    else:
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            candidate_files.append(os.path.join(script_dir, DEFAULT_ACTIVITY_RULES_FILE))
        except NameError:
            pass
        candidate_files.append(DEFAULT_ACTIVITY_RULES_FILE)

    rules = None
    source = "internal-default"
    found_file = None
    for candidate in candidate_files:
        if candidate is not None and os.path.isfile(candidate):
            found_file = candidate
            break

    if found_file is not None:
        if yaml is None:
            raise ValueError("PyYAML is required to read activity rule files. Install pyyaml or omit -rules.")
        try:
            with open(found_file, "rt") as fil:
                rules = yaml.safe_load(fil)
            source = found_file
        except Exception as exc:
            raise ValueError("Could not read activity rules file %s: %s" % (found_file, exc))
    else:
        if rules_file is not None and str(rules_file).strip() != "":
            raise ValueError("Activity rules file not found: %s" % rules_file)
        rules = copy_default_activity_rules()

    validate_activity_rules(rules)
    rules["_source"] = source
    if not silent:
        print("Activity rules:", rules.get("activity_rules_version", "unknown"), "[%s]" % source)
    return rules


def validate_activity_rules(activity_rules):
    if not isinstance(activity_rules, dict):
        raise ValueError("Activity rules must be a YAML mapping/dictionary")
    if "rules" not in activity_rules or not isinstance(activity_rules.get("rules"), list):
        raise ValueError("Activity rules must contain a 'rules' list")
    if "default_result" not in activity_rules:
        activity_rules["default_result"] = "NA"
    if "allowed_results" not in activity_rules:
        activity_rules["allowed_results"] = ["A", "N", "NA"]

    allowed = set(activity_rules.get("allowed_results", ["A", "N", "NA"]))
    if activity_rules.get("default_result") not in allowed:
        raise ValueError("default_result must be one of: %s" % ", ".join(sorted(allowed)))

    for i, rule in enumerate(activity_rules.get("rules", []), 1):
        if not isinstance(rule, dict):
            raise ValueError("Activity rule %d is not a mapping/dictionary" % i)
        if "standard_types" not in rule or not isinstance(rule.get("standard_types"), list):
            raise ValueError("Activity rule %d must contain a standard_types list" % i)
        if "relations" not in rule or not isinstance(rule.get("relations"), dict):
            raise ValueError("Activity rule %d must contain a relations mapping" % i)
        for rel in rule.get("relations", {}):
            rel_rule = rule.get("relations", {}).get(rel)
            if not isinstance(rel_rule, dict):
                raise ValueError("Relation rule '%s' in activity rule %d is not a mapping" % (rel, i))
            for result_key in ["true_result", "false_result"]:
                if rel_rule.get(result_key) not in allowed:
                    raise ValueError("%s in rule %d relation %s must be one of: %s" %
                                     (result_key, i, rel, ", ".join(sorted(allowed))))


def compare_rule_value(value, operator, threshold):
    if operator == "<":
        return value < threshold
    if operator == "<=":
        return value <= threshold
    if operator == ">":
        return value > threshold
    if operator == ">=":
        return value >= threshold
    if operator in ["=", "=="]:
        return value == threshold
    if operator in ["!=", "<>"]:
        return value != threshold
    raise ValueError("Unsupported activity-rule operator: %s" % operator)


def apply_activity_rule(act, rule, activity_rules):
    value_field = rule.get("value_field", "value")
    relation_field = rule.get("relation_field", "rel")
    units_field = rule.get("units_field", "units")

    v = safe_float(act.get(value_field))
    if v is None:
        return rule.get("missing_value_result", activity_rules.get("default_result", "NA"))

    unit = act.get(units_field)
    if unit is not None:
        unit = str(unit).strip()

    if rule.get("unit_conversion", False) == True:
        multipliers = activity_rules.get("unit_multipliers", {})
        if unit not in multipliers:
            return rule.get("missing_units_result", activity_rules.get("default_result", "NA"))
        v = v * safe_float(multipliers.get(unit))

    if "units_allowed" in rule:
        units_allowed = rule.get("units_allowed", [])
        if unit not in units_allowed:
            return rule.get("result_if_unit_mismatch", activity_rules.get("default_result", "NA"))

    rel = normalise_relation(act.get(relation_field), activity_rules=activity_rules)
    relations = rule.get("relations", {})
    if rel not in relations:
        return rule.get("unsupported_relation_result", activity_rules.get("default_result", "NA"))

    rel_rule = relations[rel]
    threshold = safe_float(rel_rule.get("threshold"))
    if threshold is None:
        return rule.get("unsupported_relation_result", activity_rules.get("default_result", "NA"))

    try:
        is_true = compare_rule_value(v, rel_rule.get("operator"), threshold)
    except ValueError:
        return rule.get("unsupported_relation_result", activity_rules.get("default_result", "NA"))

    if is_true:
        return rel_rule.get("true_result", activity_rules.get("default_result", "NA"))
    else:
        return rel_rule.get("false_result", activity_rules.get("default_result", "NA"))

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


def read_url(url, params=None, timeout=REQUEST_TIMEOUT, silent=True, retries=REQUEST_RETRIES,
             retry_delay=REQUEST_RETRY_DELAY_SECONDS, warn_immediately=False):
    """Read a JSON URL with small retry/backoff handling for public web-service glitches."""
    last_error = None
    last_status = None
    for attempt in range(int(retries) + 1):
        try:
            response = requests.get(url, params=params, timeout=timeout)
            last_status = response.status_code
            response.raise_for_status()
            return response.json()
        except requests.exceptions.Timeout as exc:
            last_error = exc
            last_status = None
        except requests.exceptions.HTTPError as exc:
            last_error = exc
            response = getattr(exc, "response", None)
            last_status = getattr(response, "status_code", last_status)
        except requests.exceptions.RequestException as exc:
            last_error = exc
            response = getattr(exc, "response", None)
            last_status = getattr(response, "status_code", last_status)
        except ValueError as exc:
            # HTTP succeeded, but the payload was not valid JSON. Retrying is unlikely
            # to help unless it was a transient proxy/service page.
            last_error = exc
            last_status = getattr(locals().get("response", None), "status_code", last_status)

        is_retriable = (last_status in RETRIABLE_HTTP_STATUS_CODES) or isinstance(last_error, requests.exceptions.Timeout)
        if is_retriable and attempt < int(retries):
            time.sleep(float(retry_delay) * (2 ** attempt))
            continue

        record_request_failure(url, params=params, status_code=last_status, error=last_error)
        if warn_immediately:
            if isinstance(last_error, requests.exceptions.Timeout):
                warn("Request timed out after retry: %s" % url, silent)
            elif isinstance(last_error, ValueError):
                warn("Could not decode JSON response from %s" % url, silent)
            else:
                warn("Request failed after retry for %s: %s" % (url, last_error), silent)
        return None

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


def _normalise_chembl_next_url(chembl_server, next_ref):
    if next_ref is None or str(next_ref).strip() == "":
        return None
    next_ref = str(next_ref)
    if next_ref.startswith("http://") or next_ref.startswith("https://"):
        return next_ref
    if next_ref.startswith("/"):
        return normalise_chembl_server(chembl_server) + next_ref
    return normalise_chembl_server(chembl_server) + "/" + next_ref


def _extract_activity_records(data, mol_id, silent=False):
    activities = []
    for act in (data or {}).get("activities", []):
        try:
            record = {"assay_id": act.get("assay_chembl_id"),
                      "target_id": act.get('target_chembl_id'),
                      "assay_type": act.get('standard_type'),
                      "units": act.get('standard_units'),
                      "value": act.get('standard_value'),
                      "rel": act.get('standard_relation'),
                      "pvalue": act.get('pchembl_value')}
            if record["target_id"] is not None:
                activities.append(record)
        except (AttributeError, KeyError) as exc:
            warn("Skipping malformed activity for %s: %s" % (mol_id, exc), silent)
    return activities


def get_activities_chembl(mol_id, chembl_server=DEFAULT_CHEMBL_SERVER, silent=False):
    """For a given ChEMBL molecule id, return activity records.

    ChEMBL activity responses are paginated. Using a normal page size and
    following page_meta['next'] is gentler on the public service than asking
    for a special all-at-once response.
    """
    if mol_id is None or mol_id is False:
        return []

    chembl_server = normalise_chembl_server(chembl_server)
    url = chembl_server + "/chembl/api/data/activity.json"
    params = {"molecule_chembl_id": mol_id, "format": "json", "limit": 1000}

    activities = []
    while url is not None:
        data = read_url(url, params=params, silent=silent)
        # Only the first request uses params; ChEMBL's page_meta.next already
        # contains all filters and offsets for subsequent pages.
        params = None
        if data is None:
            break
        activities.extend(_extract_activity_records(data, mol_id, silent=silent))
        next_ref = (data.get("page_meta", {}) or {}).get("next")
        url = _normalise_chembl_next_url(chembl_server, next_ref)

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


def check_activity(act, activity_rules=None):
    if activity_rules is None:
        activity_rules = DEFAULT_ACTIVITY_RULES

    assay_type = act.get("assay_type")
    if assay_type is None:
        return activity_rules.get("default_result", "NA")

    for rule in activity_rules.get("rules", []):
        if assay_type in rule.get("standard_types", []):
            return apply_activity_rule(act, rule, activity_rules)

    return activity_rules.get("default_result", "NA")
    
    



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




def get_activity_profile(fname=None, do_sims=False, thres=0.7, report="A", write_files=True, to_screen=False,
                         join_aggregates=False, silent=False, search="A", single_mol=None, ofield="G",
                         chembl_server=DEFAULT_CHEMBL_SERVER, rules_file=None, activity_rules=None,
                         legacy_reports=False, report_format="markdown", include_glossary=True, include_digest=True,
                         out_file=None):
    reset_request_failures()
    """Main entry point.

    The default file output is now:
       <root>_report.md     one readable Markdown report containing metadata, digest, glossary, and raw data
       <root>_results.json  the same metadata/digest plus machine-readable raw results

    Set report_format="text" to write <root>_report.txt instead.

    Set legacy_reports=True from Python, or use -legacy_reports from the command
    line, to also produce the old split aggregate/detail files.
    """
    if activity_rules is None:
        activity_rules = load_activity_rules(rules_file=rules_file, silent=silent)

    report_format, report_ext = normalise_report_format(report_format)

    # Output names. The legacy names are still created only when legacy_reports=True.
    if write_files == True:
        #if fname is not None:
        #    root = os.path.splitext(fname)[0]
        #else:
        #    root = "csanno_single_mol"
        root = out_file    
        report_file = root + "_report." + report_ext
        report_json = root + "_results.json"
        detail_0 = root + "_annotations_T0.txt"
        detail_1 = root + "_annotations_T1.txt"
        detail_A = root + "_annotations_TA.txt"
        aggregate_0 = root + "_anno_count_T0.txt"
        aggregate_1 = root + "_anno_count_T1.txt"
        full_aggregate = root + "_anno_count_full.txt"
        run_info_file = root + "_run_info.txt"
    else:
        report_file = None
        report_json = None
        detail_0 = None
        detail_1 = None
        detail_A = None
        aggregate_0 = None
        aggregate_1 = None
        full_aggregate = None
        run_info_file = None

    if single_mol is not None:
        mols = make_mol(single_mol)
    else:
        if fname is not None:
            mols = read_molecules(fname, silent=silent)
        else:
            print("No input file! Nothing to do. Exiting")
            exit()

    metadata = build_run_metadata(csanno_version=version,
                                  input_file=fname,
                                  single_mol=single_mol,
                                  mols=mols,
                                  chembl_server=chembl_server,
                                  do_sims=do_sims,
                                  thres=thres,
                                  search=search,
                                  report=report,
                                  ofield=ofield,
                                  join_aggregates=join_aggregates,
                                  activity_rules=activity_rules,
                                  report_format=report_format)

    if write_files == True and legacy_reports == True:
        write_run_info(run_info_file,
                       {"csanno_version": version,
                        "activity_rules_version": activity_rules.get("activity_rules_version", "unknown"),
                        "activity_rules_source": activity_rules.get("_source", "unknown"),
                        "chembl_server": chembl_server,
                        "input_file": fname,
                        "single_molecule_mode": single_mol is not None,
                        "molecule_count": len(mols),
                        "similarity_search": do_sims,
                        "similarity_threshold": thres,
                        "search": search,
                        "report": report,
                        "ofield": ofield,
                        "join_aggregates": join_aggregates,
                        "readable_report": report_file,
                        "json_report": report_json},
                       silent=silent)

    mids = list(mols.keys())
    icks = []
    for mid in mids:
        icks.append(mols[mid][1])

    # Tier 0: exact matches in ChEMBL.
    mol_ids, mol_acts, targets = get_targets_info_chembl(icks, silent=silent, chembl_server=chembl_server)

    mol_active_genes, target_ids = make_reports(mols, mol_ids, mol_acts, targets, search=search, ofield=ofield,
                                                activity_rules=activity_rules, activity_checker=check_activity)

    if legacy_reports == True:
        if do_sims == True:
            if "A" in report and write_files == True:
                write_report_aggregated(mols, mol_active_genes, fname=aggregate_0)
        else:
            write_reports_for_channel(mols, mol_active_genes, report, aggregate_0, detail_0,
                                      write_files and legacy_reports, False)

    # Tier 1: similar ChEMBL molecules.
    if do_sims == True:
        if not silent:
            print("Get the similars...", end=" ")
        sim_mols, sim_mols_detail = get_all_sims(mols, thres, silent, chembl_server=chembl_server)

        # Exclude exact/self matches already handled by Tier 0.
        sim_mols_x = sim_mols - set(mol_ids.values())
        if not silent:
            print("N. of similar molecules:", len(sim_mols_x))

        sim_mols_x = list(sim_mols_x)
        mol_ids_x, mol_acts_x, targets_x = get_targets_info_chembl(sim_mols_x, "chemblids", silent=silent,
                                                                    chembl_server=chembl_server)

        # Similarity reports honour the same -search option as exact-match reports.
        mol_active_genes_x, target_ids_x = make_reports_x(mol_acts_x, targets_x, search=search, ofield=ofield,
                                                          activity_rules=activity_rules, activity_checker=check_activity)

        if legacy_reports == True and "A" in report and write_files == True:
            write_report_aggregated(sim_mols_x, mol_active_genes_x, fname=aggregate_1)

        genes_T1_by_input = merge_T1(mols, sim_mols_detail, mol_active_genes_x)
        genes_T01 = merge_T0T1(mols, mol_ids, mol_active_genes, sim_mols_detail, mol_active_genes_x)

        if legacy_reports == True and "D" in report and write_files == True:
            write_report_detail(mols, mol_active_genes, fname=detail_0)
            write_report_detail(mols, genes_T1_by_input, fname=detail_1)

        if join_aggregates == True and legacy_reports == True:
            if "A" in report and write_files == True:
                write_report_aggregated(mols, genes_T01, fname=full_aggregate)
            if "D" in report and write_files == True:
                write_report_detail(mols, genes_T01, fname=detail_A)

        rep = {"mols": mols,
               "active_genes": mol_active_genes,
               "mols_x": sim_mols_x,
               "sim_mols_detail": sim_mols_detail,
               "active_genes_x": mol_active_genes_x,
               "active_genes_T1_by_input": genes_T1_by_input,
               "active_genes_T01": genes_T01,
               "target_ids": target_ids,
               "target_ids_x": target_ids_x,
               "targets": targets,
               "targets_x": targets_x}
    else:
        if join_aggregates == True:
            warn("-joint_aggs was requested but -sim was not used; no joined similarity report created", silent)

        rep = {"mols": mols,
               "active_genes": mol_active_genes,
               "mols_x": None,
               "sim_mols_detail": None,
               "active_genes_x": None,
               "active_genes_T1_by_input": None,
               "active_genes_T01": None,
               "target_ids": target_ids,
               "target_ids_x": None,
               "targets": targets,
               "targets_x": None}

    request_summary = get_request_failure_summary()
    if request_summary.get("count", 0) > 0:
        metadata["request_warnings"] = request_summary
        rep["request_warnings"] = request_summary
    else:
        metadata["request_warnings"] = {"count": 0, "by_status": {}, "by_endpoint": {}, "examples": []}
        rep["request_warnings"] = metadata["request_warnings"]

    rep["metadata"] = metadata
    rep["output_files"] = {"readable_report": report_file, "json_report": report_json}

    print_request_failure_summary(silent=silent)

    if write_files == True or to_screen == True:
        write_combined_report(rep, metadata,
                              report_fname=report_file if write_files == True else None,
                              json_fname=report_json if write_files == True else None,
                              report=report,
                              include_sims=do_sims,
                              include_joined=join_aggregates,
                              to_screen=to_screen,
                              silent=silent,
                              report_format=report_format,
                              include_glossary=include_glossary, include_digest=include_digest)

    return rep


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
    -out file_name - output file name [by default will append to the input file name]
    -report [AD]: A - aggregate report;
                  D - Detailed report (one line for each molecule)
                  Both options simultaneously are permitted and both reports are produced
    -search [ANU] A - searches for actives in the database (default);
                  N - searches non-actives in exact ChEMBL matches; 
                  U - searches unknowns in exact ChEMBL matches;
                  Similarity reports use the same search option
                  All options can be included simultaneously 
    -ofield [GUO] G - outputs the unique gene (default);
                  U - outputs the uniprot id;
                  O - outputs the organism;
                  All options can be included simultaneously
    -chembl_server URL - ChEMBL-compatible server [default: https://www.ebi.ac.uk]
    -rules YAML_file - YAML file defining active/non-active/uncertain classification rules

Output Control Options:
    By default, CSANNO writes one readable Markdown report (<root>_report.md) and one
    JSON report (<root>_results.json) containing metadata, digest, glossary, and raw data.
    -outfmt [markdown|text] - writes the readable report as .md (default) or .txt
    -legacy_reports - also writes the old split T0/T1 aggregate/detail text files
    -joint_aggs - with -sim and -legacy_reports, also produces joined T0+T1 legacy files
    -silent - no intermediate output at all
    -help - this screen
    -nofiles - no files are created
    -toscreen - writes the readable report to screen
    -pickle - writes a Python pickle for easy posterior analysis (see internal documentation)
"""
    print(version)
    print(s.strip())
    

version = "CSANNO - (C) 2019/2026 - Andre O. Falcao BioISI - DI/FCUL version 0.5.20260510"
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
    legacy_reports = False
    report_format = "markdown"
    mol = None
    input_file = None
    out_file = None
    chembl_server = DEFAULT_CHEMBL_SERVER
    rules_file = None
    digest=False
    glossary=False
    
    bin_args = ["-in", "-sim", "-search", "-report", "-mol", "-ofield", "-chembl_server", "-rules", "-outfmt", "-out"]
    all_args = bin_args[:] + ["-silent", "-help", "--help", "-nofiles", "-pickle", "-joint_aggs", 
                              "-toscreen", "-legacy_reports", "-legacy", "-digest", "-glossary"]

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
        if '-rules' in input_args:
            rules_file = input_args["-rules"]
        if '-outfmt' in input_args:
            report_format = normalise_report_format(input_args["-outfmt"])[0]
        if '-out' in input_args:
            out_file = input_args["-out"]
            writefiles = True
            to_screen = False
        
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
        if '-legacy_reports' in inp_fs or '-legacy' in inp_fs:
            legacy_reports = True
        if '-digest' in inp_fs: 
            digest = True
        if '-glossary' in inp_fs: 
            glossary = True
        
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


    if out_file is None:
        if input_file is not None: out_file = os.path.splitext(input_file)[0]
        else: out_file = "csanno_single_mol"
            
    try:
        if not silent:
            print(version)

        rep = get_activity_profile(fname=input_file, do_sims=do_sims, thres=sim_thr, report=report,
                                   search=search, join_aggregates=joint_aggs, write_files=writefiles,
                                   to_screen=to_screen, silent=silent, single_mol=mol, ofield=ofield,
                                   chembl_server=chembl_server, rules_file=rules_file,
                                   legacy_reports=legacy_reports, report_format=report_format, 
                                   include_glossary=glossary, include_digest=digest,
                                   out_file=out_file)

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
