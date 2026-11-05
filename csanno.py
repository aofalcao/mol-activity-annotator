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


from rdkit import Chem
from rdkit.Chem.inchi import InchiToInchiKey
from rdkit.Chem.MolStandardize import rdMolStandardize
import argparse
import pickle
import os
import sys
import time
import shlex

import copy

try:
    import yaml
except ImportError:
    yaml = None

from csanno_reports import (
    build_run_metadata, normalise_report_format, write_combined_report,
    make_reports, make_reports_x, merge_T0T1, merge_T1,
    write_report_aggregated, write_report_detail, write_reports_for_channel, write_run_info,
)

import chembl_access as data_access

DEFAULT_CHEMBL_SERVER = data_access.DEFAULT_SERVER
DEFAULT_ACTIVITY_RULES_FILE = "csanno_default_rules.yaml"


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

def parse_molecule(smiles):
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

    return m


def standardize_input_molecule(smiles):
    m = parse_molecule(smiles)
    try:
        m = rdMolStandardize.Cleanup(m)
        m = rdMolStandardize.FragmentParent(m)
        m = rdMolStandardize.Uncharger().uncharge(m)
        m = rdMolStandardize.IsotopeParent(m)
        Chem.SanitizeMol(m)
    except Exception as exc:
        raise ValueError("Could not standardize SMILES '%s': %s" % (smiles, exc))

    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)


def get_chem_data(smiles, standardize=False):
    if standardize == True:
        smiles = standardize_input_molecule(smiles)

    m = parse_molecule(smiles)
    inchi = Chem.inchi.MolToInchi(m)
    inchikey = InchiToInchiKey(inchi)
    if inchikey is None or inchikey == "":
        raise ValueError("Could not generate InChIKey for SMILES '%s'" % smiles)
    return inchikey


def make_mol(smiles, standardize=False):
    mols = {}
    if standardize == True:
        smiles = standardize_input_molecule(smiles)
    ick = get_chem_data(smiles)
    mols["User1"] = (smiles, ick, "NA")
    return mols
                   


def read_molecules(fname, silent=False, standardize=False):
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
                if standardize == True:
                    smiles = standardize_input_molecule(smiles)
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
    
    
def get_activity_profile(fname=None, do_sims=False, thres=0.7, report="A", write_files=True, to_screen=False,
                         join_aggregates=False, silent=False, search="A", single_mol=None, ofield="G",
                         chembl_server=DEFAULT_CHEMBL_SERVER, rules_file=None, activity_rules=None,
                         legacy_reports=False, report_format="markdown", include_glossary=True, include_digest=True,
                         out_file=None, run_command=None, standardize=False):
    """Main entry point.

    The default file output is now:
       <root>_report.md     one readable Markdown report containing metadata, digest, glossary, and raw data
       <root>_results.json  the same metadata/digest plus machine-readable raw results

    Set report_format="text" to write <root>_report.txt instead.

    Set legacy_reports=True from Python, or use -legacy_reports from the command
    line, to also produce the old split aggregate/detail files.
    """
    data_access.reset_request_failures()

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
        mols = make_mol(single_mol, standardize=standardize)
    else:
        if fname is not None:
            mols = read_molecules(fname, silent=silent, standardize=standardize)
        else:
            raise ValueError("No input file or molecule was provided")

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
                                  report_format=report_format,
                                  run_command=run_command, standardize=standardize)

    if write_files == True and legacy_reports == True:
        write_run_info(run_info_file,
                       {"csanno_version": version,
                        "run_command": run_command,
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
                        "json_report": report_json, "standardize": standardize,},
                       silent=silent)

    mids = list(mols.keys())
    icks = []
    for mid in mids:
        icks.append(mols[mid][1])

    # Tier 0: exact matches in ChEMBL.
    mol_ids, mol_acts, targets = data_access.get_targets_info(icks, silent=silent, server=chembl_server)

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
        sim_mols, sim_mols_detail = data_access.get_all_similar_molecules(mols, thres, silent, server=chembl_server)

        # Exclude exact/self matches already handled by Tier 0.
        sim_mols_x = sim_mols - set(mol_ids.values())
        if not silent:
            print("N. of similar molecules:", len(sim_mols_x))

        sim_mols_x = list(sim_mols_x)
        mol_ids_x, mol_acts_x, targets_x = data_access.get_targets_info(sim_mols_x, "chemblids", silent=silent,
                                             server=chembl_server)

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

    request_summary = data_access.get_request_failure_summary()
    if request_summary.get("count", 0) > 0:
        metadata["request_warnings"] = request_summary
        rep["request_warnings"] = request_summary
    else:
        metadata["request_warnings"] = {"count": 0, "by_status": {}, "by_endpoint": {}, "examples": []}
        rep["request_warnings"] = metadata["request_warnings"]

    rep["metadata"] = metadata
    rep["output_files"] = {"readable_report": report_file, "json_report": report_json}

    data_access.print_request_failure_summary(silent=silent)

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



def build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="csanno.py",
        add_help=False,
        description=(
            "Annotate molecules with direct ChEMBL evidence (Tier 0) and, "
            "optionally, evidence from structurally similar ChEMBL molecules (Tier 1)."
        ),
    )

    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("-in", dest="input_file", metavar="FILE",
                        help="input data set in three-column tab-separated .sar format")
    inputs.add_argument("-mol", dest="mol", metavar="SMILES",
                        help="analyse one molecule supplied as a SMILES string")

    parser.add_argument("-sim", dest="sim_thr", type=float, metavar="THRESHOLD",
                        help="enable Tier 1 similarity searching at this threshold")
    parser.add_argument("-search", default="A", metavar="ANU",
                        help="evidence classes to include: A, N and/or U (default: A)")
    parser.add_argument("-report", default="AD", choices=["A", "D", "AD"],
                        help="report sections: A, D or AD (default: AD)")
    parser.add_argument("-ofield", default="G", metavar="GUO",
                        help="target fields to output: G, U and/or O (default: G)")
    parser.add_argument("-chembl_server", default=DEFAULT_CHEMBL_SERVER, metavar="URL",
                        help="ChEMBL-compatible data server")
    parser.add_argument("-rules", dest="rules_file", metavar="YAML_FILE",
                        help="YAML file defining activity classification rules")
    parser.add_argument("-standardize", action="store_true",
                        help="standardize input molecules before Tier 0 and Tier 1 searches")

    parser.add_argument("-out", dest="out_file", metavar="NAME",
                        help="output file root")
    parser.add_argument("-outfmt", default="markdown",
                        choices=["markdown", "md", "text", "txt"],
                        help="readable report format (default: markdown)")
    parser.add_argument("-legacy_reports", "-legacy", dest="legacy_reports",
                        action="store_true",
                        help="also write the historical split report files")
    parser.add_argument("-joint_aggs", action="store_true",
                        help="with -sim and legacy reports, write joined T0+T1 reports")
    parser.add_argument("-glossary", action="store_true",
                        help="include a gene glossary")
    parser.add_argument("-digest", action="store_true",
                        help="include the report digest")
    parser.add_argument("-silent", action="store_true",
                        help="suppress intermediate output")
    parser.add_argument("-nofiles", action="store_true",
                        help="do not write output files")
    parser.add_argument("-toscreen", action="store_true",
                        help="write the readable report to the screen")
    parser.add_argument("-pickle", dest="writepickle", action="store_true",
                        help="write a Python pickle for later analysis")
    parser.add_argument("-help", "--help", action="help",
                        help="show this help message and exit")
    parser.add_argument("--version", action="version", version=version)

    return parser



def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = build_arg_parser()
    args = parser.parse_args(argv)

    run_command = shlex.join([sys.argv[0]] + list(argv))

    report = args.report
    search = args.search
    ofield = args.ofield
    do_sims = args.sim_thr is not None
    joint_aggs = args.joint_aggs
    sim_thr = args.sim_thr if args.sim_thr is not None else 0.9
    writefiles = not args.nofiles
    silent = args.silent
    writepickle = args.writepickle
    to_screen = args.toscreen
    legacy_reports = args.legacy_reports
    report_format = normalise_report_format(args.outfmt)[0]
    mol = args.mol
    input_file = args.input_file
    out_file = args.out_file
    chembl_server = data_access.normalise_server(args.chembl_server)
    rules_file = args.rules_file
    digest = args.digest
    glossary = args.glossary
    standardize = args.standardize

    if mol is not None:
        to_screen = True

    if out_file is not None:
        writefiles = True
        to_screen = False

    if input_file is None and mol is None:
        parser.error("one of -in FILE or -mol SMILES is required")

    if len(set(search) - set("AUN")) > 0 or len(search) == 0 or len(search) > 3:
        parser.error("invalid -search value %s; use A, N and/or U" % search)

    if len(set(ofield) - set("OGU")) > 0 or len(ofield) == 0 or len(ofield) > 3:
        parser.error("invalid -ofield value %s; use O, G and/or U" % ofield)

    if sim_thr <= 0.0 or sim_thr > 1.0:
        parser.error("-sim threshold must be greater than 0 and at most 1")

    if writefiles == False and to_screen == False:
        parser.error("-nofiles requires -toscreen or -mol")

    if out_file is None:
        if input_file is not None:
            out_file = os.path.splitext(input_file)[0]
        else:
            out_file = "csanno_single_mol"

    try:
        if not silent:
            print(version)

        rep = get_activity_profile(
            fname=input_file,
            do_sims=do_sims,
            thres=sim_thr,
            report=report,
            search=search,
            join_aggregates=joint_aggs,
            write_files=writefiles,
            to_screen=to_screen,
            silent=silent,
            single_mol=mol,
            ofield=ofield,
            chembl_server=chembl_server,
            rules_file=rules_file,
            legacy_reports=legacy_reports,
            report_format=report_format,
            include_glossary=glossary,
            include_digest=digest,
            out_file=out_file,
            run_command=run_command,
            standardize=standardize,
        )

        if writepickle == True:
            if input_file is not None:
                fname = input_file.replace(".sar", "") + "_annotations.pickle"
            else:
                fname = "csanno_single_mol_annotations.pickle"

            with open(fname, "wb") as fil:
                pickle.dump(rep, fil)

    except ValueError as exc:
        print("Input error: %s" % exc, file=sys.stderr)
        return 1

    except KeyboardInterrupt:
        print("Interrupted by user", file=sys.stderr)
        return 130

    return 0


version = "CSANNO - (C) 2019/2026 - Andre O. Falcao BioISI - DI/FCUL version 0.5.20260727"


if __name__ == "__main__":
    sys.exit(main())
