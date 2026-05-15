"""ChEMBL data-access layer for CSANNO.

This module is intentionally restricted to repository/data retrieval.  The main
CSANNO module should not know how ChEMBL is queried; it should only consume the
normalised dictionaries returned here.  A future local-database or alternative
repository backend can replace this file if it implements the same small API:

    get_targets_info(mol_list, data_type="inchikeys", silent=False, server=...)
    get_all_similar_molecules(mols, thres, silent, server=...)
    reset_request_failures()
    get_request_failure_summary()
    print_request_failure_summary()
    normalise_server(server)

The activity records returned by get_activities() are normalised to the schema
expected by CSANNO's rule engine:
    assay_id, target_id, assay_type, units, value, rel, pvalue
"""

import json
import os
import pickle
import time

import requests


DEFAULT_SERVER = "https://www.ebi.ac.uk"
DEFAULT_CHEMBL_SERVER = DEFAULT_SERVER  # backwards-compatible alias

REQUEST_TIMEOUT = 30
REQUEST_RETRIES = 2
REQUEST_RETRY_DELAY_SECONDS = 0.75
RETRIABLE_HTTP_STATUS_CODES = set([429, 500, 502, 503, 504])

TARGET_CACHE_JSON = "./btargets_cache.json"
TARGET_CACHE_PICKLE = "./btargets.pickle"

# Collected during one CSANNO run. These are reported once, instead of printing
# one warning line for every transient ChEMBL failure.
REQUEST_FAILURES = []


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


def normalise_server(server):
    if server is None or str(server).strip() == "":
        return DEFAULT_SERVER
    return str(server).rstrip("/")


def normalise_chembl_server(chembl_server):
    """Backward-compatible name used by older CSANNO code."""
    return normalise_server(chembl_server)


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


def _normalise_next_url(server, next_ref):
    if next_ref is None or str(next_ref).strip() == "":
        return None
    next_ref = str(next_ref)
    if next_ref.startswith("http://") or next_ref.startswith("https://"):
        return next_ref
    if next_ref.startswith("/"):
        return normalise_server(server) + next_ref
    return normalise_server(server) + "/" + next_ref


def get_molecule_id(inchikey, server=DEFAULT_SERVER, silent=False):
    """Return the ChEMBL molecule id for an InChIKey, or False when absent."""
    server = normalise_server(server)
    url = server + "/chembl/api/data/molecule.json"
    params = {"molecule_structures__standard_inchi_key": inchikey, "format": "json"}
    data = read_url(url, params=params, silent=silent)
    if data is None:
        return False

    molecules = data.get("molecules", [])
    if len(molecules) == 1:
        return molecules[0].get("molecule_chembl_id", False)
    return False


def _extract_activity_records(data, mol_id, silent=False):
    activities = []
    for act in (data or {}).get("activities", []):
        try:
            record = {"assay_id": act.get("assay_chembl_id"),
                      "target_id": act.get("target_chembl_id"),
                      "assay_type": act.get("standard_type"),
                      "units": act.get("standard_units"),
                      "value": act.get("standard_value"),
                      "rel": act.get("standard_relation"),
                      "pvalue": act.get("pchembl_value")}
            if record["target_id"] is not None:
                activities.append(record)
        except (AttributeError, KeyError) as exc:
            warn("Skipping malformed activity for %s: %s" % (mol_id, exc), silent)
    return activities


def get_activities(molecule_id, server=DEFAULT_SERVER, silent=False):
    """Return normalised activity records for a repository molecule id.

    ChEMBL activity responses are paginated. Using a normal page size and
    following page_meta['next'] is gentler on the public service than asking
    for a special all-at-once response.
    """
    if molecule_id is None or molecule_id is False:
        return []

    server = normalise_server(server)
    url = server + "/chembl/api/data/activity.json"
    params = {"molecule_chembl_id": molecule_id, "format": "json", "limit": 1000}

    activities = []
    while url is not None:
        data = read_url(url, params=params, silent=silent)
        # Only the first request uses params; ChEMBL's page_meta.next already
        # contains all filters and offsets for subsequent pages.
        params = None
        if data is None:
            break
        activities.extend(_extract_activity_records(data, molecule_id, silent=silent))
        next_ref = (data.get("page_meta", {}) or {}).get("next")
        url = _normalise_next_url(server, next_ref)

    return activities


def get_target_components(target_id, server=DEFAULT_SERVER, silent=False):
    """Return biological target-component information for a ChEMBL target id."""
    if target_id is None or target_id is False:
        return []

    server = normalise_server(server)
    url = server + "/chembl/api/data/target.json"
    params = {"target_chembl_id": target_id, "format": "json"}
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


def get_target_info(inchikey, server=DEFAULT_SERVER, silent=False):
    molecule_id = get_molecule_id(inchikey, server=server, silent=silent)
    if not molecule_id:
        return {}
    acts = get_activities(molecule_id, server=server, silent=silent)
    targets = {}
    for act in acts:
        targs = get_target_components(act.get("target_id"), server=server, silent=silent)
        if len(targs) >= 1:
            targets[act["target_id"]] = targs
    return targets


def get_targets_info(mol_list, data_type="inchikeys", silent=False, server=DEFAULT_SERVER):
    """
    For each molecule in mol_list, retrieve activity records and target metadata.

    Returns:
       mol_ids: dictionary with inchikey as key and database id as value when
                data_type == "inchikeys"; empty for data_type == "chemblids".
       mol_acts: dictionary with database molecule id as key and activity list as value.
       btargets: dictionary with target id as key and target-component tuples
                 (organism, gene, uniprot) as value.
    """
    mol_ids = {}
    mol_acts = {}
    if not silent:
        print("Processing molecules...")

    if data_type == "inchikeys":
        i = 0
        for ick in mol_list:
            molecule_id = get_molecule_id(ick, server=server, silent=silent)
            if molecule_id:
                mol_ids[ick] = molecule_id
                acts = get_activities(molecule_id, server=server, silent=silent)
                mol_acts[molecule_id] = acts
                i += 1
                if i % 100 == 0 and not silent:
                    print("\tProcessed %4d molecules" % i)
    elif data_type == "chemblids":
        i = 0
        for molecule_id in mol_list:
            acts = get_activities(molecule_id, server=server, silent=silent)
            mol_acts[molecule_id] = acts
            i += 1
            if i % 100 == 0 and not silent:
                print("\tProcessed %4d molecules" % i)
    else:
        print("Wrong data type. Exiting")
        return {}, {}, {}

    ctargets = set()
    for mid in mol_acts:
        acts = mol_acts[mid]
        for act in acts:
            target_id = act.get("target_id")
            if target_id is not None:
                ctargets.add(target_id)

    if not silent:
        print("Getting the Chembl target data...")

    btargets = load_target_cache(silent=silent)
    update_cache = False
    i = 0
    for target_id in ctargets:
        if target_id in btargets:
            continue

        update_cache = True
        if not silent:
            print("\tthis one is new! ->", target_id)
        btargs = get_target_components(target_id, server=server, silent=silent)
        btargets[target_id] = []
        for bt in btargs:
            btargets[target_id].append((bt.get("organism", ""), bt.get("gene", ""), bt.get("uniprot", "")))
        i += 1
        if i % 100 == 0 and not silent:
            print("\tProcessed %4d targets" % i)

    if update_cache == True:
        save_target_cache(btargets, silent=silent)
    if not silent:
        print("Done!")

    return mol_ids, mol_acts, btargets


def get_similar_molecules(smiles, sim_thr=0.9, server=DEFAULT_SERVER, silent=False):
    """Return repository molecule ids similar to the supplied SMILES."""
    sim_mols = []
    ssim = str(int(sim_thr * 100))
    server = normalise_server(server)

    url = server + "/chembl/api/data/similarity.json"
    params = {"smiles": smiles, "limit": 0, "similarity": ssim}
    data = read_url(url, params=params, silent=silent)
    if data is None:
        return []

    the_mols = data.get("molecules", [])
    page_meta = data.get("page_meta", {})
    while page_meta.get("next") is not None:
        next_url = _normalise_next_url(server, page_meta.get("next"))
        data = read_url(next_url, silent=silent)
        if data is None:
            break
        the_mols += data.get("molecules", [])
        page_meta = data.get("page_meta", {})

    for mol in the_mols:
        molecule_id = mol.get("molecule_chembl_id")
        if molecule_id is not None:
            sim_mols.append(molecule_id)
    return sim_mols


def get_all_similar_molecules(mols, thres, silent, server=DEFAULT_SERVER):
    sim_mols = []
    sim_mols_detail = {}
    i = 0
    for mid in mols:
        smi = mols[mid][0]
        sims = get_similar_molecules(smi, thres, server=server, silent=silent)
        sim_mols += sims
        sim_mols_detail[mid] = sims
        i += 1
        if i % 100 == 0 and not silent:
            print("\tProcessed %4d similarities" % i)
    return set(sim_mols), sim_mols_detail


# ---------------------------------------------------------------------------
# Backwards-compatible wrappers for older scripts that imported these names.


def get_chembl_id(inchikey, chembl_server=DEFAULT_SERVER, silent=False):
    return get_molecule_id(inchikey, server=chembl_server, silent=silent)


def get_activities_chembl(mol_id, chembl_server=DEFAULT_SERVER, silent=False):
    return get_activities(mol_id, server=chembl_server, silent=silent)


def get_targets_chembl(tid, chembl_server=DEFAULT_SERVER, silent=False):
    return get_target_components(tid, server=chembl_server, silent=silent)


def get_target_info_chembl(inchikey, chembl_server=DEFAULT_SERVER, silent=False):
    return get_target_info(inchikey, server=chembl_server, silent=silent)


def get_targets_info_chembl(mol_list, data_type="inchikeys", silent=False, chembl_server=DEFAULT_SERVER):
    return get_targets_info(mol_list, data_type=data_type, silent=silent, server=chembl_server)


def get_sims(smiles, sim_thr=0.9, chembl_server=DEFAULT_SERVER, silent=False):
    return get_similar_molecules(smiles, sim_thr=sim_thr, server=chembl_server, silent=silent)


def get_all_sims(mols, thres, silent, chembl_server=DEFAULT_SERVER):
    return get_all_similar_molecules(mols, thres=thres, silent=silent, server=chembl_server)
