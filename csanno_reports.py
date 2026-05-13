"""Readable text, Markdown, legacy text, and JSON reporting for CSANNO.

This module deliberately contains no RDKit or ChEMBL calls. It receives
already-computed CSANNO dictionaries and handles:
  * Tier 0/Tier 1 target extraction from activity records;
  * legacy aggregate/detail text files;
  * one human-readable report as Markdown (.md) or plain text (.txt);
  * one JSON file with metadata, digest, sections, glossary, and raw results.
"""

from collections import Counter
from datetime import datetime
import json
import math
import os
import re

try:
    import requests
except ImportError:  # pragma: no cover - requests is already required by CSANNO proper
    requests = None


SEARCH_ORDER = "ANU"

SEARCH_TEXT = {
    "A": "active / positive evidence",
    "N": "non-active / negative evidence",
    "U": "uncertain or not classifiable evidence",
}

REPORT_TEXT = {
    "A": "aggregate target-frequency tables",
    "D": "molecule-by-molecule detail tables",
}

OFIELD_TEXT = {
    "O": "organism",
    "G": "gene symbol",
    "U": "UniProt accession",
}

REPORT_FORMATS = {
    "markdown": ("markdown", "md"),
    "md": ("markdown", "md"),
    "text": ("text", "txt"),
    "txt": ("text", "txt"),
}

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
GENECARDS_URL = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s"
UNIPROT_ENTRY_URL = "https://www.uniprot.org/uniprotkb/%s"
DEFAULT_GLOSSARY_CACHE_JSON = "./csanno_gene_glossary_cache.json"
GLOSSARY_REQUEST_TIMEOUT = 4

COMMON_NON_GENE_TOKENS = {
    "HOMO", "SAPIENS", "MUS", "MUSCULUS", "RATTUS", "NORVEGICUS", "DANIO", "RERIO",
    "BOS", "TAURUS", "CANIS", "FAMILIARIS", "SUS", "SCROFA", "DROSOPHILA", "MELANOGASTER",
    "ARABIDOPSIS", "THALIANA", "ESCHERICHIA", "COLI", "STAPHYLOCOCCUS", "AUREUS",
    "SARS", "COV", "INFLUENZA", "VIRUS", "ORGANISM", "UNKNOWN", "NA", "NONE",
}


# ---------------------------------------------------------------------------
# Small generic helpers


def _ordered_letters(value, allowed):
    value = str(value or "")
    return "".join([letter for letter in allowed if letter in value])


def normalise_report_format(report_format):
    """Return (canonical_name, extension) for markdown/text report formats."""
    key = str(report_format or "markdown").strip().lower()
    if key not in REPORT_FORMATS:
        raise ValueError("Invalid output format: %s. Use markdown/md or text/txt." % report_format)
    return REPORT_FORMATS[key]


def _json_safe(obj):
    """Convert tuples, sets, numpy scalars, etc. into ordinary JSON objects."""
    if obj is None or isinstance(obj, (str, int, float, bool)):
        if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
            return None
        return obj
    if isinstance(obj, dict):
        return {str(key): _json_safe(value) for key, value in obj.items()}
    if isinstance(obj, set):
        return sorted([_json_safe(value) for value in obj], key=lambda x: json.dumps(x, sort_keys=True))
    if isinstance(obj, (list, tuple)):
        return [_json_safe(value) for value in obj]
    try:
        return obj.item()  # numpy scalar
    except Exception:
        return str(obj)


def _now_iso():
    return datetime.now().astimezone().isoformat(timespec="seconds")


def _target_list(value):
    return sorted([str(v) for v in value if str(v).strip() != ""])


def _format_percent(value):
    return "%.1f%%" % (100.0 * float(value))


def _format_request_warning_summary(metadata):
    warnings = (metadata or {}).get("request_warnings", {}) or {}
    count = int(warnings.get("count", 0) or 0)
    if count <= 0:
        return "none"
    by_status = warnings.get("by_status", {}) or {}
    pieces = ["%s=%s" % (status, by_status[status]) for status in sorted(by_status)]
    if pieces:
        return "%d failed after retry (%s)" % (count, ", ".join(pieces))
    return "%d failed after retry" % count


def _escape_md(value):
    value = "" if value is None else str(value)
    return value.replace("|", "\\|").replace("\n", "<br>")


def _plain_url_gene(gene):
    return GENECARDS_URL % str(gene).strip()


def _markdown_gene_link(gene):
    gene = str(gene).strip()
    return "[%s](%s)" % (_escape_md(gene), GENECARDS_URL % gene)


# ---------------------------------------------------------------------------
# Metadata and option explanation


def explain_options(metadata):
    """Return English explanations for the options stored in metadata."""
    opts = metadata.get("options", {})
    search = _ordered_letters(opts.get("search", "A"), SEARCH_ORDER)
    report = _ordered_letters(opts.get("report", "A"), "AD")
    ofield = _ordered_letters(opts.get("ofield", "G"), "OGU")

    explanations = []
    if search:
        explanations.append(
            "Search evidence included: " + "; ".join(SEARCH_TEXT.get(letter, letter) for letter in search) + "."
        )
    if report:
        explanations.append(
            "Report sections requested: " + "; ".join(REPORT_TEXT.get(letter, letter) for letter in report) + "."
        )
    if ofield:
        explanations.append(
            "Target labels shown as: " + " + ".join(OFIELD_TEXT.get(letter, letter) for letter in ofield) + "."
        )

    report_format = opts.get("report_format")
    if report_format:
        explanations.append("Human-readable report format: %s." % report_format)

    if opts.get("similarity_search"):
        explanations.append(
            "Tier 1 enabled: similar ChEMBL molecules were retrieved and self matches were excluded; "
            "the similarity threshold was %.3f." % float(opts.get("similarity_threshold", 0.0))
        )
    else:
        explanations.append("Tier 1 disabled: only exact ChEMBL matches for the input molecules were used.")

    if opts.get("join_aggregates"):
        explanations.append("Joined T0+T1 summaries were also calculated for JSON/raw output.")

    rules = metadata.get("activity_rules", {})
    if rules:
        explanations.append(
            "Activity classification used rule profile '%s' from %s."
            % (rules.get("version", "unknown"), rules.get("source", "unknown"))
        )

    return explanations


def build_run_metadata(csanno_version, input_file=None, single_mol=None, mols=None,
                       chembl_server=None, do_sims=False, thres=None, search="A",
                       report="A", ofield="G", join_aggregates=False,
                       activity_rules=None, report_format="markdown"):
    """Create the metadata block required in the report and JSON."""
    mols = mols or {}
    activity_rules = activity_rules or {}
    fmt_name, _ = normalise_report_format(report_format)
    source = {
        "mode": "smiles" if single_mol is not None else "file",
        "file_name": input_file,
        "smiles": single_mol,
    }
    metadata = {
        "csanno_version": csanno_version,
        "query_datetime": _now_iso(),
        "source": source,
        "molecule_count": len(mols),
        "chembl_server": chembl_server,
        "options": {
            "similarity_search": bool(do_sims),
            "similarity_threshold": thres,
            "search": _ordered_letters(search, SEARCH_ORDER),
            "report": _ordered_letters(report, "AD"),
            "ofield": _ordered_letters(ofield, "OGU"),
            "join_aggregates": bool(join_aggregates),
            "report_format": fmt_name,
        },
        "activity_rules": {
            "version": activity_rules.get("activity_rules_version", "unknown"),
            "source": activity_rules.get("_source", "unknown"),
            "profile_name": activity_rules.get("profile_name", "unknown"),
            "profile_description": activity_rules.get("profile_description", ""),
        },
    }
    metadata["options_explained"] = explain_options(metadata)
    return metadata


# ---------------------------------------------------------------------------
# Target annotation extraction moved out of the main CSANNO script


def _activity_filter(search):
    act_filter = []
    if "A" in str(search or ""):
        act_filter.append("A")
    if "U" in str(search or ""):
        act_filter.append("NA")
    if "N" in str(search or ""):
        act_filter.append("N")
    return act_filter


def _classify_activity(act, activity_rules=None, activity_checker=None):
    if activity_checker is None:
        raise ValueError("make_reports requires activity_checker=check_activity from the main CSANNO module")
    return activity_checker(act, activity_rules=activity_rules)


def _single_biological_target_label(targets, tid, ofield):
    if tid not in targets:
        return None
    # Keep the historical CSANNO constraint: only one biological component.
    if len(targets[tid]) != 1:
        return None
    organism, gene, uniprot = targets[tid][0]
    fields = []
    if "O" in ofield:
        fields.append(str(organism).upper())
    if "G" in ofield:
        fields.append(str(gene).upper())
    if "U" in ofield:
        fields.append(str(uniprot).upper())
    output = "_".join(fields)
    if output.strip() == "":
        return None
    return output


def make_reports(mols, mol_ids, mol_acts, targets, search="A", ofield="G", activity_rules=None,
                 activity_checker=None):
    """Return annotations for Tier 0 exact matches and their target IDs.

    This preserves the historical CSANNO return shape:
        (mol_active_genes, good_targs)
    where good_targs is a set of (chembl_target_id, rendered_target_label).
    """
    good_targs = set()
    act_filter = _activity_filter(search)
    mol_active_genes = {}

    for mid in list((mols or {}).keys()):
        ick = mols[mid][1]
        if ick not in mol_ids:
            continue
        db_id = mol_ids[ick]
        for act in mol_acts.get(db_id, []):
            is_active = _classify_activity(act, activity_rules=activity_rules, activity_checker=activity_checker)
            if is_active not in act_filter:
                continue
            tid = act.get("target_id")
            output = _single_biological_target_label(targets, tid, ofield)
            if output is None:
                continue
            mol_active_genes.setdefault(mid, set()).add(output)
            good_targs.add((tid, output))

    return mol_active_genes, good_targs


def make_reports_x(mol_acts, targets, search="A", ofield="G", activity_rules=None, activity_checker=None):
    """Return annotations for Tier 1 ChEMBL analogue molecules and their target IDs."""
    good_targs = set()
    act_filter = _activity_filter(search)
    mol_active_genes = {}

    for mid in mol_acts or {}:
        for act in mol_acts[mid]:
            is_active = _classify_activity(act, activity_rules=activity_rules, activity_checker=activity_checker)
            if is_active not in act_filter:
                continue
            tid = act.get("target_id")
            output = _single_biological_target_label(targets, tid, ofield)
            if output is None:
                continue
            mol_active_genes.setdefault(mid, set()).add(output)
            good_targs.add((tid, output))

    return mol_active_genes, good_targs


def merge_T0T1(mols, mols_t0, mol_genes_T0, mols_t1, mol_genes_T1):
    """Return T0+T1 annotations grouped by original input molecule."""
    mol_genes_T01 = {}
    for mid in mols:
        genes = []
        if mid in mol_genes_T0:
            genes = list(mol_genes_T0[mid])
        if mid in mols_t1:
            for analogue_id in mols_t1[mid]:
                if analogue_id in mol_genes_T1:
                    genes = genes + list(mol_genes_T1[analogue_id])
        mol_genes_T01[mid] = set(genes)
    return mol_genes_T01


def merge_T1(mols, mols_t1, mol_genes_T1):
    """Return Tier 1 annotations only, grouped by original input molecule."""
    mol_genes_T1_by_input = {}
    for mid in mols:
        genes = []
        if mid in mols_t1:
            for analogue_id in mols_t1[mid]:
                if analogue_id in mol_genes_T1:
                    genes = genes + list(mol_genes_T1[analogue_id])
        mol_genes_T1_by_input[mid] = set(genes)
    return mol_genes_T1_by_input


# ---------------------------------------------------------------------------
# Legacy raw text reports moved out of the main CSANNO script


def gene_counting(mols, mol_active_genes):
    genes = []
    uniqs = set()
    for mid in mols:
        if mid in mol_active_genes:
            for gene in mol_active_genes[mid]:
                genes.append(gene)
            uniqs = uniqs | mol_active_genes[mid]
    gene_counts = []
    for gene in uniqs:
        gene_counts.append((genes.count(gene), gene))
    gene_counts.sort(key=lambda item: (-item[0], item[1]))
    return gene_counts


def write_report_aggregated(mols, mol_active_genes, fname=None):
    fil = open(fname, "wt") if fname is not None else None
    gene_counts = gene_counting(mols, mol_active_genes)
    denominator = len(mols)
    if denominator == 0:
        if fil is not None:
            fil.close()
        return

    for count, gene in gene_counts:
        if len(str(gene).strip()) == 0:
            continue
        line = "%s\t%5d\t%6.4f" % (gene, count, count / denominator)
        if fil is None:
            print(line)
        else:
            fil.write(line + "\n")
    if fil is not None:
        fil.close()


def write_report_detail(mols, mol_active_genes, fname=None, append=False):
    fmode = "at" if append else "wt"
    prefix = "_" if append else ""
    fil = open(fname, fmode) if fname is not None else None

    for mid in mols:
        line = "%s:" % (prefix + mid)
        if mid in mol_active_genes:
            for gene in sorted(mol_active_genes[mid]):
                if len(str(gene).strip()) > 0:
                    line += "\t%s" % gene
        if fil is None:
            print(line)
        else:
            fil.write(line + "\n")

    if fil is not None:
        fil.close()


def write_run_info(fname, params, silent=False, warn_func=None):
    if fname is None:
        return
    try:
        with open(fname, "wt") as fil:
            for key in sorted(params.keys()):
                fil.write("%s: %s\n" % (key, params[key]))
    except IOError as exc:
        if warn_func is not None:
            warn_func("Could not write run information file: %s" % exc, silent)
        elif not silent:
            print("WARNING: Could not write run information file: %s" % exc)


def write_reports_for_channel(mols, mol_active_genes, report, fname_aggregate, fname_detail,
                              write_files, to_screen):
    if "A" in report:
        if write_files:
            write_report_aggregated(mols, mol_active_genes, fname=fname_aggregate)
        if to_screen:
            write_report_aggregated(mols, mol_active_genes, fname=None)

    if "D" in report:
        if write_files:
            write_report_detail(mols, mol_active_genes, fname=fname_detail)
        if to_screen:
            write_report_detail(mols, mol_active_genes, fname=None)


# ---------------------------------------------------------------------------
# Structured payload used by Markdown/text/JSON outputs


def _molecule_records(mols):
    records = []
    for mol_id in mols:
        smiles, inchikey, input_activity = mols[mol_id]
        records.append({
            "molecule_id": mol_id,
            "smiles": smiles,
            "inchikey": inchikey,
            "input_activity": input_activity,
        })
    return records


def aggregate_rows(molecule_ids, annotations_by_molecule):
    molecule_ids = list(molecule_ids or [])
    denominator = len(molecule_ids)
    counts = Counter()
    for molecule_id in molecule_ids:
        for target in _target_list((annotations_by_molecule or {}).get(molecule_id, [])):
            counts[target] += 1

    rows = []
    for target, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
        rows.append({
            "target": target,
            "count": count,
            "fraction": (float(count) / denominator) if denominator else 0.0,
        })
    return rows


def detail_rows(molecule_ids, mols, annotations_by_molecule):
    rows = []
    for molecule_id in list(molecule_ids or []):
        mol_entry = (mols or {}).get(molecule_id)
        if mol_entry is None:
            smiles, inchikey, input_activity = None, None, None
        else:
            smiles, inchikey, input_activity = mol_entry
        targets = _target_list((annotations_by_molecule or {}).get(molecule_id, []))
        rows.append({
            "molecule_id": molecule_id,
            "smiles": smiles,
            "inchikey": inchikey,
            "input_activity": input_activity,
            "annotation_count": len(targets),
            "annotations": targets,
        })
    return rows


def _median(values):
    values = sorted(values)
    n = len(values)
    if n == 0:
        return 0.0
    mid = n // 2
    if n % 2:
        return float(values[mid])
    return (float(values[mid - 1]) + float(values[mid])) / 2.0


def _detail_stats(rows):
    counts = [row.get("annotation_count", 0) for row in rows]
    if not counts:
        return {
            "molecule_count": 0,
            "molecules_with_annotations": 0,
            "average_annotations_per_molecule": 0.0,
            "median_annotations_per_molecule": 0.0,
            "max_annotations_per_molecule": 0,
            "most_annotated_molecules": [],
        }
    most = sorted(rows, key=lambda row: (-row.get("annotation_count", 0), row.get("molecule_id", "")))[:5]
    return {
        "molecule_count": len(rows),
        "molecules_with_annotations": sum(1 for count in counts if count > 0),
        "average_annotations_per_molecule": sum(counts) / float(len(counts)),
        "median_annotations_per_molecule": _median(counts),
        "max_annotations_per_molecule": max(counts),
        "most_annotated_molecules": [
            {"molecule_id": row["molecule_id"], "annotation_count": row.get("annotation_count", 0)}
            for row in most if row.get("annotation_count", 0) > 0
        ],
    }


def _top_targets(rows, limit=10):
    return rows[:limit]


def _evidence_language(search):
    search = _ordered_letters(search, SEARCH_ORDER)
    if search == "A":
        return {
            "annotation_phrase": "active target annotation",
            "annotation_plural": "active target annotations",
            "tier0_action": "direct ChEMBL active evidence",
            "tier1_action": "analogue-derived probable target evidence",
            "breadth_name": "promiscuity",
            "common_t0": "most common targets",
            "common_t1": "most common probable targets",
            "tier1_note": "Tier 1 findings are analogue-derived suggestions, not direct activity claims for the input molecules.",
        }
    if search == "N":
        return {
            "annotation_phrase": "negative / non-active annotation",
            "annotation_plural": "negative / non-active annotations",
            "tier0_action": "direct ChEMBL negative evidence",
            "tier1_action": "analogue-derived negative evidence",
            "breadth_name": "breadth of negative evidence",
            "common_t0": "most frequent negative annotations",
            "common_t1": "most frequent analogue-derived negative annotations",
            "tier1_note": "Tier 1 negative findings mean similar molecules had reported non-active evidence; they should not be read as predicted absence of activity for the input molecules.",
        }
    if search == "U":
        return {
            "annotation_phrase": "uncertain annotation",
            "annotation_plural": "uncertain annotations",
            "tier0_action": "direct ChEMBL uncertain evidence",
            "tier1_action": "analogue-derived uncertain evidence",
            "breadth_name": "breadth of uncertain evidence",
            "common_t0": "most frequent uncertain annotations",
            "common_t1": "most frequent analogue-derived uncertain annotations",
            "tier1_note": "Tier 1 uncertain findings identify where similar molecules had evidence that could not be classified as active or negative.",
        }
    return {
        "annotation_phrase": "selected annotation",
        "annotation_plural": "selected annotations",
        "tier0_action": "direct ChEMBL evidence matching the selected search classes (%s)" % search,
        "tier1_action": "analogue-derived evidence matching the selected search classes (%s)" % search,
        "breadth_name": "annotation breadth",
        "common_t0": "most frequent selected annotations",
        "common_t1": "most frequent analogue-derived selected annotations",
        "tier1_note": "Mixed search classes are shown as selected annotations; interpret them according to the A/N/U classes requested.",
    }


def _format_target_row(row, denominator_label="molecules"):
    return "%s (%d %s, %s)" % (
        row.get("target", ""),
        row.get("count", 0),
        denominator_label,
        _format_percent(row.get("fraction", 0.0)),
    )


def _render_top_targets(rows, denominator_label):
    if not rows:
        return "none found"
    return "; ".join(_format_target_row(row, denominator_label) for row in rows)


def _render_top_molecules(rows):
    if not rows:
        return "none"
    return "; ".join("%s (%d)" % (row["molecule_id"], row["annotation_count"]) for row in rows)


def target_reference_records(target_ids, targets, tier):
    """Create explicit target metadata from (target_id, rendered_label) sets."""
    records = []
    seen = set()
    for ref in sorted(list(target_ids or []), key=lambda value: str(value)):
        if not isinstance(ref, (list, tuple)) or len(ref) < 2:
            continue
        target_id, label = ref[0], ref[1]
        comps = (targets or {}).get(target_id, [])
        organism = gene = uniprot = ""
        if len(comps) == 1:
            organism, gene, uniprot = comps[0]
        key = (tier, target_id, label, organism, gene, uniprot)
        if key in seen:
            continue
        seen.add(key)
        records.append({
            "tier": tier,
            "chembl_target_id": target_id,
            "label": label,
            "organism": organism,
            "gene": str(gene).upper() if gene is not None else "",
            "uniprot": str(uniprot).upper() if uniprot is not None else "",
        })
    return records


def build_target_metadata(results):
    records = []
    records.extend(target_reference_records(results.get("target_ids"), results.get("targets"), "T0"))
    records.extend(target_reference_records(results.get("target_ids_x"), results.get("targets_x"), "T1"))
    # Deduplicate by gene/uniprot/label across tiers, but preserve tier membership.
    merged = {}
    for rec in records:
        key = (rec.get("chembl_target_id"), rec.get("label"), rec.get("gene"), rec.get("uniprot"))
        if key not in merged:
            merged[key] = dict(rec)
            merged[key]["tiers"] = [rec.get("tier")]
        else:
            if rec.get("tier") not in merged[key]["tiers"]:
                merged[key]["tiers"].append(rec.get("tier"))
    return sorted(merged.values(), key=lambda rec: (rec.get("gene", ""), rec.get("label", "")))


# ---------------------------------------------------------------------------
# Gene glossary retrieval


def _load_glossary_cache(cache_file=DEFAULT_GLOSSARY_CACHE_JSON):
    if cache_file is None or not os.path.isfile(cache_file):
        return {}
    try:
        with open(cache_file, "rt") as fil:
            data = json.load(fil)
        return data if isinstance(data, dict) else {}
    except (IOError, ValueError):
        return {}


def _save_glossary_cache(cache, cache_file=DEFAULT_GLOSSARY_CACHE_JSON):
    if cache_file is None:
        return
    try:
        with open(cache_file, "wt") as fil:
            json.dump(cache, fil, indent=2, sort_keys=True)
            fil.write("\n")
    except IOError:
        pass


def _extract_uniprot_protein_name(result):
    description = result.get("proteinDescription", {}) if isinstance(result, dict) else {}
    for key in ["recommendedName", "submissionNames"]:
        value = description.get(key)
        if isinstance(value, dict):
            full = value.get("fullName", {})
            if isinstance(full, dict) and full.get("value"):
                return full.get("value")
        elif isinstance(value, list):
            for item in value:
                full = item.get("fullName", {}) if isinstance(item, dict) else {}
                if isinstance(full, dict) and full.get("value"):
                    return full.get("value")
    # Fallback to a compact protein existence/name if UniProt changes fields.
    if description.get("recommendedName"):
        return str(description.get("recommendedName"))
    return ""


def _fetch_uniprot_gene_summary(gene, timeout=GLOSSARY_REQUEST_TIMEOUT):
    gene = str(gene or "").strip().upper()
    if not gene or requests is None:
        return None

    queries = [
        "(gene_exact:%s) AND (reviewed:true) AND (organism_id:9606)" % gene,
        "(gene_exact:%s) AND (reviewed:true)" % gene,
        "gene_exact:%s" % gene,
    ]
    for query in queries:
        try:
            response = requests.get(
                UNIPROT_SEARCH_URL,
                params={
                    "query": query,
                    "fields": "accession,protein_name,gene_names,organism_name",
                    "format": "json",
                    "size": 1,
                },
                timeout=timeout,
            )
            response.raise_for_status()
            data = response.json()
        except Exception:
            # Network/service failures should not make every gene wait for a timeout.
            return {"_fetch_failed": True}
        results = data.get("results", []) if isinstance(data, dict) else []
        if not results:
            continue
        result = results[0]
        accession = result.get("primaryAccession") or ""
        organism = ""
        if isinstance(result.get("organism"), dict):
            organism = result.get("organism", {}).get("scientificName", "")
        protein_name = _extract_uniprot_protein_name(result)
        if accession or protein_name:
            return {
                "gene": gene,
                "description": protein_name or "description unavailable",
                "source": "UniProt",
                "uniprot_accession": accession,
                "organism": organism,
                "uniprot_url": (UNIPROT_ENTRY_URL % accession) if accession else None,
                "genecards_url": GENECARDS_URL % gene,
            }
    return None


def _fallback_gene_summary(gene):
    gene = str(gene or "").strip().upper()
    return {
        "gene": gene,
        "description": "description unavailable",
        "source": "not found",
        "uniprot_accession": None,
        "organism": None,
        "uniprot_url": None,
        "genecards_url": GENECARDS_URL % gene,
    }


def _looks_like_gene_symbol(token):
    token = str(token or "").strip().upper()
    if token in COMMON_NON_GENE_TOKENS:
        return False
    if len(token) < 2 or len(token) > 18:
        return False
    if re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$", token):  # common UniProt accession pattern
        return False
    if re.match(r"^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$", token):
        return False
    return re.match(r"^[A-Z][A-Z0-9-]*[0-9A-Z]$", token) is not None


def _collect_gene_symbols_from_labels(sections):
    genes = set()
    labels = []
    for section in sections or []:
        for row in section.get("rows", []):
            if section.get("kind") == "aggregate":
                labels.append(row.get("target", ""))
            elif section.get("kind") == "detail":
                labels.extend(row.get("annotations", []))
    for label in labels:
        for token in re.split(r"[_\s,/;:|]+", str(label or "").upper()):
            if _looks_like_gene_symbol(token):
                genes.add(token)
    return genes


def build_gene_glossary(target_metadata, sections, cache_file=DEFAULT_GLOSSARY_CACHE_JSON,
                        fetch=True, silent=False):
    genes = set()
    for rec in target_metadata or []:
        gene = str(rec.get("gene") or "").strip().upper()
        if gene and _looks_like_gene_symbol(gene):
            genes.add(gene)
    if not genes:
        genes = _collect_gene_symbols_from_labels(sections)

    cache = _load_glossary_cache(cache_file)
    glossary = []
    changed = False
    remote_lookup_disabled = False
    for gene in sorted(genes):
        record = cache.get(gene)
        if not isinstance(record, dict):
            record = None
        if record is None and fetch and not remote_lookup_disabled:
            record = _fetch_uniprot_gene_summary(gene)
            if isinstance(record, dict) and record.get("_fetch_failed"):
                remote_lookup_disabled = True
                record = _fallback_gene_summary(gene)
            elif record is None:
                record = _fallback_gene_summary(gene)
            cache[gene] = record
            changed = True
        elif record is None:
            record = _fallback_gene_summary(gene)
        glossary.append(record)

    if changed:
        _save_glossary_cache(cache, cache_file)
    return glossary


def build_report_payload(results, metadata, report="A", include_sims=False, include_joined=False,
                         include_glossary=True, glossary_cache_file=DEFAULT_GLOSSARY_CACHE_JSON):
    """Build a structured report payload before rendering text/Markdown/JSON."""
    results = results or {}
    metadata = dict(metadata or {})
    options = metadata.get("options", {})
    report = _ordered_letters(report or options.get("report", "A"), "AD") or "A"
    search = options.get("search", "A")

    mols = results.get("mols") or {}
    input_ids = list(mols.keys())
    t0_annotations = results.get("active_genes") or {}
    t0_agg = aggregate_rows(input_ids, t0_annotations)
    t0_detail = detail_rows(input_ids, mols, t0_annotations)

    sections = []
    if "A" in report:
        sections.append({
            "id": "tier0_aggregate",
            "title": "Tier 0 aggregate - exact ChEMBL matches",
            "kind": "aggregate",
            "denominator": len(input_ids),
            "denominator_label": "input molecules",
            "rows": t0_agg,
        })
    if "D" in report:
        sections.append({
            "id": "tier0_detail",
            "title": "Tier 0 detail - exact ChEMBL matches by input molecule",
            "kind": "detail",
            "rows": t0_detail,
        })

    tier1 = None
    if include_sims and results.get("mols_x") is not None:
        analogue_ids = list(results.get("mols_x") or [])
        t1_pool_annotations = results.get("active_genes_x") or {}
        t1_by_input = results.get("active_genes_T1_by_input") or {}
        t1_joined = results.get("active_genes_T01") or {}
        t1_pool_agg = aggregate_rows(analogue_ids, t1_pool_annotations)
        t1_by_input_agg = aggregate_rows(input_ids, t1_by_input)
        t1_detail = detail_rows(input_ids, mols, t1_by_input)
        tier1 = {
            "analogue_count": len(analogue_ids),
            "analogue_pool_aggregate": t1_pool_agg,
            "by_input_aggregate": t1_by_input_agg,
            "by_input_detail": t1_detail,
            "joined_T0_T1_by_input_detail": detail_rows(input_ids, mols, t1_joined),
            "joined_T0_T1_by_input_aggregate": aggregate_rows(input_ids, t1_joined),
        }
        if "A" in report:
            sections.append({
                "id": "tier1_aggregate",
                "title": "Tier 1 aggregate — similar ChEMBL analogue pool",
                "kind": "aggregate",
                "denominator": len(analogue_ids),
                "denominator_label": "similar ChEMBL molecules",
                "rows": t1_pool_agg,
            })
        if "D" in report:
            sections.append({
                "id": "tier1_detail",
                "title": "Tier 1 detail — analogue-derived annotations grouped by input molecule",
                "kind": "detail",
                "rows": t1_detail,
            })

    lang = _evidence_language(search)
    t0_stats = _detail_stats(t0_detail)
    digest = {
        "evidence_language": lang,
        "tier0": {
            "stats": t0_stats,
            "top_targets": _top_targets(t0_agg, 10),
        },
        "tier1": None,
    }

    if tier1 is not None:
        t1_stats = _detail_stats(tier1["by_input_detail"])
        analogue_rows = detail_rows(list(results.get("mols_x") or []), {}, results.get("active_genes_x") or {})
        analogue_stats = _detail_stats(analogue_rows)
        t0_targets = set(row["target"] for row in t0_agg)
        t1_targets = set(row["target"] for row in tier1["by_input_aggregate"])
        digest["tier1"] = {
            "stats_by_input": t1_stats,
            "stats_by_analogue_pool": analogue_stats,
            "top_probable_targets_by_input": _top_targets(tier1["by_input_aggregate"], 10),
            "top_targets_in_analogue_pool": _top_targets(tier1["analogue_pool_aggregate"], 10),
            "new_targets_vs_tier0": sorted(t1_targets - t0_targets),
        }

    target_metadata = build_target_metadata(results)
    glossary = build_gene_glossary(target_metadata, sections, cache_file=glossary_cache_file) if include_glossary else []

    raw_results = {
        "molecules": _molecule_records(mols),
        "target_metadata": target_metadata,
        "tier0": {
            "aggregate": t0_agg,
            "detail": t0_detail,
        },
        "tier1": tier1,
        "internal_result_payload": _json_safe(results),
    }

    payload = {
        "metadata": _json_safe(metadata),
        "digest": _json_safe(digest),
        "glossary": _json_safe(glossary),
        "sections": _json_safe(sections),
        "raw_results": _json_safe(raw_results),
    }
    return payload


# ---------------------------------------------------------------------------
# Text rendering


def _render_text_glossary(glossary):
    lines = []
    lines.append("Gene glossary")
    lines.append("-------------")
    if not glossary:
        lines.append("No gene symbols could be identified for glossary lookup.")
        return lines
    for rec in glossary:
        gene = rec.get("gene", "")
        desc = rec.get("description") or "description unavailable"
        accession = rec.get("uniprot_accession")
        if accession:
            source = "from UniProt entry %s" % accession
        elif rec.get("source") == "UniProt":
            source = "from UniProt"
        else:
            source = "not found in UniProt"
        lines.append("%s - %s (%s). GeneCards: %s" % (gene, desc, source, rec.get("genecards_url") or _plain_url_gene(gene)))
    return lines


def render_text_report(payload):
    metadata = payload.get("metadata", {})
    digest = payload.get("digest", {})
    lang = digest.get("evidence_language", _evidence_language("A"))
    source = metadata.get("source", {})

    lines = []
    lines.append("CSANNO readable report")
    lines.append("======================")
    lines.append("")
    lines.append("Metadata")
    lines.append("--------")
    lines.append("Query date/time: %s" % metadata.get("query_datetime", "unknown"))
    if source.get("mode") == "smiles":
        lines.append("Input SMILES: %s" % source.get("smiles"))
    else:
        lines.append("Input file: %s" % source.get("file_name"))
    lines.append("Molecule count: %s" % metadata.get("molecule_count", "unknown"))
    lines.append("CSANNO version: %s" % metadata.get("csanno_version", "unknown"))
    lines.append("ChEMBL server: %s" % metadata.get("chembl_server", "unknown"))
    lines.append("ChEMBL request warnings: %s" % _format_request_warning_summary(metadata))
    lines.append("")
    lines.append("CSANNO options, in English")
    lines.append("--------------------------")
    for explanation in metadata.get("options_explained", []):
        lines.append("- " + explanation)
    lines.append("")
    lines.append("Comprehensive digest")
    lines.append("--------------------")

    t0 = digest.get("tier0", {})
    t0_stats = t0.get("stats", {})
    n0 = t0_stats.get("molecule_count", 0)
    lines.append(
        "Tier 0 (%s): %d/%d input molecules had at least one %s. "
        "The median was %.1f annotations per molecule, with a maximum of %d."
        % (
            lang.get("tier0_action", "direct evidence"),
            t0_stats.get("molecules_with_annotations", 0),
            n0,
            lang.get("annotation_phrase", "annotation"),
            t0_stats.get("median_annotations_per_molecule", 0.0),
            t0_stats.get("max_annotations_per_molecule", 0),
        )
    )
    lines.append(
        "Tier 0 %s: %s."
        % (lang.get("breadth_name", "annotation breadth"), _render_top_molecules(t0_stats.get("most_annotated_molecules", [])))
    )
    lines.append(
        "Tier 0 %s: %s."
        % (lang.get("common_t0", "most frequent annotations"), _render_top_targets(t0.get("top_targets", []), "input molecules"))
    )

    t1 = digest.get("tier1")
    if t1 is None:
        lines.append("Tier 1: not run; no similarity search was requested.")
    else:
        t1_input = t1.get("stats_by_input", {})
        t1_pool = t1.get("stats_by_analogue_pool", {})
        lines.append(
            "Tier 1 (%s): %d similar ChEMBL molecules were evaluated; %d of them had at least one %s."
            % (
                lang.get("tier1_action", "analogue evidence"),
                t1_pool.get("molecule_count", 0),
                t1_pool.get("molecules_with_annotations", 0),
                lang.get("annotation_phrase", "annotation"),
            )
        )
        lines.append(
            "Grouped back to the input file, %d/%d input molecules had at least one Tier 1 %s; "
            "the median was %.1f annotations per molecule."
            % (
                t1_input.get("molecules_with_annotations", 0),
                t1_input.get("molecule_count", 0),
                lang.get("annotation_phrase", "annotation"),
                t1_input.get("median_annotations_per_molecule", 0.0),
            )
        )
        lines.append(
            "Tier 1 %s by input molecule: %s."
            % (lang.get("common_t1", "frequent analogue-derived annotations"), _render_top_targets(t1.get("top_probable_targets_by_input", []), "input molecules"))
        )
        lines.append(
            "Tier 1 strongest signals in the analogue pool: %s."
            % _render_top_targets(t1.get("top_targets_in_analogue_pool", []), "similar molecules")
        )
        new_targets = t1.get("new_targets_vs_tier0", [])
        if new_targets:
            lines.append(
                "Tier 1 annotations not seen in Tier 0: %s."
                % ", ".join(new_targets[:25])
                + (" ..." if len(new_targets) > 25 else "")
            )
        lines.append(lang.get("tier1_note", ""))

    lines.append("")
    lines.extend(_render_text_glossary(payload.get("glossary", [])))

    lines.append("")
    lines.append("Raw data")
    lines.append("--------")
    for section in payload.get("sections", []):
        lines.append("")
        lines.append(section.get("title", section.get("id", "Section")))
        lines.append("~" * len(lines[-1]))
        if section.get("kind") == "aggregate":
            lines.append("Target\tMolecule count\tFraction")
            for row in section.get("rows", []):
                lines.append("%s\t%d\t%.4f" % (row.get("target", ""), row.get("count", 0), row.get("fraction", 0.0)))
            if not section.get("rows"):
                lines.append("No annotations found.")
        elif section.get("kind") == "detail":
            lines.append("Molecule\tAnnotation count\tAnnotations")
            for row in section.get("rows", []):
                annotations = ", ".join(row.get("annotations", []))
                lines.append("%s\t%d\t%s" % (row.get("molecule_id", ""), row.get("annotation_count", 0), annotations))
            if not section.get("rows"):
                lines.append("No molecules found.")

    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Markdown rendering


def _render_md_top_targets(rows, denominator_label):
    if not rows:
        return "none found"
    parts = []
    for row in rows:
        parts.append("**%s** (%d %s, %s)" % (
            _escape_md(row.get("target", "")), row.get("count", 0), denominator_label,
            _format_percent(row.get("fraction", 0.0))))
    return "; ".join(parts)


def _render_md_glossary(glossary):
    lines = []
    lines.append("## Gene glossary")
    lines.append("")
    if not glossary:
        lines.append("No gene symbols could be identified for glossary lookup.")
        lines.append("")
        return lines
    for rec in glossary:
        gene = rec.get("gene", "")
        desc = rec.get("description") or "description unavailable"
        accession = rec.get("uniprot_accession")
        gene_link = _markdown_gene_link(gene)
        if accession:
            uniprot_link = "[UniProt entry %s](%s)" % (_escape_md(accession), UNIPROT_ENTRY_URL % accession)
            source = "from %s" % uniprot_link
        elif rec.get("source") == "UniProt":
            source = "from UniProt"
        else:
            source = "not found in UniProt"
        lines.append("- %s - %s (%s)" % (gene_link, _escape_md(desc), source))
    lines.append("")
    return lines


def _render_markdown_table(headers, rows):
    lines = []
    lines.append("| " + " | ".join(_escape_md(h) for h in headers) + " |")
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        lines.append("| " + " | ".join(_escape_md(cell) for cell in row) + " |")
    return lines


def render_markdown_report(payload):
    metadata = payload.get("metadata", {})
    digest = payload.get("digest", {})
    lang = digest.get("evidence_language", _evidence_language("A"))
    source = metadata.get("source", {})

    lines = []
    lines.append("# CSANNO readable report")
    lines.append("")
    lines.append("## Metadata")
    lines.append("")
    metadata_rows = [
        ("Query date/time", metadata.get("query_datetime", "unknown")),
        ("Input SMILES" if source.get("mode") == "smiles" else "Input file", source.get("smiles") if source.get("mode") == "smiles" else source.get("file_name")),
        ("Molecule count", metadata.get("molecule_count", "unknown")),
        ("CSANNO version", metadata.get("csanno_version", "unknown")),
        ("ChEMBL server", metadata.get("chembl_server", "unknown")),
        ("ChEMBL request warnings", _format_request_warning_summary(metadata)),
    ]
    lines.extend(_render_markdown_table(["Field", "Value"], metadata_rows))
    lines.append("")

    lines.append("## CSANNO options, in English")
    lines.append("")
    for explanation in metadata.get("options_explained", []):
        lines.append("- " + _escape_md(explanation))
    lines.append("")

    lines.append("## Comprehensive digest")
    lines.append("")
    t0 = digest.get("tier0", {})
    t0_stats = t0.get("stats", {})
    n0 = t0_stats.get("molecule_count", 0)
    lines.append(
        "**Tier 0 (%s).** %d/%d input molecules had at least one %s. "
        "The median was %.1f annotations per molecule, with a maximum of %d."
        % (
            _escape_md(lang.get("tier0_action", "direct evidence")),
            t0_stats.get("molecules_with_annotations", 0),
            n0,
            _escape_md(lang.get("annotation_phrase", "annotation")),
            t0_stats.get("median_annotations_per_molecule", 0.0),
            t0_stats.get("max_annotations_per_molecule", 0),
        )
    )
    lines.append("")
    lines.append(
        "**Tier 0 %s:** %s."
        % (_escape_md(lang.get("breadth_name", "annotation breadth")), _escape_md(_render_top_molecules(t0_stats.get("most_annotated_molecules", []))))
    )
    lines.append("")
    lines.append(
        "**Tier 0 %s:** %s."
        % (_escape_md(lang.get("common_t0", "most frequent annotations")), _render_md_top_targets(t0.get("top_targets", []), "input molecules"))
    )
    lines.append("")

    t1 = digest.get("tier1")
    if t1 is None:
        lines.append("**Tier 1:** not run; no similarity search was requested.")
        lines.append("")
    else:
        t1_input = t1.get("stats_by_input", {})
        t1_pool = t1.get("stats_by_analogue_pool", {})
        lines.append(
            "**Tier 1 (%s).** %d similar ChEMBL molecules were evaluated; %d of them had at least one %s."
            % (
                _escape_md(lang.get("tier1_action", "analogue evidence")),
                t1_pool.get("molecule_count", 0),
                t1_pool.get("molecules_with_annotations", 0),
                _escape_md(lang.get("annotation_phrase", "annotation")),
            )
        )
        lines.append("")
        lines.append(
            "Grouped back to the input file, %d/%d input molecules had at least one Tier 1 %s; "
            "the median was %.1f annotations per molecule."
            % (
                t1_input.get("molecules_with_annotations", 0),
                t1_input.get("molecule_count", 0),
                _escape_md(lang.get("annotation_phrase", "annotation")),
                t1_input.get("median_annotations_per_molecule", 0.0),
            )
        )
        lines.append("")
        lines.append(
            "**Tier 1 %s by input molecule:** %s."
            % (_escape_md(lang.get("common_t1", "frequent analogue-derived annotations")), _render_md_top_targets(t1.get("top_probable_targets_by_input", []), "input molecules"))
        )
        lines.append("")
        lines.append(
            "**Tier 1 strongest signals in the analogue pool:** %s."
            % _render_md_top_targets(t1.get("top_targets_in_analogue_pool", []), "similar molecules")
        )
        lines.append("")
        new_targets = t1.get("new_targets_vs_tier0", [])
        if new_targets:
            lines.append(
                "**Tier 1 annotations not seen in Tier 0:** %s."
                % _escape_md(", ".join(new_targets[:25]) + (" ..." if len(new_targets) > 25 else ""))
            )
            lines.append("")
        lines.append(_escape_md(lang.get("tier1_note", "")))
        lines.append("")

    lines.extend(_render_md_glossary(payload.get("glossary", [])))

    lines.append("## Raw data")
    lines.append("")
    for section in payload.get("sections", []):
        lines.append("### " + _escape_md(section.get("title", section.get("id", "Section"))))
        lines.append("")
        if section.get("kind") == "aggregate":
            rows = [[row.get("target", ""), row.get("count", 0), "%.4f" % row.get("fraction", 0.0)]
                    for row in section.get("rows", [])]
            if rows:
                lines.extend(_render_markdown_table(["Target", "Molecule count", "Fraction"], rows))
            else:
                lines.append("No annotations found.")
        elif section.get("kind") == "detail":
            rows = [[row.get("molecule_id", ""), row.get("annotation_count", 0), ", ".join(row.get("annotations", []))]
                    for row in section.get("rows", [])]
            if rows:
                lines.extend(_render_markdown_table(["Molecule", "Annotation count", "Annotations"], rows))
            else:
                lines.append("No molecules found.")
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# File writer


def write_combined_report(results, metadata, report_fname=None, json_fname=None, report="A",
                          include_sims=False, include_joined=False, to_screen=False, silent=False,
                          report_format="markdown", txt_fname=None,
                          glossary_cache_file=DEFAULT_GLOSSARY_CACHE_JSON):
    """Write the combined report and JSON payload. Returns (rendered_report, payload)."""
    fmt_name, _ = normalise_report_format(report_format)
    # Backward compatibility with the previous function argument name.
    if report_fname is None and txt_fname is not None:
        report_fname = txt_fname

    payload = build_report_payload(results, metadata, report=report, include_sims=include_sims,
                                   include_joined=include_joined, include_glossary=True,
                                   glossary_cache_file=glossary_cache_file)
    if fmt_name == "markdown":
        rendered = render_markdown_report(payload)
    else:
        rendered = render_text_report(payload)

    if report_fname is not None:
        with open(report_fname, "wt", encoding="utf8") as fil:
            fil.write(rendered)
    if json_fname is not None:
        with open(json_fname, "wt") as fil:
            json.dump(payload, fil, indent=2, sort_keys=True)
            fil.write("\n")
    if to_screen:
        print(rendered)
    return rendered, payload
