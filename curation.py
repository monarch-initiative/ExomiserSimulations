import pandas as pd

def group_consequence(element):
    """Map CONSEQUENCE to friedlier values."""
    if isinstance(element, pd.Series):
        csq = element['cs']
    else:
        csq = element
        
    if isinstance(csq, float):
        return "None"    
    elif csq == "Alternative/cryptic 5' splice site":
        return "5_CSS"
    elif csq == "Alternative/cryptic 3' splice site":
        return "3_CSS"
    elif csq in ("Exon skipping", "Multiple exons skipped"):
        return "EXON_SKIP"
    elif csq in ("Pseudoexon inclusion"):
        return "PSEUDOEXON"
    return csq


def group_pathomechanism(element):
    """Map pathomechanism to friendlier values."""
    if isinstance(element, pd.Series):
        pth = element['pm']
    else:
        pth = element
    
    #if pth.startswith("splicing|SRE"):
    #    return "splicing|SRE"
    if pth.startswith("splicing|PPT"):
        return "splicing|3ss|disrupted"
    elif pth.startswith("splicing|branchpoint"):
        return "splicing|3ss|disrupted"
    elif pth.startswith("coding"):
        return "coding"
    return pth


def group_pathomechanism_for_threes_evaluation(entry):
    """Map pathomechanism to friendlier values."""
    # The pathomechanism can consist of two semicolon-separated entries
    # e.g. `splicing|3ss|disrupted;coding|missense`
    # We are interested only in those that begin with `splicing`.
    tokens = [k for k in entry.split(";") if k.startswith("splicing")]
    simplified = set([simplify_pathomechanism(k) for k in tokens])
    prioritized = prioritize(simplified)
    return ";".join(prioritized)


def simplify_pathomechanism(element):
    """Discretize pathomechanism labels."""
    if element in ("splicing|branchpoint", "splicing|PPT|disrupted"):
        return "PPT"
    if element.startswith("splicing|SRE"):
        return "SRE"
    else:
        return element[9:]


def prioritize(entries):
    if len(entries) == 1:
        # nothing to prioritize
        return entries
    elif len(entries) > 2:
        # unable to prioritize 3 and more entries
        return entries
    else:
        # we have 2 entries
        we_have_sre = any(["SRE" in e for e in entries])
        we_do_not_have_sre = any(["SRE" not in e for e in entries])
        if we_have_sre and we_do_not_have_sre:
            for e in entries:
                if not e.startswith("SRE"):
                    return e
        else:
            return entries
