import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_style('whitegrid')


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

    # if pth.startswith("splicing|SRE"):
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


def load_threes_ranks_table(fpath):
    df = pd.read_csv(fpath, sep="\t")

    # melt rows with multiple PATHOMECHANISM entries
    df['PATHOMECHANISM'] = df['PATHOMECHANISM'].str.split(";")
    df = df.explode("PATHOMECHANISM")

    # remove lines with non-splicing pathomechanism
    df = df.loc[df.PATHOMECHANISM.str.startswith("splicing"), :]

    df['PATHOGRP'] = df.PATHOMECHANISM.map(group_pathomechanism_for_threes_evaluation)
    return df


def load_threes_phenopacket_scores_table(fpath):
    df = pd.read_csv(fpath, sep="\t")

    df = df.loc[df.PATHOMECHANISM.str.startswith("splicing"), :]
    df['PATHOGRP'] = df.PATHOMECHANISM.apply(simplify_pathomechanism)
    # df = df.rename(columns={'PHENOPACKET':'CASE'})
    # df['CASE'] = df['CASE'].apply(lambda x: x[:x.find('.json')])
    return df


def plot_single_color(w_splicing, wo_splicing, ax=None):
    if not ax:
        fig, ax = plt.subplots(figsize=(12, 8))

    max_val = max([w_splicing.max(), wo_splicing.max()])
    line_points = np.linspace(0, max_val, 10)
    line = ax.plot(line_points, line_points, "--", color="grey", alpha=0.5)
    sct = ax.scatter(w_splicing, wo_splicing, color='red', alpha=0.5)

    xl = ax.set_xlabel("WITH SPLICING")
    yl = ax.set_ylabel("WITHOUT SPLICING")
    title = ax.set_title("Causal gene ranks with or without\nSPLICING pathogenicity source",
                         size="xx-large", fontweight="bold")


def plot_multicolor(w_splicing, wo_splicing, hue, ho, ax=None):
    if not ax:
        fig, ax = plt.subplots(figsize=(12, 8), dpi=150)

    max_val = max([w_splicing.max(), wo_splicing.max()])
    line_points = np.linspace(0, max_val, 10)
    line = ax.plot(line_points, line_points, "--", color="grey", alpha=0.5)
    sct = sns.scatterplot(w_splicing, wo_splicing, hue=hue, hue_order=ho, alpha=0.5, ax=ax)

    xl = ax.set_xlabel("WITH SPLICING")
    yl = ax.set_ylabel("WITHOUT SPLICING")
    title = ax.set_title("Causal gene ranks with or without\nSPLICING pathogenicity source",
                         size="large", fontweight="bold")


def load_spliceai_table(fpath: str):
    """Load SpliceAI result table created by `spliceai-benchmark` code.
    
    - create the DataFrame
    - remove records representing variant impact to different gene (variant close to two overlapping transcripts 
    of different genes will be annotated with respect to both transcripts)
    - remove records not representing `splicing|*` pathomechanism
    """
    dtypes = {'DS_AG': np.float, 'DS_AL': np.float, 'DS_DG': np.float, 'DS_DL': np.float}
    sd = pd.read_csv(fpath, sep='\t', dtype=dtypes, na_values={'.', np.nan})
    # remove records regarding transcript of a noncausal gene. The `CASE` column looks like `PMID:11828341-Ishii-2002-GLA-proband_1`
    # it contains gene symbol `GLA`.
    sd = sd.loc[sd.apply(lambda x: x.SYMBOL in x.CASE, axis=1), :]

    # remove non-splicing entries
    sd = sd.loc[sd.PATHOMECHANISM.str.startswith("splicing"), :]

    sd['PATHOGRP'] = sd.PATHOMECHANISM.apply(simplify_pathomechanism)
    print("Loaded {} records".format(len(sd)))
    return sd.reset_index(drop=True)
