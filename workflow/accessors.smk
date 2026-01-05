import os

###########################################################################
# Unified Accessors for Assemblies, Purged Assemblies, and Fastq Files
###########################################################################

def normalize_sample(wc):
    return wc.sample

# FASTQs from config
def get_fq1_abs(wc):
    return sample_fq1[normalize_sample(wc)]

def get_fq2_abs(wc):
    return sample_fq2[normalize_sample(wc)]

# Raw assemblies
def get_asm_hifiasm_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"Assemblies/{s}/HIFIASM/{s}.fasta")

def get_asm_lja_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"Assemblies/{s}/LJA/assembly.fasta")

# Sorted purged assemblies
def get_sorted_purged_hifiasm_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"Assemblies/{s}/PURGE_DUPS/purged-hifiasm.sorted.fa")

def get_sorted_purged_lja_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"Assemblies/{s}/PURGE_DUPS/purged-lja.sorted.fa")

###########################################################################
# Unified Purged Assembly Selection (prefers RagTag > sorted > raw)
###########################################################################

def get_purged_hifiasm_or_original(wc):
    s = normalize_sample(wc)
    ragtag = outpath(f"Assemblies/{s}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.reforder.fasta")
    if s in RAGTAG_SAMPLES and os.path.exists(ragtag):
        return ragtag
    sorted_fa = get_sorted_purged_hifiasm_abs(wc)
    if os.path.exists(sorted_fa):
        return sorted_fa
    return get_asm_hifiasm_abs(wc)

def get_purged_lja_or_original(wc):
    s = normalize_sample(wc)
    ragtag = outpath(f"Assemblies/{s}/RAGTAG_LJA/purged-lja_ragtag/ragtag.scaffold.reforder.fasta")
    if s in RAGTAG_SAMPLES and os.path.exists(ragtag):
        return ragtag
    sorted_fa = get_sorted_purged_lja_abs(wc)
    if os.path.exists(sorted_fa):
        return sorted_fa
    return get_asm_lja_abs(wc)

###########################################################################
# Ploidy Plot Accessor
###########################################################################

def get_purged_fa(wc):
    if wc.asm == "hifiasm":
        return get_purged_hifiasm_or_original(wc)
    elif wc.asm == "lja":
        return get_purged_lja_or_original(wc)
    else:
        raise ValueError("Unknown ASM type in ploidy plot")
