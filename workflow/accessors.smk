import os

###########################################################################
# Unified Accessors for Assemblies, Purged Assemblies, Final Assemblies,
# and Fastq Files
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


# Purged assemblies before mito integration

def get_purged_hifiasm_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"Assemblies/{s}/PURGE_DUPS/purged.hifiasm.fa")


def get_purged_lja_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"Assemblies/{s}/PURGE_DUPS_LJA/purged.lja.fa")


# Sorted purged assemblies before mito integration

def get_sorted_purged_hifiasm_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"Assemblies/{s}/PURGE_DUPS/purged-hifiasm.sorted.fa")


def get_sorted_purged_lja_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"Assemblies/{s}/PURGE_DUPS_LJA/purged-lja.sorted.fa")


# RagTag-renamed final assemblies before mito integration

def get_ragtag_final_hifiasm_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"Assemblies/{s}/RAGTAG/FINAL_{s}/{s}.purged.final.fasta")


def get_ragtag_final_lja_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"Assemblies/{s}/RAGTAG_LJA/FINAL_{s}/{s}.purged.final.fasta")


###########################################################################
# Purged assembly selection before mito integration
# (for rules that should still run on the pre-final assembly state)
###########################################################################

def get_purged_hifiasm_or_original(wc):
    s = normalize_sample(wc)
    if s in RAGTAG_SAMPLES:
        return get_ragtag_final_hifiasm_abs(wc)
    return get_sorted_purged_hifiasm_abs(wc)


def get_purged_lja_or_original(wc):
    s = normalize_sample(wc)
    if s in RAGTAG_SAMPLES:
        return get_ragtag_final_lja_abs(wc)
    return get_sorted_purged_lja_abs(wc)


###########################################################################
# Final assemblies after mito integration
# These are the assemblies that downstream QC should evaluate.
###########################################################################

def get_final_hifiasm_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"{s}_FINAL_ASSEMBLIES/hifiasm/{s}.purged.final.fasta")


def get_final_lja_abs(wc):
    s = normalize_sample(wc)
    return outpath(f"{s}_FINAL_ASSEMBLIES/lja/{s}.purged.final.fasta")


###########################################################################
# Ploidy Plot / generic accessor helper
###########################################################################

def get_purged_fa(wc):
    if wc.asm == "hifiasm":
        return get_final_hifiasm_abs(wc)
    elif wc.asm == "lja":
        return get_final_lja_abs(wc)
    else:
        raise ValueError("Unknown ASM type in ploidy plot")
