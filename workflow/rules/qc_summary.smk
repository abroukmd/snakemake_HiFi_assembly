from types import SimpleNamespace


def qc_summary_records():
    records = []
    for sample in SAMPLES:
        wc = SimpleNamespace(sample=sample)
        for asm in enabled_assemblers(wc):
            for stage in ("raw", "final"):
                raw = stage == "raw"
                assembly = ((get_asm_hifiasm_abs if asm == "hifiasm" else get_asm_lja_abs)(wc) if raw
                            else (get_final_hifiasm_abs if asm == "hifiasm" else get_final_lja_abs)(wc))
                busco_label = f"{sample}-{asm}" if raw else f"{sample}-purged-{asm}"
                busco_branch = f"BUSCO-{asm}" if raw else f"BUSCO-purged-{asm}"
                records.append({
                    "sample": sample, "assembler": asm, "stage": stage, "assembly": assembly,
                    "busco": outpath(f"Assemblies/{sample}/QC/{busco_branch}/{busco_label}"),
                    "merqury": outpath(f"Assemblies/{sample}/QC/merqury/{stage}/{asm}"),
                    "mapping": None if raw else outpath(f"Assemblies/{sample}/QC/Mapping/purged-{asm}.flagstat.txt"),
                    "tidk": outpath(f"Assemblies/{sample}/QC/tidk/{asm}/assembly_summary.tsv") if not raw and RUN_TIDK else None,
                    "heterozygosity": outpath(f"Assemblies/{sample}/QC/heterozygosity/{asm}/summary.json") if not raw and RUN_HETEROZYGOSITY else None,
                })
    return records


QC_SUMMARY_RECORDS = qc_summary_records()


rule qc_summary:
    input:
        source_script=workflow.source_path("../scripts/qc_summary.py"),
        assemblies=[row["assembly"] for row in QC_SUMMARY_RECORDS],
        busco=[row["busco"] for row in QC_SUMMARY_RECORDS],
        qv=[row["merqury"] + "/merqury.qv" for row in QC_SUMMARY_RECORDS],
        completeness=[row["merqury"] + "/completeness.stats" for row in QC_SUMMARY_RECORDS],
        merqury_directories=[row["merqury"] for row in QC_SUMMARY_RECORDS],
        mapping=[row["mapping"] for row in QC_SUMMARY_RECORDS if row["mapping"]],
        tidk=[row["tidk"] for row in QC_SUMMARY_RECORDS if row["tidk"]],
        tidk_config=[str(Path(row["tidk"]).parent / "run_config.json") for row in QC_SUMMARY_RECORDS if row["tidk"]],
        tidk_receipts=[str(Path(row["tidk"]).parent / "complete.json") for row in QC_SUMMARY_RECORDS if row["tidk"]],
        heterozygosity=[row["heterozygosity"] for row in QC_SUMMARY_RECORDS if row["heterozygosity"]],
        heterozygosity_receipts=[str(Path(row["heterozygosity"]).parent / "receipt.json") for row in QC_SUMMARY_RECORDS if row["heterozygosity"]]
    output:
        table=outpath("Assemblies/qc_summary/assembly_qc.tsv"),
        provenance=outpath("Assemblies/qc_summary/provenance.json"),
        pdf=outpath("Assemblies/qc_summary/assembly_qc.pdf"),
        svg=outpath("Assemblies/qc_summary/assembly_qc.svg"),
        png=outpath("Assemblies/qc_summary/assembly_qc.png")
    params:
        records=QC_SUMMARY_RECORDS
    log:
        outpath("Assemblies/qc_summary/qc_summary.log")
    threads: 1
    resources:
        mem_mb=4000
    conda:
        "../envs/qc_summary_env.yaml"
    script:
        "../scripts/qc_summary.py"
