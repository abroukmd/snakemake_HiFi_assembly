"""Draw the release QC schematic, preserving exact stage and input scopes."""

from diagram_style import C, Diagram


def draw(path):
    d = Diagram(
        2200, 1600,
        "HiFi assembly workflow: quality-control evidence",
        "Connected quality-control schematic for the selected hifiasm and LJA branches. "
        "HiFi reads feed NanoPlot and TeloSearchLR; pooled HiFi and both Illumina mates feed "
        "KMC and exploratory GenomeScope2 fits. Raw and final core assemblies feed BUSCO, "
        "reference-free QUAST, Illumina-k-mer Merqury, and reference-based Chromeister. "
        "Additional haplotig Chromeister comparisons are manual targets. "
        "Final assemblies and Illumina reads produce matching BAMs; nQuire selects one "
        "preferred final BAM per sample. Optional TIDK and depth-callable diploid-model "
        "heterozygosity have separate inputs and outputs. Heterozygosity also uses OATK mito "
        "FASTA for sequence-verified mitochondrial exclusion. Barrnap is restricted to pre-mito "
        "hifiasm. Aggregate QC has defined metric sources; MultiQC is an independent dashboard "
        "of supported native reports. MGA is outside these core QC branches."
    )
    d.text(70, 62, "QUALITY CONTROL", 16, "input", True)
    d.text(70, 112, "Follow the evidence", 42, bold=True)
    d.text(70, 148, "Exact inputs, explicit assembly stages and independent complementary analyses.",
           22, "muted")

    # Light outlines provide local reading order without replacing the graph.
    d.group(400, 175, 1255, 365, "Read profiles and exploratory models", "input")
    d.group(400, 560, 1255, 355, "Each selected hifiasm / LJA assembly · raw and final", "qc")
    d.group(400, 925, 1255, 485, "Final-assembly and alignment evidence", "qc")
    d.group(1720, 175, 420, 720, "Optional report summaries", "mga")
    d.group(1720, 925, 420, 485, "Stage-specific annotation", "mito")

    # HiFi read branches. The KMC junction pools HiFi with BOTH Illumina mates.
    d.route([(315, 280), (370, 280), (370, 255), (480, 255)], "input")
    d.route([(370, 280), (370, 325), (900, 325), (900, 290)], "input")
    d.route([(350, 280), (350, 350), (440, 350), (440, 410), (480, 410)], "input")
    d.route([(315, 460), (440, 460), (440, 410)], "input", arrow=False)
    d.route([(710, 410), (785, 410)], "hif")
    d.route([(315, 460), (340, 460), (340, 510), (1035, 510),
             (1035, 450), (1090, 450)], "input")
    # The Illumina-only k-mer database is a separate Merqury dependency.
    d.route([(1205, 485), (1205, 605)], "input")
    d.route([(1205, 685), (1205, 730)], "input")

    d.input(180, 280, 135, 49, ["PacBio HiFi", "reads"], size=21)
    d.input(180, 460, 135, 49, ["Illumina", "R1 + R2"], size=21)
    d.node(480, 220, 230, 70, "NanoPlot", "Read length / yield", size=25)
    d.node(785, 220, 230, 70, "TeloSearchLR", "Read-repeat discovery", size=24)
    d.node(480, 375, 230, 70, "KMC", "Pooled read k-mers", size=25)
    d.node(785, 375, 230, 70, "GenomeScope2", "Exploratory fits", size=24)
    d.node(1090, 415, 230, 70, "meryl", "Illumina only", size=25)
    d.label(570, 487, "HiFi + both Illumina mates", "muted", 18)

    # Raw/final FASTAs feed all four tools through one stage-aware bus.
    d.route([(315, 755), (370, 755), (370, 870), (1495, 870)], "hif", arrow=False)
    for cx in [595, 900, 1205, 1495]:
        d.route([(cx, 870), (cx, 805)], "hif")
    d.route([(595, 685), (595, 730)], "input")
    d.route([(1495, 685), (1495, 730)], "input")
    d.input(180, 755, 135, 59, ["Raw + final", "FASTA", "hifiasm / LJA"], size=20)
    d.input(595, 645, 130, 40, ["BUSCO", "lineage DB"], size=19)
    d.input(1205, 645, 140, 40, ["Illumina", "read k-mers"], size=19)
    d.input(1495, 645, 130, 40, ["External", "reference"], size=19)
    d.node(480, 730, 230, 75, "BUSCO", "Gene completeness", size=25)
    d.node(785, 730, 230, 75, "QUAST", "Assembly statistics", size=25)
    d.node(1090, 730, 230, 75, "Merqury", "QV + completeness", size=25)
    d.node(1380, 730, 230, 75, "Chromeister", "Reference dotplots", size=24)
    d.label(900, 842, "No external reference", "muted", 17)
    d.label(1010, 899, "Keep sample × assembler × stage identities", "muted", 19)
    d.label(1495, 899, "Haplotig plots: manual", "muted", 16)

    # Matching final FASTAs and Illumina reads generate one BAM per core branch.
    d.route([(315, 1070), (480, 1070)], "input")
    d.route([(315, 1205), (340, 1205), (340, 1145), (595, 1145),
             (595, 1105)], "input")
    d.route([(710, 1070), (775, 1070)], "hif")
    d.route([(1005, 1070), (1040, 1070), (1040, 1000), (1090, 1000)], "hif")
    d.route([(1005, 1070), (1040, 1070), (1040, 1170), (1090, 1170)],
            "mga", dash=True)
    # Heterozygosity also consumes the exact assembly used to map the BAM.
    d.route([(350, 1070), (370, 1070), (370, 1350), (1240, 1350),
             (1240, 1215)], "mga", dash=True)
    d.route([(370, 1070), (370, 1270), (480, 1270)], "mga", dash=True)
    d.route([(1515, 1230), (1515, 1170), (1390, 1170)], "input")
    d.input(180, 1070, 135, 50, ["Final FASTA", "hifiasm / LJA"], size=20)
    d.input(180, 1205, 135, 43, ["Illumina", "R1 + R2"], size=20)
    d.node(480, 1035, 230, 70, "BWA / SAMtools", "Map to each final FASTA", size=23)
    d.input(890, 1070, 115, 47, ["Matching", "final BAM"], size=20)
    d.node(1090, 965, 300, 70, "nQuire", "One preferred BAM / sample", size=25)
    d.label(1240, 1062, "hifiasm if selected; otherwise LJA", "muted", 17)
    d.node(1090, 1125, 300, 90, "Heterozygosity", "Depth-callable diploid model",
           color="mga", optional=True, size=25)
    d.node(480, 1235, 230, 70, "TIDK", "Exact repeats at ends",
           color="mga", optional=True, size=25)
    d.input(1515, 1275, 115, 45, ["OATK mito", "FASTA"], size=19)
    d.text(1515, 1350, "Verify mito exclusion", 16, "muted", anchor="middle")
    d.text(760, 1278, "Candidate arrays / end coordinates", 19, "muted")
    d.label(1030, 1376, "Exact final FASTA + matching BAM", "muted", 18)
    d.text(890, 1128, "Proper-pairing evidence", 17, "muted", anchor="middle")

    # These two summaries have different source sets and no implied merge.
    d.route([(1930, 320), (1930, 355)], "mga", dash=True)
    d.route([(1930, 667), (1930, 720)], "mga", dash=True)
    d.input(1930, 270, 180, 50, ["Supported native", "QC reports"], size=21)
    d.node(1790, 355, 280, 80, "MultiQC", "Independent HTML dashboard",
           color="mga", optional=True, size=25)
    d.text(1930, 478, "Supported modules only", 18, "muted", anchor="middle")
    d.input(1930, 600, 190, 67,
            ["FASTA · BUSCO · Merqury", "final mapping", "± enabled TIDK / heterozygosity"], size=18)
    d.node(1790, 720, 280, 80, "Aggregate QC", "TSV · figures · provenance",
           color="mga", optional=True, size=25)
    d.text(1930, 842, "sample × assembler × stage", 18, "muted", anchor="middle")

    # Barrnap is deliberately separated from final-FASTA QC.
    d.route([(1930, 1135), (1930, 1205)], "mito")
    d.input(1930, 1080, 180, 55, ["Pre-mito hifiasm", "purged / RagTag stage"], size=20)
    d.node(1790, 1205, 280, 80, "Barrnap", "rRNA GFF3", size=25)
    d.text(1930, 1332, "hifiasm branch only", 19, "muted", anchor="middle")

    # Compact line legend; scientific caveats remain visible below it.
    for x, color, label, dash in [
        (70, "input", "Input / reference", False),
        (590, "hif", "Core QC / mapping", False),
        (1110, "mga", "Optional module / summary", True),
        (1715, "mito", "Stage-specific annotation", False),
    ]:
        d.route([(x, 1470), (x + 66, 1470)], color, dash=dash, arrow=False)
        d.text(x + 83, 1477, label, 20)
    d.text(1100, 1530, "Model fits do not establish ploidy; repeat candidates do not establish T2T.",
           20, "muted", anchor="middle")
    d.text(1100, 1565, "MGA remains a separate first-phase candidate outside these hifiasm / LJA QC branches.",
           20, "muted", anchor="middle")
    d.save(path)
