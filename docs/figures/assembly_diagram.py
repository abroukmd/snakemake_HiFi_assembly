"""Compact assembly route map with original genomic artwork."""

from diagram_style import Diagram


def draw(path):
    d = Diagram(2080, 1280, "HiFi assembly workflow — assembly route map",
                "Full HiFi reads feed independent hifiasm and LJA routes, parallel OATK mitochondrial "
                "assembly, and optional MGA. Core assemblies undergo overlap-only purging, optional "
                "RagTag or its sorted-purged bypass, then mitochondrial integration. Raw and final QC "
                "checkpoints preserve assembler identity. MGA ends at its separate polished first-phase "
                "FASTA and native graphs. Required data are blue ellipses.")
    d.text(64, 65, "HiFi assembly workflow", 36, bold=True)
    d.text(64, 104, "Independent assemblies, refinement and QC", 20, "muted", max_width=580)

    d.route([(295, 470), (345, 470), (345, 368), (390, 368)], "input", width=6)
    d.route([(345, 470), (345, 608), (390, 608)], "input", width=6)
    d.route([(345, 608), (345, 868), (390, 868)], "mga", dash=True, width=5)
    d.route([(175, 420), (175, 216), (1000, 216)], "input", width=5)
    d.route([(935, 110), (960, 110), (960, 216), (1000, 216)], "input", width=5)
    d.route([(175, 520), (175, 770)], "input", width=5)

    # Reference is required by Chromeister even when RagTag is disabled.
    d.route([(1780, 155), (1780, 280), (1100, 280), (1100, 330)], "input", width=3.8)
    d.route([(1238, 280), (1238, 520), (1100, 520), (1100, 570)], "input", width=3.8)
    d.route([(1910, 110), (2020, 110), (2020, 1022), (1850, 1022), (1850, 977)], "input", width=3.8)
    d.route([(1790, 1090), (1790, 1045), (1740, 1045), (1740, 977)], "input", width=3.8)

    # One shared organellar branch supplies both core integration steps.
    d.route([(1200, 216), (1300, 216)], "mito", width=6)
    d.route([(1417, 252), (1417, 330)], "mito", width=5)
    d.route([(1535, 216), (1590, 216), (1590, 525), (1417, 525), (1417, 570)], "mito", width=5)

    for y, color in [(368, "hif"), (608, "lja")]:
        d.route([(600, y), (700, y)], color, width=6)
        d.route([(915, y), (1000, y)], color, width=6)
        d.route([(1200, y), (1300, y)], color, width=6)
        d.route([(1535, y), (1650, y)], color, width=6)
        d.route([(940, y), (950, y), (950, y + 82), (1417, y + 82), (1417, y + 38)],
                "bypass", dash=True, width=4.5)
        d.label(1170, y + 74, "ragtag: N", color="bypass", size=17)
        d.label(650, y - 21, "raw", size=17)
        d.label(960, y - 21, "primary", size=16)
        d.route([(650, y), (650, y + 50)], "qc", width=4)
        d.checkpoint(650, y + 66)
        d.route([(1910, y), (1980, y), (1980, 921), (1970, 921)], "qc", width=4)
    d.route([(600, 868), (700, 868)], "mga", width=6)

    # Nodes are layered after the routes, keeping their interiors clean.
    d.input(175, 470, 120, 50, ["PacBio HiFi", "reads"], size=22)
    d.input(810, 110, 125, 45, ["OATK database", "organellar models"], size=19)
    d.input(1780, 110, 130, 45, ["External reference", "genome"], size=19)
    d.input(1790, 1135, 130, 45, ["Illumina", "R1 / R2"], size=22)
    d.node(1000, 180, 200, 72, "OATK", "mito assembly", size=25)
    d.node(1300, 180, 235, 72, "Mito contigs", "OATK FASTA", "mito", output=True, size=23)

    for y, color, name in [(330, "hif", "hifiasm"), (570, "lja", "LJA")]:
        d.node(390, y, 210, 76, name, size=27)
        d.node(700, y, 215, 76, "purge_dups", "overlap only", size=23)
        d.node(1000, y, 200, 76, "RagTag", "optional scaffolding", "bypass", optional=True, size=25)
        d.node(1300, y, 235, 76, "Mito integration", "Coverage + identity", size=23)
        d.node(1650, y, 260, 76, "Final candidate", f"{name} FASTA", color, output=True, size=24)
        d.text(807, y + 127, "Primary + hap FASTAs retained", 17, "muted", anchor="middle")

    d.node(55, 770, 240, 100, "Read profiles", ["NanoPlot", "TeloSearchLR"], size=23)
    d.node(390, 830, 210, 76, "MGA", "optional first phase", "mga", optional=True, size=27)
    d.node(700, 818, 370, 100, "MGA consensus", "Polished FASTA + native graphs", "mga", output=True, size=25)
    d.text(1100, 856, "Ends after native polishing.", 21, "mga", bold=True)
    d.text(1100, 887, "Separate from core refinement and QC.", 18, "muted")
    d.node(1650, 865, 320, 112, "Assembly QC", ["BUSCO · QUAST", "Merqury · Chromeister"], "qc", output=True, size=26)
    d.text(1810, 836, "Raw + final checkpoints", 19, "muted", anchor="middle")

    for x, text in [(495, "ASSEMBLE"), (807, "PURGE"), (1100, "SCAFFOLD"),
                    (1417, "INTEGRATE MITO"), (1780, "CANDIDATE FASTA")]:
        d.text(x, 307, text, 14, "muted", bold=True, anchor="middle")
    d.legend(60, 941, [
        ("hif", "hifiasm route", "Processed independently from LJA"),
        ("lja", "LJA route", "Processed independently from hifiasm"),
        ("mga", "Optional MGA branch", "run_mga: Y in sample sheet", True),
        ("mito", "Mitochondrial contigs", "OATK supplies both core branches"),
        ("bypass", "Skip reference scaffolding", "ragtag: N; continue from sorted purged contigs", True),
    ], gap=51)
    d.motif(1260, 1000, scale=1.1)
    d.checkpoint(1703, 1235)
    d.text(1730, 1241, "Assembly QC checkpoint", 16, "muted")
    d.text(64, 1261, "Blue ellipses: required data. Dashed nodes/routes: optional modules. See the QC map for model, final-mapping and optional-analysis scopes.",
           16, "muted", max_width=1550)
    d.save(path)
