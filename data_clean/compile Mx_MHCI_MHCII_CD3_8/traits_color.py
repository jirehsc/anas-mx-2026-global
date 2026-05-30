from Bio import SeqIO

INPUT_FASTA  = "trimmed mx_mhcI_mhcII_cd3_cd8.fasta"
OUTPUT_NEX   = "haplotype_mx_mhcI_mhcII_cd3_cd8_traits.NEX"

gene_map = {
    "KT340673.1":"Mx","NM_001310409.1":"Mx","MW246836.1":"Mx",
    "Z21549.1":"Mx","KR025554.1":"Mx",
    "AF393511.1":"MHC_I","KF319119.1":"MHC_I","KF319118.1":"MHC_I",
    "KF319120.1":"MHC_I","KF319121.1":"MHC_I",
    "HQ317493.1":"MHC_II","HQ317494.1":"MHC_II","AY905539.1":"MHC_II",
    "AY905541.1":"MHC_II","HM775322.1":"MHC_II","NM_001310815.1":"MHC_II",
    "NM_001310784.2":"MHC_II","JK310770.1":"MHC_II","JK311149.1":"MHC_II",
    "JK311170.1":"MHC_II","JK311229.1":"MHC_II","DR764108.1":"MHC_II",
    "DR764248.1":"MHC_II","DR764385.1":"MHC_II","DR764528.1":"MHC_II",
    "DR764750.1":"MHC_II","DR764988.1":"MHC_II","DR766023.1":"MHC_II",
    "DR766328.1":"MHC_II","DQ490139.1":"MHC_II","AF390589.1":"MHC_II",
    "FJ524847.1":"CD3","FJ524848.1":"CD3",
    "AF378373.1":"CD8","EF205151.1":"CD8","AY738733.1":"CD8",
    "FJ527828.1":"CD8","FJ528914.1":"CD8",
}

haplotypes = {
    "Hap_1": ["KT340673.1","NM_001310409.1","MW246836.1","Z21549.1","KR025554.1"],
    "Hap_2": ["AF393511.1","KF319119.1"],
    "Hap_3": ["KF319118.1"],
    "Hap_4": ["KF319120.1","KF319121.1"],
    "Hap_5": ["HQ317493.1","HQ317494.1","AY905539.1"],
    "Hap_6": ["AY905541.1"],
    "Hap_7": ["HM775322.1"],
    "Hap_8": ["NM_001310815.1","DR764108.1","DQ490139.1"],
    "Hap_9": ["NM_001310784.2","DR764248.1","DR766023.1"],
    "Hap_10":["JK310770.1","JK311149.1","JK311170.1","JK311229.1"],
    "Hap_11":["DR764385.1","DR764528.1"],
    "Hap_12":["DR764750.1","AF390589.1"],
    "Hap_13":["DR764988.1"],
    "Hap_14":["DR766328.1"],
    "Hap_15":["FJ524847.1","FJ524848.1"],
    "Hap_16":["AF378373.1","EF205151.1","AY738733.1","FJ528914.1"],
    "Hap_17":["FJ527828.1"],
}

records = list(SeqIO.parse(INPUT_FASTA, "fasta"))
n_seq = len(records)
seq_len = len(records[0].seq)

with open(OUTPUT_NEX, "w") as f:
    f.write("#NEXUS\n\n")
    f.write("BEGIN TAXA;\n")
    f.write(f"  Dimensions NTAX={n_seq};\n")
    f.write("  TAXLABELS\n")
    for rec in records:
        f.write(f"    {rec.id}\n")
    f.write("  ;\nEND;\n\n")

    f.write("BEGIN CHARACTERS;\n")
    f.write(f"  Dimensions NCHAR={seq_len};\n")
    f.write("  Format DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;\n")
    f.write("  Matrix\n")
    for rec in records:
        f.write(f"    {rec.id:<30} {str(rec.seq)}\n")
    f.write("  ;\nEND;\n\n")

    f.write("BEGIN TRAITS;\n")
    f.write("  Dimensions NTRAITS=1;\n")
    f.write("  Format labels=yes missing=? separator=Comma;\n")
    f.write("  TraitLabels Gene;\n")
    f.write("  Matrix\n")
    for hap, members in haplotypes.items():
        gene = gene_map[members[0]]
        for acc in members:
            f.write(f"    {acc:<30} {gene}\n")
    f.write("  ;\nEND;\n")

print(f"Done! Saved to: {OUTPUT_NEX}")
print(f"\nGene-Haplotype Summary:")
for hap, members in haplotypes.items():
    gene = gene_map[members[0]]
    print(f"  {hap} ({len(members)} seq): {gene}")
