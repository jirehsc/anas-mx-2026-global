import re
with open("haplotype mx_mhcI_mhcII_cd3_cd8.nex", "r") as f:
    content = f.read()

GENE_TAGS = {
    "Mx":     ["KT340673","NM_001310409","MW246836","Z21549","KR025554"],
    "MHC_I":  ["AF393511","KF319118","KF319120","KF319121","KF319119"],
    "MHC_II": ["HQ317493","HQ317494","AY905541","HM775322","AY905539",
               "NM_001310815","NM_001310784","JK310770","JK311149","JK311170",
               "JK311229","DR764108","DR764248","DR764385","DR764528",
               "DR764750","DR764988","DR766023","DR766328","DQ490139","AF390589"],
    "CD3":    ["FJ524847","FJ524848"],
    "CD8":    ["AF378373","EF205151","AY738733","FJ527828","FJ528914"],
}

def get_gene(seq_id):
    for gene, ids in GENE_TAGS.items():
        if any(tag in seq_id for tag in ids):
            return gene
    return "Unknown"

hap_gene_map = {}
for match in re.finditer(r'\[Hap_(\d+):\s+\d+\s+(.+?)\]', content):
    hap_num = f"Hap_{match.group(1)}"
    seq_ids = match.group(2).strip().split()
    if any("." in s for s in seq_ids):
        gene = get_gene(seq_ids[0])
        hap_gene_map[hap_num] = gene
        print(f"  {hap_num}: -> {gene}")

genes_list = list(GENE_TAGS.keys())
traits_block = "\nBEGIN TRAITS;\n"
traits_block += f"    Dimensions NTRAITS={len(genes_list)};\n"
traits_block += "    Format labels=yes missing=? separator=Comma;\n"
traits_block += f"    TraitLabels {' '.join(genes_list)};\n"
traits_block += "    Matrix\n"
for hap, gene in hap_gene_map.items():
    values = ["1" if g == gene else "0" for g in genes_list]
    traits_block += f"    {hap}    {','.join(values)}\n"
traits_block += "    ;\nEND;\n"

with open("haplotype mx_mhcI_mhcII_cd3_cd8_traits.NEX", "w", encoding="utf-8") as f:
    f.write(content + traits_block)

print("\nSuccess! Saved to: haplotype mx_mhcI_mhcII_cd3_cd8_traits.NEX")
