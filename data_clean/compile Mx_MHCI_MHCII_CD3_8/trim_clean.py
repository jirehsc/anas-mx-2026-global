from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

INPUT_FILE    = "aligned mx_mhcI_mhcII_cd3_cd8.FAS"
OUTPUT_FILE   = "mx_mhcI_mhcII_cd3_cd8_trimmed.fasta"
GAP_THRESHOLD = 0.50

print("Reading aligned sequences...")
alignment = AlignIO.read(INPUT_FILE, "fasta")
n_seq   = len(alignment)
aln_len = alignment.get_alignment_length()

print(f"  Sequences  : {n_seq}")
print(f"  Alignment  : {aln_len} bp\n")

print("Analyzing gaps...")
cols_to_keep = [
    i for i in range(aln_len)
    if alignment[:, i].count("-") / n_seq <= GAP_THRESHOLD
]

if not cols_to_keep:
    print("ERROR: No columns survived filter.")
    raise SystemExit(1)

trimmed_records = []
for rec in alignment:
    seq_str = "".join(rec.seq[i] for i in cols_to_keep)
    trimmed_records.append(
        SeqRecord(Seq(seq_str), id=rec.id, description=rec.description)
    )

trimmed_len = len(cols_to_keep)
SeqIO.write(trimmed_records, OUTPUT_FILE, "fasta")

removed  = aln_len - trimmed_len
pct_kept = trimmed_len / aln_len * 100

print(f"Trimming Results:")
print(f"  Original alignment : {aln_len} bp")
print(f"  Trimmed alignment  : {trimmed_len} bp ({pct_kept:.1f}%)")
print(f"  Columns removed    : {removed} bp ({100-pct_kept:.1f}%)")
print(f"\n√ Results saved to: {OUTPUT_FILE}")
