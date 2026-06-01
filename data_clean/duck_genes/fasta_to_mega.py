import os
import time
import io
import ssl
import urllib.request
import urllib.parse
from Bio import SeqIO
from Bio.Seq import Seq

input_file = "Mx.fasta"
output_file = "Mx_cleaned.meg"
raw_cds_records = []

print("🔄 Starting the automated pipeline...")

# 1. DELETE OLD FILE TO PREVENT STALE CACHED DATA
if os.path.exists(output_file):
    try:
        os.remove(output_file)
        print("🗑️ Deleted old output file to guarantee fresh results.")
    except PermissionError:
        print("❌ ERROR: 'Mx_cleaned.meg' is locked by another program!")
        print("👉 Please CLOSE DnaSP completely, then run this script again.")
        exit()

def find_true_cds_by_atg(sequence_str):
    ungapped = sequence_str.replace("-", "").upper()
    start_indices = [i for i in range(len(ungapped) - 2) if ungapped[i:i+3] == "ATG"]
    
    if not start_indices:
        return ungapped
        
    best_cds = ""
    max_protein_len = -1
    
    for start in start_indices:
        test_seq = Seq(ungapped[start:])
        trim_len = (len(test_seq) // 3) * 3
        test_seq = test_seq[:trim_len]
        protein = str(test_seq.translate())
        
        if "*" in protein:
            stop_idx = protein.index("*")
            # Slice BEFORE the stop codon so DnaSP never sees it and stops complaining
            cds_dna = str(test_seq[:stop_idx * 3])
            protein_len = stop_idx
        else:
            cds_dna = str(test_seq)
            protein_len = len(protein)
            
        if protein_len > max_protein_len:
            max_protein_len = protein_len
            best_cds = cds_dna
            
    return best_cds

if not os.path.exists(input_file):
    print(f"❌ Error: Could not find your input file '{input_file}' in this folder.")
    exit()

print("Step 1: Stripping UTRs and isolating true coding frames...")
for record in SeqIO.parse(input_file, "fasta"):
    clean_id = record.id.replace('"', '').replace(" ", "_").replace("'", "")
    cleaned_dna = find_true_cds_by_atg(str(record.seq))
    raw_cds_records.append(f">{clean_id}\n{cleaned_dna}")

fasta_payload = "\n".join(raw_cds_records)

print("Step 2: Submitting sequences to EMBL-EBI Clustal Omega for alignment...")
submit_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
data = {
    'email': 'duck_researcher@gmail.com',
    'title': 'duck_mx_alignment',
    'sequence': fasta_payload
}

encoded_data = urllib.parse.urlencode(data).encode('utf-8')
req = urllib.request.Request(submit_url, data=encoded_data)

# Create an unverified context to totally bypass Windows SSL certificate blocks
unverified_context = ssl._create_unverified_context()

try:
    with urllib.request.urlopen(req, context=unverified_context) as response:
        job_id = response.read().decode('utf-8').strip()
    print(f"   🚀 Job accepted by server! (ID: {job_id})")
except Exception as e:
    print(f"❌ Network Error: Failed to connect to alignment server. Details: {e}")
    exit()

status_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"
while True:
    try:
        with urllib.request.urlopen(status_url, context=unverified_context) as response:
            status = response.read().decode('utf-8').strip()
    except:
        status = "RUNNING"
        
    print(f"   Alignment status: {status}")
    if status in ["FINISHED", "FAILURE", "NOT_FOUND"]:
        break
    time.sleep(3)

if status != "FINISHED":
    print("❌ Alignment failed on the remote server side.")
    exit()

print("Step 3: Downloading perfectly aligned matrix...")
result_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-fasta"
with urllib.request.urlopen(result_url, context=unverified_context) as response:
    aligned_fasta = response.read().decode('utf-8')

print("Step 4: Compiling clean alignment matrix into MEGA format...")
with open(output_file, "w") as f:
    f.write("#mega\n")
    f.write("!Title Duck_Genes_Aligned;\n")
    f.write("!Format DataType=DNA;\n\n")
    
    for record in SeqIO.parse(io.StringIO(aligned_fasta), "fasta"):
        f.write(f"#{record.id}\n{str(record.seq)}\n")

print(f"\n🎉 SUCCESS! Fresh, aligned data saved directly to: {output_file}")
