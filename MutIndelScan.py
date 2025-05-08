import os
import csv
from Bio import SeqIO
import pandas as pd

def find_indels_and_mutations(seq1, seq2):
    """
    Find indels and mutations between two aligned sequences.
    """
    insertions = []
    deletions = []
    mutations = []
    
    gap_char = '-'

    assert len(seq1) == len(seq2), "Sequences are not aligned or of different lengths."
    
    in_insertion = False
    in_deletion = False
    insertion_start = None
    deletion_start = None
    
    for i, (a, b) in enumerate(zip(seq1.upper(), seq2.upper())):
        if a == gap_char and b != gap_char:  
            if not in_insertion:
                insertion_start = i
                in_insertion = True
            in_deletion = False
        elif b == gap_char and a != gap_char:  
            if not in_deletion:
                deletion_start = i
                in_deletion = True
            in_insertion = False
        elif a != b:  
            mutations.append((i, a, b))
            in_insertion = False
            in_deletion = False
        else:
            if in_insertion:
                insertions.append((insertion_start, i))
                in_insertion = False
            if in_deletion:
                deletions.append((deletion_start, i))
                in_deletion = False
    
    if in_insertion:
        insertions.append((insertion_start, len(seq1)))
    if in_deletion:
        deletions.append((deletion_start, len(seq1)))
    
    return insertions, deletions, mutations

def process_fasta(input_fasta, output_file):
    print(f"Processing {input_fasta}...")
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    if len(sequences) < 2:
        print(f"Not enough sequences to compare in {input_fasta}.")
        return
    
    ref_seq = sequences[0]
    file_name = os.path.splitext(os.path.basename(input_fasta))[0]
    
    rows = []
    for seq in sequences[1:]:
        insertions, deletions, mutations = find_indels_and_mutations(str(ref_seq.seq), str(seq.seq))
        
        for start, end in insertions:
            rows.append({
                'File_Name': file_name,
                'Ref_ID': ref_seq.id,
                'Seq_ID': seq.id,
                'Type': 'Insertion',
                'Start': start + 1,
                'End': end,
                'Indel_Length': end - start - str(seq.seq[start:end]).count('-'),  # Exclude gaps from the length
                'Sequence': str(seq.seq[start:end]).upper(),
                'Ref_Base': '',
                'Seq_Base': ''
            })
        
        for start, end in deletions:
            rows.append({
                'File_Name': file_name,
                'Ref_ID': ref_seq.id,
                'Seq_ID': seq.id,
                'Type': 'Deletion',
                'Start': start + 1,
                'End': end,
                'Indel_Length': end - start - str(ref_seq.seq[start:end]).count('-'),  # Exclude gaps from the length
                'Sequence': str(ref_seq.seq[start:end]).upper(),
                'Ref_Base': '',
                'Seq_Base': ''
            })
        
        for pos, ref_base, seq_base in mutations:
            rows.append({
                'File_Name': file_name,
                'Ref_ID': ref_seq.id,
                'Seq_ID': seq.id,
                'Type': 'Mutation',
                'Start': pos + 1,
                'End': pos + 1,
                'Indel_Length': '',
                'Sequence': '',
                'Ref_Base': ref_base.upper(),
                'Seq_Base': seq_base.upper()
            })
    
    with open(output_file, "a", newline='') as csvfile:
        fieldnames = ['File_Name', 'Ref_ID', 'Seq_ID', 'Type', 'Start', 'End', 'Indel_Length', 'Sequence', 'Ref_Base', 'Seq_Base']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        if os.stat(output_file).st_size == 0:
            writer.writeheader()
        writer.writerows(rows)
    
    print(f"Finished processing {input_fasta}. Output appended to {output_file}.")

def calculate_gc(sequence):
    gc_count = sequence.upper().count("G") + sequence.upper().count("C")
    gc_content = (gc_count / len(sequence.replace('-', ''))) * 100  # Exclude gaps from the length
    return round(gc_content, 2)

def main():
    print("Starting batch processing...")
    output_file = "Indels&Mutations.csv"
    
    if os.path.exists(output_file):
        os.remove(output_file)

    fasta_files = [file for file in os.listdir() if file.endswith('.fasta')]
    
    for fasta_file in fasta_files:
        process_fasta(fasta_file, output_file)

    gc_info = {}
    length_info = {}

    for fasta_file in fasta_files:
        file_name = os.path.splitext(fasta_file)[0]
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                sequence_id = record.id.split('_')[0]  
                sequence = str(record.seq)
                gc_content = calculate_gc(sequence)
                seq_length = len(sequence.replace('-', ''))
                if file_name not in gc_info:
                    gc_info[file_name] = {}
                    length_info[file_name] = {}
                gc_info[file_name][sequence_id] = gc_content
                length_info[file_name][sequence_id] = seq_length

    gc_df = pd.DataFrame(gc_info).transpose()
    gc_df.index.name = 'file_name'
    
    length_df = pd.DataFrame(length_info).transpose()
    length_df.index.name = 'file_name'

    with pd.ExcelWriter("GCcontent%Size.xlsx") as writer:
        gc_df.to_excel(writer, sheet_name='GC Content')
        length_df.to_excel(writer, sheet_name='Sequence Length')

    print("GC content and sequence length information saved to GCcontent%Size.xlsx")

    print("Batch processing completed.")
    print(f"Combined output saved to {output_file}")

if __name__ == "__main__":
    main()
