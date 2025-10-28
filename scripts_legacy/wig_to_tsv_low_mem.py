import argparse
import numpy as np
import re

parser = argparse.ArgumentParser(description="Process SAM file and reference FASTA file to produce a TSV file with enhanced single molecule information.")
parser.add_argument("-i", help="Input SAM file path", required=True)
parser.add_argument("-r", help="Reference FASTA file path", required=True)
parser.add_argument("-o", help="Output TSV file path", required=True)

args = parser.parse_args()

in_file = args.i
ref_file = args.r
out_file = args.o

def enhance_wig(data_path, ref_path, out_file):
	'''
	Main function to process SAM file and reference FASTA file to produce a TSV file with single molecule information.
	'''
	# Read in reference fasta
	ref_dict = read_reference_fasta(ref_path)
	# Read SAM file and write TSV output
	read_sam_file(data_path, ref_dict, out_file)

def read_reference_fasta(ref_path):
	# Read in reference fasta
	try:
		with open(ref_path) as f:
			lines = f.readlines()
	except OSError as e:
		print(f"Error reading file {ref_path}: {e}")
		return {}
	
	ref_dict = {}
	name = ""
	seq = ""
	for line in lines:
		if line.startswith(">"):
			if name: 
				ref_dict[name] = seq  # Update the dictionary once per sequence
			name = line.split()[0].strip(">")
			seq = ""  # Reset sequence for the new entry
		elif line.strip():
			seq += line.strip()  # Accumulate sequence lines
	if name:
		ref_dict[name] = seq  # Update the dictionary for the last sequence
	return ref_dict

def read_sam_file(data_path, ref_dict, out_file):
	try:
		with open(data_path) as f:
			sam_lines = f.readlines()
	except OSError as e:
		print(f"Error reading file {data_path}: {e}")
		return

	# Initialize output file with header
	with open(out_file, "w") as f:
		f.write("read_id\tgene\tposition\tbase\tdeletion\tinsertion\tmutation\tA\tC\tG\tT\tN\n")

	with open(out_file, "a") as f:
		for line in sam_lines:
			if line.startswith("@"):
				continue

			parts = line.strip().split("\t")
			read_id = parts[0]
			gene = parts[2]
			pos = int(parts[3])
			seq = parts[9]
			
			cigartuples = [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', parts[5])]
			gene_seq = ref_dict.get(gene, "")

			current_pos = pos
			read_index = 0

			for length, operation in cigartuples:
				if operation == "M":  # Match/Mismatch
					for i in range(length):
						base = seq[read_index]
						ref_base = gene_seq[current_pos-1] if current_pos-1 < len(gene_seq) else ""
						mutation = 1 if base != ref_base else 0
						f.write(f"{read_id}\t{gene}\t{current_pos}\t{base}\t0\t0\t{mutation}\t0\t0\t0\t0\t0\n")
						current_pos += 1
						read_index += 1
				elif operation == "I":  # Insertion
					for i in range(length):
						base = seq[read_index]
						f.write(f"{read_id}\t{gene}\t{current_pos}\t{base}\t0\t1\t0\t0\t0\t0\t0\t0\n")
						read_index += 1
				elif operation == "D":  # Deletion
					for i in range(length):
						f.write(f"{read_id}\t{gene}\t{current_pos}\t-\t1\t0\t0\t0\t0\t0\t0\t0\n")
						current_pos += 1
				else:
					read_index += length if operation in ("H", "S") else 0  # Hard and soft clip increments read index

if __name__ == "__main__":
	enhance_wig(in_file, ref_file, out_file)


