from Bio import SeqIO
from statistics import median
from Bio.Seq import Seq
import argparse
import sys

def main(argv, out):
# Write help message and set arguments
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', help='File to process (required)')
	parser.add_argument('--isolate', help='Sample ID (required)')
	parser.add_argument('--gap', help="Minimum gap length to be considered a scaffold (optional) [2]", default=2)
	parser.add_argument('--output', help="Output file name (optional) ['asmstats.txt']", default="asmstats")
	args = parser.parse_args(argv)
# Run each function and select relevant outputs
	count, s_n50, gap_sum, total_asm, s_n90, gc_cont, gap_perc = calculate_scaffold_stats(args.fasta)
	c_count, c_n50, gap_count, c_n90 = calculate_contig_stats(args.fasta, args.gap)
# Write outputs to file
	write_stats_to_file(args.isolate, args.output+".txt", count, s_n50, c_count, c_n50, gap_count, gap_sum, total_asm, s_n90, c_n90, gc_cont, gap_perc)

# Calculate common scaffold statistics
def calculate_scaffold_stats(fasta_file):
	records = list(SeqIO.parse(fasta_file, "fasta"))
	count = len(records) # Number of scaffolds 
	total_asm = sum(len(record.seq) for record in records) # Total assembly size (including gaps)
	lengths = [len(record.seq) for record in records] # Scaffold lengths
	lengths.sort(reverse=True)
	gc_count = sum(record.seq.upper().count("G") + record.seq.upper().count("C") for record in records) # GC count
	all_count = sum(record.seq.upper().count("G") + record.seq.upper().count("C") + record.seq.upper().count("T") + record.seq.upper().count("A") for record in records) # Total bases (excluding gaps)
	gc_cont = (gc_count/all_count)*100 # GC content
	s_n50 = next((length for length in lengths if sum(lengths[:lengths.index(length) + 1]) >= total_asm * 0.5), None) # Scaffold N50
	s_n90 = next((length for length in lengths if sum(lengths[:lengths.index(length) + 1]) >= total_asm * 0.9), None) # Scaffold N90
	gap_sum = sum(record.seq.count("N") for record in records) # Total gap count
	gap_perc = round((gap_sum/total_asm)*100, 2) # Total gap count
	return count, s_n50, gap_sum, total_asm, s_n90, gc_cont, gap_perc

# Calculate common contig statistics, splitting the scaffolds on gaps of N's of user defined length
def calculate_contig_stats(fasta_file, gap_size):
	records = [record for record in SeqIO.parse(fasta_file, "fasta")]
	new_records = []
	for record in records:
		seq = str(record.seq)
		contigs = seq.split("N" * int(gap_size))
		for contig in contigs:
			if len(contig) > 0 and "N" not in contig:
				new_record = record[:]
				new_record.seq = Seq(contig)
				new_records.append(new_record)
	c_count = len(new_records) # Number of sequences
	gap_count = len(new_records) - len(records)
	total_len = sum(len(record.seq) for record in new_records) # Total length of all sequences
	lengths = [len(record.seq) for record in new_records] # Length of each sequence
	lengths.sort(reverse=True)
	c_n50 = next((length for length in lengths if sum(lengths[:lengths.index(length) + 1]) >= total_len * 0.5), None) # Contig N50
	c_n90 = next((length for length in lengths if sum(lengths[:lengths.index(length) + 1]) >= total_len * 0.9), None) # Contig N90
	return c_count, c_n50, gap_count, c_n90

def write_stats_to_file(isolate, file_name, count, s_n50, c_count, c_n50, gap_count, gap_sum, total_asm, s_n90, c_n90, gc_cont, gap_perc):
	with open(file_name, "w") as f:
		f.write(isolate+",assembly_length_bp,{},\n".format(total_asm))
		f.write(isolate+",scaffold_count,{},\n".format(count))
		f.write(isolate+",scaffold_N50_bp,{},\n".format(s_n50))
		f.write(isolate+",scaffold_N90_bp,{},\n".format(s_n90))
		f.write(isolate+",contig_count,{},\n".format(c_count))
		f.write(isolate+",contig_N50_bp,{},\n".format(c_n50))
		f.write(isolate+",contig_N90_bp,{},\n".format(c_n90))
		f.write(isolate+",GC_perc,{:.2f},\n".format(gc_cont))
		f.write(isolate+",gaps_count,{},\n".format(gap_count))
		f.write(isolate+",gaps_sum_bp,{},\n".format(gap_sum))
		f.write(isolate+",gaps_perc,{},\n".format(gap_perc))

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
