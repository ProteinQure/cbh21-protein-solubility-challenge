import argparse
from Bio import SeqIO

# set up argument parsing (make sure these match those in config.yml)
parser = argparse.ArgumentParser()
parser.add_argument("--infile", type=str, required=True)
args = parser.parse_args()

# output file, this will be used for benchmarking
predictions_outfile = "predictions.csv"

# process input
fasta_dict = dict()
with open(args.infile, 'r') as fh:
    for record in SeqIO.parse(fh, 'fasta'):
        fasta_dict[record.id] = str(record.seq)

# save predictions to file
with open(predictions_outfile, 'w') as fh:
    fh.write("name,adenine_count\n")
    for key, value in fasta_dict.items():
        adenine_prediction = int(0.3 * len(value))
        fh.write(f"{key},{adenine_prediction}\n")

# print predictions to screen
print(open(predictions_outfile, 'r').read())
