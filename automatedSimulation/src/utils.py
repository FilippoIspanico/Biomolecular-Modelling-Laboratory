from Bio import SeqIO


def pdb2fasta(pdb_file: str) -> str:
    with open(pdb_file, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
           return str(record.seq)

