import pandas as pd
from Bio import Blast
import xml.etree.ElementTree as ET
import yaml
from pathlib import Path
import subprocess
import os
from src.utils import pdb2fasta
import glob


class PositionScanner:

    def __init__(self, yaml_configuration_path: str, pdb_name: str, eps: float = -0.4):

        config_file = open(yaml_configuration_path, "r")
        self.yaml_configuration = yaml.safe_load(config_file)

        self.foldX_program: Path = Path(self.yaml_configuration["foldX"]["program_path"])

        if not self.foldX_program.exists():
            print(f"Can't find foldX at the indicated directory: {self.foldX_program}")
        else:
            print(f"Found foldX at the indicated directory: {self.foldX_program}")

        self.pdb_name: str = pdb_name

        self.temperature: int = self.yaml_configuration["foldX"]["temperature"]

        self.log_path: Path = Path('log/foldX/')
        self.log_path.mkdir(parents=True, exist_ok=True)

        self.log_file_path = self.log_path / 'foldX.log'

        self.output_dir: Path = self.log_path / 'output/'
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.rotabase_dir = self.output_dir / "rotabase.txt"

        self.pdb_dir: Path = Path('proteins')
        self.pdb_dir.mkdir(parents=True, exist_ok=True)

        self.eps = eps

    def foldX_scan(self, positions: list[int]):
        """
        This function uses foldex PositionScan to scan the input positions and returns the best ones.
        FoldX only work if the protein to be mutated resides in the same folder of such program. 
        For this reason in this function we will first copy the current pdb (specified in the yaml file) into such folder, 
        and then run PositionScan. We will also need to modify the configuration file. Due to the large number of files produced by FoldX, 
        each time the program is called (which is b-a) the produced files will be removed. Once the best position has been identified 
        """

        fasta_sequence: str = pdb2fasta(f"proteins/{self.pdb_name}.pdb")
        mutations: str = ""

        for position in positions:
            mutations = mutations + f"{fasta_sequence[position]}A{position+1}a,"

        if len(positions) > 0:

            mutations = mutations[:-1]  # we remove the hanging comma

            f = open(self.log_file_path, 'w')

            subprocess.run(
                [self.foldX_program,
                 "--command=PositionScan",
                 f"--pdb-dir={self.pdb_dir}",
                 f"--output-dir={self.output_dir}",
                 f"--rotabaseLocation={self.rotabase_dir}",
                 f"--pdb={self.pdb_name}.pdb",
                 f"--positions={mutations}",
                 f"--temperature={self.temperature}",
                 "--screen=false"  # this option suppress the verbose mode

                 ],
                stdout=f,
                stderr=subprocess.STDOUT
            )

    def parse_scan(self):

        aa3to1 = {
            'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M',
            'ILE': 'I', 'LEU': 'L', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
            'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H',
            'CYS': 'C', 'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G'
        }

        aa1to3 = {v: k for k, v in aa3to1.items()}

        output_file: Path = self.output_dir / f"PS_{self.pdb_name}_scanning_output.txt"
        df = pd.read_table(output_file, header=None)
        df.columns = ["mutation", "energy"]
        potential_mutations = df[df["energy"] < self.eps]

        potential_mutations.to_csv('results/foldX.csv', mode='a', header=False)

        # first 4 char indicates the initial aminoacid
        # the last char indicates the mutation amino acid
        # the remaining part indicates the position


        for mutation in potential_mutations["mutation"]:
            amino_acid = mutation[-1]
            position = mutation[4:-1]
            mutation_filename = f"{aa1to3[amino_acid]}{position}_{self.pdb_name}.pdb"
            if os.path.exists(mutation_filename):
                os.rename(mutation_filename, f"proteins/{mutation_filename}")

        pattern = "*.pdb"
        for file_path in glob.glob(pattern):
            os.remove(file_path)

        pattern = "*.txt"
        for file_path in glob.glob(pattern):
            os.remove(file_path)


class BLAST:

    def __init__(self, yaml_configuration_path: str, pdb_name: str):

        config_file = open(yaml_configuration_path, "r")
        self.yaml_configuration = yaml.safe_load(config_file)



        self.fasta_string: str = pdb2fasta(f"proteins/{pdb_name}.pdb")

        self.alignments = self.yaml_configuration["blast"]["alignments"]

        log_folder: Path = Path("log/BLAST/")
        self.blast_output_filename: Path = log_folder / f"blast_{self.alignments}.xml"

        self.sequences: list[str] = []
        Blast.email = self.yaml_configuration["blast"]["email"]

    def execute(self):
        print(f"""Executing blast query. 
        Number of expected results: {self.alignments}
        """)

        print("Waiting for blast query ...")

        result_stream = Blast.qblast("blastp", "nr", self.fasta_string, hitlist_size=self.alignments)

        with open(self.blast_output_filename, "wb") as out_stream:
            out_stream.write(result_stream.read())

        result_stream.close()

        print(f"Written blast query result at: {self.blast_output_filename}")

    def get_sequences(self):

        field = "Hsp_hseq"

        try:
            # Parse the XML file
            tree = ET.parse(self.blast_output_filename)
            root = tree.getroot()

            # Find all Hsp elements
            for hsp in root.findall(".//Hsp"):
                # Extract the specified field value
                element = hsp.find(field)
                if element is not None:
                    self.sequences.append(element.text)
        except ET.ParseError as e:
            print(f"Error parsing XML: {e}")
        except FileNotFoundError as e:
            print(f"File not found: {e}")

        return self.sequences

    @staticmethod
    def scores(protein_sequence: str, similar_sequences: list[str]):
        """
        This function returns the score of every position of the protein, based on the alignments specified in the input file.

        :param protein_sequence: the sequence of the protein, for which position you want to compute the scores.
        :param similar_sequences: the sequences found by a BLAST query
        :return: an array containing the scores.
        """

        input_sequence_length = len(protein_sequence)
        if input_sequence_length == 0:
            raise ValueError('No input sequence found. Cannot compute score')

        if len(similar_sequences) == 0:
            raise ValueError('No similar sequences found. Cannot compute score')

        score = []

        for position in range(input_sequence_length):
            reference_amino_acid = similar_sequences[0][position]
            reference_amino_acid_occurrences = 0

            for sequence in similar_sequences:
                if position < len(sequence) and sequence[position] == reference_amino_acid:
                    reference_amino_acid_occurrences += 1

            score.append(1 - reference_amino_acid_occurrences / len(similar_sequences))

        return score

    def get_scores(self):
        self.get_sequences()
        return self.scores(protein_sequence=self.fasta_string, similar_sequences=self.sequences)

    def write_scores(self):
        write_folder: Path = Path(self.yaml_configuration["blast"]["write_folder"])
        scores_file: Path = write_folder / "positions_scores.csv"

        scores = self.get_scores()

        fout = open(scores_file, 'w')
        fout.write("idx,score\n")

        for idx, score in enumerate(scores):
            fout.write(f"{idx},{score}\n")

        fout.close()
