import re
import pandas as pd
from Bio import Blast
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import yaml
from pathlib import Path


class PositionScanner:

    def __init__(self, yaml_configuration_path: str = "config.yaml"):

        try:
            config_file = open(yaml_configuration_path, "r")
            self.yaml_configuration = yaml.safe_load(config_file)
            self.pdb: str = self.yaml_configuration["pdb"]

            print("Inizialing Position Scanner: ")
            print(f"Current pdb: {self.pdb}")

        except:
            print(f'Cannot open input file: {yaml_configuration_path}')

    class Mutation:

        """
        This class represents a mutation: an aminoacid in a specific position.
        """

        def __init__(self, position: int, amino_acid: str):
            if len(amino_acid) != 1:
                raise ValueError('Not valid amino_acid input: please use one letter code for specifing the aminoacid')

            self.position = position

            self.amino_acid = amino_acid

        def __init__(self, foldX_name: str):
            self.position, self.amino_acid = self.toSequence(foldX_name)

        def toFoldX(self):
            """
            Return the mutation in FoldX format i.e. amino_acid + "A" + position
            """
            return self.amino_acid + "A" + str(self.position)

        def __repr__(self):
            return f"({self.position}, {self.amino_acid})"

        def toSequence(self, foldX_name):
            amino_acid = foldX_name[-1]
            position = re.findall(r'\d+', foldX_name)

            return position[0], amino_acid

    def foldX_scan(self, a: int, b: int):
        """
        This function uses foldex PositionScan to scan the position in [a, b] and returns the best ones.
        FoldX only work if the protein to be mutated resides in the same folder of such program. 
        For this reason in this function we will first copy the current pdb (specified in the yaml file) into such folder, 
        and then run PositionScan. We will also need to modify the configuration file. Due to the large number of files produced by FoldX, 
        each time the program is called (which is b-a) the produced files will be removed. Once the best position has been identified 
        """

    def parse_FXPositionScan_output(self, epsilon: float):
        """
        This method parses PositionsCan output file. It returns mutations that:
        1. have a negative energy gain. In particular energy_gain<=epsilon 
        2. have a negative energy gain but epislon<energy_gain<0 and the proposed aminoacid in that position is the most common one (for that position). 
        """

        scanning_output_name = "/PS_" + self.pdb + "_scanning_output.txt"
        scanning_output_name = self.yaml_configuration["foldX_installation_folder"] + scanning_output_name

        # scanning_output = open(scanning_output_name)

        scanning_output = pd.read_csv(scanning_output_name, header=None, sep="\t")
        scanning_output.columns = ["Mutation", "EnergyGain"]

        mask = scanning_output["EnergyGain"] < epsilon

        for idx, foldXmutation in enumerate(scanning_output["Mutation"]):
            print(idx, foldXmutation)
            if foldXmutation[-1] == self.yaml_configuration["most_common_sequence"][idx]:
                mask[idx] = True

        result = []
        for mutation in scanning_output[mask]["Mutation"]:
            result.append(self.Mutation(foldX_name=mutation))

        print(result)


class BLAST:

    def __init__(self, yaml_configuration_path: str):

        config_file = open(yaml_configuration_path, "r")
        self.yaml_configuration = yaml.safe_load(config_file)

        self.fasta_string = self.yaml_configuration["blast"]["fasta_query"]
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
            fout.write(f"{idx+1},{score}\n")


        fout.close()

