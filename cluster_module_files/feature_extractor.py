# K-mer feature extraction
import numpy as np
import pandas as pd 
from Bio import SeqIO

def kmer_count(sequence, k=3):
    kmer_dict = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1
    return kmer_dict

def generate_kmer_matrix(sequences, k=3):
    kmer_matrices = []
    for sequence in sequences:
        kmer_matrices.append(kmer_count(sequence, k))
    return kmer_matrices

class PseudoAminoAcidComposition:
    def __init__(self, fasta_data, lambda_=30, w=0.05):
        self.lambda_ = lambda_
        self.w = w
        self.physicochemical_properties = {
            'hydrophobicity': {'A': 0.62, 'C': 0.29, 'D': -0.90, 'E': -0.74, 'F': 1.19, 'G': 0.48, 'H': -0.40, 'I': 1.38, 'K': -1.50, 'L': 1.06, 'M': 0.64, 'N': -0.78, 'P': 0.12, 'Q': -0.85, 'R': -2.53, 'S': -0.18, 'T': -0.05, 'V': 1.08, 'W': 0.81, 'Y': 0.26},
            'hydrophilicity': {'A': -0.5, 'C': -1.0, 'D': 3.0, 'E': 3.0, 'F': -2.5, 'G': 0.0, 'H': -0.5, 'I': -1.8, 'K': 3.0, 'L': -1.8, 'M': -1.3, 'N': 0.2, 'P': 0.0, 'Q': 0.2, 'R': 3.0, 'S': 0.3, 'T': -0.4, 'V': -1.5, 'W': -3.4, 'Y': -2.3},
            'mass': {'A': 15.0, 'C': 47.0, 'D': 59.0, 'E': 73.0, 'F': 91.0, 'G': 1.0, 'H': 82.0, 'I': 57.0, 'K': 73.0, 'L': 57.0, 'M': 75.0, 'N': 58.0, 'P': 42.0, 'Q': 72.0, 'R': 101.0, 'S': 31.0, 'T': 45.0, 'V': 43.0, 'W': 130.0, 'Y': 107.0}
        }
        self.amino_acids = {idx + 1: aa for idx, aa in enumerate("ACDEFGHIKLMNPQRSTVWY")}
        self.collect_features(fasta_data)

    def collect_features(self, fasta_data):
        features = []
        features.append(['#'] + list(self.amino_acids.values()) + [f'λ{i}' for i in range(1, self.lambda_ + 1)])
        for i in range(0, len(fasta_data), 2):
            accession = fasta_data[i].strip()
            sequence = fasta_data[i + 1].strip()
            if len(sequence) < self.lambda_:
                print(f"Sequence {accession} skipped due to insufficient length.")
                continue
            paac = self.calculate_paac(sequence)
            features.append([accession] + paac)
        self.df = pd.DataFrame(features[1:], columns=features[0])

    def calculate_paac(self, sequence):
        lower_theta = self.calculate_lower_theta(sequence)
        denominator = 1 + (self.w * sum(lower_theta.values()))
        paac = []
        for i in range(1, 21 + self.lambda_):
            if i <= 20:
                numerator = sequence.count(self.amino_acids[i]) / len(sequence)
                paac.append(numerator / denominator)
            else:
                numerator = self.w * lower_theta[i - 20]
                paac.append(numerator / denominator)
        return paac

    def calculate_lower_theta(self, sequence):
        lower_theta = {}
        for i in range(1, self.lambda_ + 1):
            if len(sequence) <= i:
                lower_theta[i] = 0
            else:
                lower_theta[i] = (1 / (len(sequence) - i)) * self.calculate_upper_theta(sequence, i)
        return lower_theta

    def calculate_upper_theta(self, sequence, i):
        upper_theta = []
        for j in range(len(sequence) - i):
            diff = [
                (self.physicochemical_properties[prop][sequence[j]] -
                 self.physicochemical_properties[prop][sequence[j + i]]) ** 2
                for prop in self.physicochemical_properties
            ]
            upper_theta.append(sum(diff) / len(self.physicochemical_properties))
        return sum(upper_theta)
    
def load_sequences(file_path):
    return [str(record.seq) for record in SeqIO.parse(file_path, "fasta")]

def combine_features(sequences, k=3, lambda_=30, w=0.05):
    kmer_features = generate_kmer_matrix(sequences, k)

    fasta_data = []
    for idx, seq in enumerate(sequences):
        if len(seq) >= lambda_:
            fasta_data.extend([f">seq_{idx}", seq])
        else:
            print(f"Skipping sequence seq_{idx} (length < {lambda_})")

    paac = PseudoAminoAcidComposition(fasta_data, lambda_, w)
    paac_features = paac.df.iloc[:, 1:].values.tolist()

    combined_features = []
    max_length = max(len(list(kmer.values()) + paac) for kmer, paac in zip(kmer_features, paac_features))

    for kmer, paac in zip(kmer_features, paac_features):
        feature = list(kmer.values()) + paac
        feature += [0] * (max_length - len(feature))  # Pad with zeros
        combined_features.append(feature)

    return np.array(combined_features, dtype=float)
