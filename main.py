import csv
import re
from src import ProteinToGene
from src import ProteinOnDNA
from src import DownstreamAnalyzer
from Bio import SeqIO, SeqRecord
from Bio.Data.CodonTable import TranslationError
from tqdm import tqdm

__author__ = "ronmor"
NEWLINE = '\n'
N = 1000


def find_accessions_number(subject_protein, output_handle):
    """
    :param subject_protein:
    :type subject_protein: str
    :param output_handle: file handle open for writing
    :type output_handle: _io.TextIOWrapper
    :return:
    """
    # Find accession number and db type
    proteins_accessions = {}
    manual_extraction = None
    for protein in re.findall("(\w{2,3})\|([A-Z_0-9]*\.\d)", subject_protein):
        db_type, accession_number = protein
        # refseq and SwissProt have different parsing
        if db_type in ('ref', 'sp'):
            if db_type == 'ref':
                manual_extraction = accession_number
                output_handle.write(manual_extraction + NEWLINE)
        else:
            proteins_accessions[accession_number] = db_type
    if not proteins_accessions and manual_extraction:
        pass
    return proteins_accessions


def main(input_csv_path):
    proteins_accessions_by_row = []
    dict_reader_keys = ("Source protein", "Subject protein", "Identities", "Positives", "Query match length",
                        "Unknown1", "Unknown2", "Query start", "Query end", "Subject start", "Subject length",
                        "E-value", "Bit score")
    out_handle = open('output_fasta', 'w')
    old_features_seqs = set()
    proteins_on_dna = set()
    with open(input_csv_path) as csvf:
        reader = csv.DictReader(csvf, fieldnames=dict_reader_keys)
        non_supported_accessions_file = open(r'manual_extraction.txt', 'w')
        for row in reader:
            accessions = row['Subject protein']
            proteins_accessions_by_row.append(find_accessions_number(accessions, non_supported_accessions_file))
    for row in proteins_accessions_by_row:
        for protein_accession in row:
            other_keys = list(row.keys())
            other_keys.remove(protein_accession)
            ProteinToGene(protein_accession, row[protein_accession], other_keys)
    print("fetching all!")
    ProteinToGene.fetchall()
    print("finished fetching")
    for prot_accession in tqdm(ProteinToGene.get_accession_to_object()):
        prot = ProteinToGene.get_accession_to_object()[prot_accession]
        assert isinstance(prot, ProteinToGene)
        prot.set_genome_accession_no()
        location = prot.extract_location_on_gene()
        if prot.genome_accession_no is None:
            continue
        if location is not None:
            protein_on_dna = ProteinOnDNA(prot.genome_accession_no, location, prot.accession_no)
        else:
            protein_on_dna = ProteinOnDNA(prot.genome_accession_no, protein_accession_no=prot.accession_no)
        proteins_on_dna.add(protein_on_dna)
    for mediator in tqdm(proteins_on_dna):
        assert isinstance(mediator, ProteinOnDNA)
        if mediator.protein_location is None:
            prot_location = mediator.find_protein_location()
            mediator.protein_location = prot_location
        next_nuc = mediator.get_next_n_nucleotide_of_feature_location(N)
        is_complementing = True if mediator.protein_location.strand == -1 else False
        downstreamer = DownstreamAnalyzer(next_nuc, 0, mediator.genbank_file, is_complementing)
        features_generator = downstreamer.generate_downstream_cdss()
        for feature in features_generator:
            try:
                translated_prot = downstreamer.translate_feature(feature)
            except TranslationError:
                continue
            if str(translated_prot) not in old_features_seqs:
                old_features_seqs.add(str(translated_prot))
                protein_id = feature.qualifiers.get('protein_id')
                protein_id = protein_id[0] if protein_id else protein_id
                desc = feature.qualifiers.get('product')
                desc = desc[0] if desc else desc
                SeqIO.write(SeqRecord.SeqRecord(translated_prot, id=str(protein_id), description=str(desc)), out_handle, 'fasta')


if __name__ == '__main__':
    main(r'AimR phi3t blast results NCBI.csv')
