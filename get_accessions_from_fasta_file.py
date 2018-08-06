import re
ACCESSION_REGEX = re.compile('^>([A-z0-9_]+\.\d)', re.M)


def make_accessions_list(fasta_file_path):
    with open(fasta_file_path) as f:
        fasta_text = f.read()
    return ACCESSION_REGEX.findall(fasta_text)

if __name__ == '__main__':
    fasta_path = r"Z:\ronmor\mature_file"
    accessions_list = make_accessions_list(fasta_path)
    for accession in accessions_list:
        print(accession)
