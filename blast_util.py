#!/usr/bin/python2.7
"""

Performs a web BLAST query against NR using DNA sequences read from a fasta
file, storing the results in a SQLite database.

Dependencies: BioPython
    BioPython may be installed using the Python package index:
    $ sudo easy_install -f http://biopython.org/DIST/ biopython

usage: blast_util [-h] --input-file INPUT_FILE [--database DATABASE]
                 [--output-folder OUTPUT_FOLDER] [--limit LIMIT]
                 [--matrix MATRIX] [--e-value E_VALUE]

optional arguments:
  -h, --help            show this help message and exit
  --input-file INPUT_FILE
                        path to an input file containing fasta sequences
  --database DATABASE   name of the BLAST database you wish to query against
                        (default is NR)
  --output-folder OUTPUT_FOLDER
                        path to folder for output SQLite databases (default is
                        working directory)
  --limit LIMIT         limit your BLAST query to N results
  --matrix MATRIX       choose an alternative substitution matrix, selecting
                        from [PAM30, PAM70, BLOSUM80, BLOSUM45]
  --e-value E_VALUE     choose an e-value threshold for the BLAST search.
                        Default is 10.0

@author: Erik Edlund <eedlun@gmail.com>
"""


import argparse
import os
import sqlite3
import sys

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from operator import attrgetter


class BLASTInterface:
    """
    Utility class for interacting with NCBI's BLAST over the web.
    """
    def __init__(self, blast_db, limit=100, matrix=None, e_value_threshold=10.0):
        """
        blast_db: database to query (e.g. NR)
        limit: limit N results returned from NCBI
        matrix: scoring matrix (default: let NCBI choose)
        e_value_threshold: an expect value cutoff
        """
        self.blast_db = blast_db
        self.limit = limit
        self.matrix = matrix
        self.e_value_threshold = e_value_threshold

    def blast(self, fasta_seq):
        """
        Perform a web BLAST with a fasta sequence, returning a BLASTResult object
        """
        response_xml = NCBIWWW.qblast("blastn",
                                              database=self.blast_db,
                                              sequence=fasta_seq.seq,
                                              hitlist_size=self.limit,
                                              matrix_name=self.matrix,
                                              expect=self.e_value_threshold)
        blast_response = NCBIXML.read(response_xml)

        return self.parse(blast_response, fasta_seq)

    def parse(self, result, seq):
        return BLASTResult(result, seq)

class BLASTResult:
    """
    An instance stores relevant data from a BLAST search and provides
    a means of writing them to a SQLite database.
    """
    def __init__(self, record, seq):

        self.descriptions = record.descriptions
        self.alignments = record.alignments
        self.database = record.database
        self.database_length = record.database_length
        self.fasta_seq = seq

    def write_to_db(self, db):
        conn = sqlite3.connect(db)
        conn.execute('''CREATE TABLE IF NOT EXISTS results (sequence_id text, description text, percent_id real, e_value real)''')
        #append results
        for result_align in self.alignments:
            #record data for the highest scoring HSP
            hsp_best = max(result_align.hsps, key=attrgetter('identities'))
            percent_id = (100 * hsp_best.identities / len(self.fasta_seq.seq))
            conn.execute("INSERT INTO results VALUES ('%s', '%s', '%s', '%s')" %(self.fasta_seq.id, result_align.title, percent_id, float(hsp_best.expect)))
        #commit the changes
        conn.commit()
        conn.close()



def create_arg_parser():
    """
    Create an argument parser for the script.

    """
    d = """Performs a web BLAST query against NR using DNA sequences read from a fasta file, storing the results in a
    SQLite database.

    Dependencies: BioPython
        BioPython may be installed using the Python package index:
        $ sudo easy_install -f http://biopython.org/DIST/ biopython
    """
    p = argparse.ArgumentParser(prog="blast_util", description=d)
    p.add_argument("--input-file", help="path to an input file containing fasta sequences", required=True )
    p.add_argument("--database", help="name of the BLAST database you wish to query against (default is NR)", default="nr" )
    p.add_argument("--output-folder", help="path to folder for output SQLite databases (default is working directory)", default=None)
    p.add_argument("--limit", help="limit your BLAST query to N results", type=int, default=100)
    p.add_argument("--matrix", help="choose an alternative substitution matrix, selecting from [PAM30, PAM70, BLOSUM80, BLOSUM45]", default=None)
    p.add_argument("--e-value", help="choose an e-value threshold for the BLAST search. Default is 10.0", default=10.0)
    return p

def main(argv):
    try:
        args = create_arg_parser().parse_args(argv)

        input_file = args.input_file

        blaster = BLASTInterface(args.database, args.limit, args.matrix, args.e_value)

        output_dir = args.output_folder if args.output_folder else os.getcwd()
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        db_file = os.path.join(output_dir, '%s_blast_results.db' % os.path.basename(input_file))
        if not os.path.exists(input_file): raise Exception('Could not find file: %s' % input_file )

        #Run BLAST on each fasta record one at a time, recording the results in the DB
        for fasta_record in SeqIO.parse(input_file, 'fasta'):
            print 'Running BLAST for sequence %s...' % fasta_record.id
            record_results = blaster.blast(fasta_record)
            print 'Writing results to %s...' % db_file
            record_results.write_to_db(db_file)

        print 'Done.'
        return 0
    except KeyboardInterrupt:
        print '\nBLAST search aborted.'
        return 0
    except Exception as e:
        sys.stderr.write('Operation failed:\n')
        sys.stderr.write(str(e) + '\n')
        return 1

if __name__ == '__main__':
    #Pass along all relevant program args.
    p = os.path.basename( os.path.abspath(__file__) )
    x = [x for x in sys.argv if x.endswith(p)].pop(0)
    argv = sys.argv[sys.argv.index(x) + 1:]
    sys.exit( main(argv) )
