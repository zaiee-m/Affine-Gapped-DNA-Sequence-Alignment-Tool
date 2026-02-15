Affine Gapped Sequence Alignment Tool
A Python implementation of the Smith-Waterman-Gotoh algorithm for local nucleotide sequence alignment. This tool uses Affine Gap Penalties, which distinguish between opening a gap and extending one, providing more biologically relevant results for sequences with long insertions or deletions (indels).

Features
Local Alignment: Finds the most similar region between two sequences.

Affine Gap Scoring: Supports independent penalties for Gap Open and Gap Extend.

FASTA Support: Directly parses standard .fasta or .fa files.

Visual Output: Generates a clear alignment map highlighting matches, mismatches, and indels.

Installation
Ensure you have Python 3.x installed.

Clone this repository or save the script as align_tool.py.

No external libraries are required (uses argparse and os from the Python Standard Library).

Usage
The tool is executed via the command line. You must provide two FASTA files as positional arguments.

python align_tool.py <file1.fasta> <file2.fasta> [options]

Command Line Arguments
Argument | Description | Default
file1 | Path to the first FASTA file (Reference) | Required
file2 | Path to the second FASTA file (Query) | Required
--match | Score for a base match | 5.0
--mismatch | Penalty for a mismatch (negative) | -4.0
--open | Penalty for opening a new gap (negative) | -12.0
--extend | Penalty for extending an existing gap (negative) | -4.0

Example Command
python align_tool.py seq1.fasta seq2.fasta --match 5 --mismatch -4 --open -12 --extend -4

Understanding the Output
The tool provides a visual representation of the alignment using a three-line format:

| : Match (Identical nucleotides)

. : Mismatch (Substitution)

: Indel (Insertion or Deletion)

Example Visualization (Single Indel):
Seq1 (Ref): GCTAGTCA--GATCGACCGT
||||||||  ||||||||||
Seq2 (Alt): GCTAGTCACAGATCGACCGT

How it Works (Affine Logic)
Unlike linear gap penalties where every gap character costs the same, this tool uses the Gotoh Algorithm. It maintains three separate scoring matrices to track the state of the alignment:

M Matrix: Score for matching/mismatching base i and base j.

X Matrix: Score for a gap in sequence A (Insertion).

Y Matrix: Score for a gap in sequence B (Deletion).

The cost of a gap of length k is calculated as:
Score = Gap_Open + (k * Gap_Extend)

Testing
To validate the results, you can compare the output of this tool with the EBI LALIGN web server using the same scoring parameters.