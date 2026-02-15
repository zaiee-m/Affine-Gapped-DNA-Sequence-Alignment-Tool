import sys
import argparse
import os

def read_fasta(file_path):
    """
    Reads a FASTA file and returns the first sequence found.
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        sequence = []
        header = ""
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                if not header:
                    header = line[1:]
                # If we already found a sequence and hit another header, stop (taking 1st seq only)
                elif sequence: 
                    break
            else:
                sequence.append(line)
        
        return "".join(sequence).upper(), header
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)

def affine_local_alignment(seq1, seq2, match_score, mismatch_pen, gap_open, gap_extend):
    """
    Performs Smith-Waterman Local Alignment with Affine Gap Penalties.
    
    Matrices:
    M: Best score ending in a match/mismatch
    X: Best score ending in a gap in Seq1 (insertion in Seq2)
    Y: Best score ending in a gap in Seq2 (insertion in Seq1)
    """
    n = len(seq1)
    m = len(seq2)
    
    # Initialize matrices with -infinity (or very low number) for calculation
    # We use float('-inf') to ensure valid paths are chosen
    # But for Local Alignment, 0 is the floor.
    
    # DP Tables
    M = [[0.0] * (m + 1) for _ in range(n + 1)]
    X = [[0.0] * (m + 1) for _ in range(n + 1)]
    Y = [[0.0] * (m + 1) for _ in range(n + 1)]
    
    # Traceback Tables: 1=Diag, 2=Up, 3=Left, 0=Stop
    # We need separate traceback pointers for each state to handle affine jumps correctly
    # However, for simplicity in a single script, we often just recalculate or store limited info.
    # Here we will store specific pointers.
    # M_tb stores where M came from: 0:Stop, 1:M(i-1,j-1), 2:X(i-1,j-1), 3:Y(i-1,j-1)
    M_tb = [[0] * (m + 1) for _ in range(n + 1)]
    
    max_score = 0
    max_pos = (0, 0)
    
    # Fill Matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s1 = seq1[i-1]
            s2 = seq2[j-1]
            score_func = match_score if s1 == s2 else mismatch_pen
            
            # --- Update X (Gap in seq1 / Insertion in seq2) ---
            # Extend existing gap (X[i-1][j]) or Open new gap from M (M[i-1][j])
            # Cost model: Gap_Cost(k) = Open + k * Extend
            
            # Case 1: Extend gap in X
            x_ext_score = X[i-1][j] + gap_extend
            
            # Case 2: Open gap from M
            x_open_score = M[i-1][j] + gap_open + gap_extend
            
            X[i][j] = max(x_ext_score, x_open_score)
            # For local alignment, gaps can't inherently drop below 0 if M resets to 0, 
            # but strictly X/Y track gap states. We floor them at 0 only if M floors at 0.
            # Actually, X and Y components shouldn't be floored individually in strict affine 
            # unless we allow starting a local alignment with a gap (uncommon).
            # We will keep them distinct.
            
            # --- Update Y (Gap in seq2 / Insertion in seq1) ---
            y_ext_score = Y[i][j-1] + gap_extend
            y_open_score = M[i][j-1] + gap_open + gap_extend
            Y[i][j] = max(y_ext_score, y_open_score)
            
            # --- Update M (Match/Mismatch) ---
            # Match comes from diagonal of M, X, or Y
            m_from_m = M[i-1][j-1] + score_func
            m_from_x = X[i-1][j-1] + score_func
            m_from_y = Y[i-1][j-1] + score_func
            
            best_m = max(m_from_m, m_from_x, m_from_y, 0) # 0 is for Local Alignment (Start new)
            M[i][j] = best_m
            
            # Store Traceback Info for M
            if best_m == 0:
                M_tb[i][j] = 0
            elif best_m == m_from_m:
                M_tb[i][j] = 1 # Match
            elif best_m == m_from_x:
                M_tb[i][j] = 2 # From X
            else:
                M_tb[i][j] = 3 # From Y
                
            # Track Global Max
            if best_m > max_score:
                max_score = best_m
                max_pos = (i, j)

    # --- Traceback ---
    align1 = []
    align2 = []
    i, j = max_pos
    
    # We assume we end in M state (best local alignments usually end on a match/mismatch)
    # If the max score was in X or Y, we'd need to handle that, but M[i][j] checks X and Y diagonals.
    
    curr_state = 1 # Start assuming we trace back from M
    
    while i > 0 and j > 0 and M[i][j] > 0:
        if curr_state == 1: # In M state
            type_src = M_tb[i][j]
            if type_src == 0: break # Stop
            
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            
            if type_src == 1: curr_state = 1
            elif type_src == 2: curr_state = 2
            elif type_src == 3: curr_state = 3
            
            i -= 1
            j -= 1
            
        elif curr_state == 2: # In X state (gap in seq1, i decreases, j stays)
            # Did we extend X or open X?
            # Recalculate to decide path (memory saving vs complexity trade-off)
            if round(X[i][j] - (X[i-1][j] + gap_extend), 5) == 0:
                # Extended X
                curr_state = 2
            else:
                # Opened X from M
                curr_state = 1
            
            align1.append(seq1[i-1])
            align2.append("-")
            i -= 1
            
        elif curr_state == 3: # In Y state (gap in seq2, j decreases, i stays)
            if round(Y[i][j] - (Y[i][j-1] + gap_extend), 5) == 0:
                curr_state = 3
            else:
                curr_state = 1
                
            align1.append("-")
            align2.append(seq2[j-1])
            j -= 1

    return "".join(reversed(align1)), "".join(reversed(align2)), max_score

def print_alignment(seq1_aligned, seq2_aligned, score):
    """
    Pretty prints the alignment with match bars.
    """
    
    match_line = []
    matches = 0
    gaps = 0
    length = len(seq1_aligned)
    
    for k in range(length):
        s1 = seq1_aligned[k]
        s2 = seq2_aligned[k]
        
        if s1 == s2 and s1 != '-':
            match_line.append("|") # Match
            matches += 1
        elif s1 == '-' or s2 == '-':
            match_line.append(" ") # Indel (Gap)
            gaps += 1
        else:
            match_line.append(".") # Mismatch
            
    match_str = "".join(match_line)
    
    print("\n" + "="*60)
    print("FINAL ALIGNMENT VISUALIZATION")
    print("="*60)
    print(f"Score: {score:<10} Length: {length}")
    print(f"Identity: {matches}/{length} ({matches/length*100:.1f}%)")
    print(f"Gaps:     {gaps}/{length} ({gaps/length*100:.1f}%)")
    print("-" * 60)
    
    # Print in blocks for readability
    block_size = 60
    for k in range(0, length, block_size):
        end = min(k + block_size, length)
        print(f"Seq1 (Ref): {seq1_aligned[k:end]}")
        print(f"            {match_str[k:end]}")
        print(f"Seq2 (Alt): {seq2_aligned[k:end]}")
        print()

def main():
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        description="Nucleotide Affine Gap Local Alignment Tool (Smith-Waterman-Gotoh)"
    )
    
    # Positional arguments (Required)
    parser.add_argument("file1", help="Path to the first FASTA file")
    parser.add_argument("file2", help="Path to the second FASTA file")
    
    # Optional arguments (Scoring parameters with defaults)
    parser.add_argument("--match", type=float, default=5.0, 
                        help="Score for a match (default: 5.0)")
    parser.add_argument("--mismatch", type=float, default=-4.0, 
                        help="Penalty for a mismatch (default: -4.0)")
    parser.add_argument("--open", type=float, default=-12.0, 
                        help="Penalty for opening a gap (default: -12.0)")
    parser.add_argument("--extend", type=float, default=-4.0, 
                        help="Penalty for extending a gap (default: -4.0)")

    # Parse arguments
    args = parser.parse_args()

    # Validate file existence before processing
    if not os.path.exists(args.file1):
        print(f"Error: File '{args.file1}' not found.")
        sys.exit(1)
    if not os.path.exists(args.file2):
        print(f"Error: File '{args.file2}' not found.")
        sys.exit(1)

    print("--- Nucleotide Affine Gap Local Alignment Tool ---")
    print(f"Parameters: Match={args.match}, Mismatch={args.mismatch}, "
          f"Open={args.open}, Extend={args.extend}")

    # Read Sequences
    # Note: We pass the filenames from the args namespace
    s1, h1 = read_fasta(args.file1)
    s2, h2 = read_fasta(args.file2)
    
    print(f"\nAligning:\n> {h1} ({len(s1)} bp)\n> {h2} ({len(s2)} bp)")
    
    # Perform Alignment
    # We pass the scoring parameters from the args namespace
    a1, a2, score = affine_local_alignment(
        s1, s2, 
        args.match, 
        args.mismatch, 
        args.open, 
        args.extend
    )
    
    # Output
    print_alignment(a1, a2, score)

if __name__ == "__main__":
    main()