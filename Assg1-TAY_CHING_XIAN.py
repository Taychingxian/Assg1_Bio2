import streamlit as st
import numpy as np
import pandas as pd

# --- 1. Core Alignment Algorithms (Needleman-Wunsch & Smith-Waterman) ---

def initialize_matrix(seq1_len, seq2_len, gap_penalty, is_local=False):
    """Initializes the score matrix and traceback matrix."""
    rows = seq1_len + 1
    cols = seq2_len + 1
    # Initialize score matrix with zeros
    matrix = np.zeros((rows, cols), dtype=int)
    # Initialize traceback matrix (stores pointers for alignment path)
    # 0: Stop/End (Local), 1: Diagonal (Match/Mismatch), 2: Up (Gap in seq2), 3: Left (Gap in seq1)
    traceback = np.zeros((rows, cols), dtype=int)

    if not is_local:
        # Global Alignment (Needleman-Wunsch) Initialization
        # Fill first row and column with gap penalties [cite: 134, 163]
        for i in range(1, rows):
            matrix[i, 0] = matrix[i-1, 0] + gap_penalty
            traceback[i, 0] = 2  # Arrow up
        for j in range(1, cols):
            matrix[0, j] = matrix[0, j-1] + gap_penalty
            traceback[0, j] = 3  # Arrow left
    
    # For Local Alignment (Smith-Waterman), first row and column remain 0 [cite: 346]
    return matrix, traceback

def perform_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, is_local=False):
    """
    Performs the sequence alignment (Global or Local) using Dynamic Programming.
    
    Global: Needleman-Wunsch [cite: 64, 132]
    Local: Smith-Waterman [cite: 51, 314]
    """
    m, n = len(seq1), len(seq2)
    matrix, traceback = initialize_matrix(m, n, gap_penalty, is_local)
    
    max_score = 0
    max_pos = (0, 0) # Used for Smith-Waterman traceback start [cite: 374]

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            
            # 1. Diagonal Move (Match/Mismatch) [cite: 203, 354]
            current_score = match_score if seq1[i-1] == seq2[j-1] else mismatch_score
            score_diag = matrix[i-1, j-1] + current_score
            
            # 2. Up Move (Gap in seq2) [cite: 204, 355]
            score_up = matrix[i-1, j] + gap_penalty
            
            # 3. Left Move (Gap in seq1) [cite: 205, 356]
            score_left = matrix[i, j-1] + gap_penalty
            
            # Find the maximum score among the three moves [cite: 206, 357]
            scores = [score_diag, score_up, score_left]
            
            if is_local:
                # Local Alignment: Score cannot be negative [cite: 346]
                scores.append(0) 
            
            max_score_cell = max(scores)
            matrix[i, j] = max_score_cell

            # Update Traceback Matrix (Store which direction gave the max score)
            # Preference order: Diagonal (1), Up (2), Left (3), Stop (0 - Local only)
            if is_local and max_score_cell == 0:
                traceback[i, j] = 0
            elif max_score_cell == score_diag:
                traceback[i, j] = 1 # Diagonal
            elif max_score_cell == score_up:
                traceback[i, j] = 2 # Up
            elif max_score_cell == score_left:
                traceback[i, j] = 3 # Left
            
            # For Local Alignment, track the overall maximum score and position [cite: 360, 361]
            if is_local and matrix[i, j] > max_score:
                max_score = matrix[i, j]
                max_pos = (i, j)

    final_score = max_score if is_local else matrix[m, n]
    return final_score, matrix, traceback, max_pos

def traceback_alignment(seq1, seq2, matrix, traceback, max_pos, is_local=False):
    """
    Traces back through the matrix to reconstruct the optimal alignment path.
    """
    aligned_seq1 = ""
    aligned_seq2 = ""
    
    # Global Alignment starts at bottom-right corner [cite: 210]
    # Local Alignment starts at the position of the max score [cite: 374]
    i, j = max_pos if is_local else (len(seq1), len(seq2))
    
    path = []
    
    while i > 0 or j > 0:
        if is_local and traceback[i, j] == 0:
            break # Stop at zero score for local alignment [cite: 374]

        current_path_cell = (i, j)
        path.append(current_path_cell)
        
        move = traceback[i, j]
        
        if move == 1: # Diagonal Move [cite: 285]
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif move == 2: # Up Move (Gap in seq2) [cite: 284]
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        elif move == 3: # Left Move (Gap in seq1) [cite: 284]
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        elif not is_local and (i == 0 or j == 0):
            # Handles the end of the global alignment when one index hits 0
            if i > 0:
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                i -= 1
            elif j > 0:
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                j -= 1
        else:
            # Should not happen in correctly initialized matrix
            break
    
    path.append((i, j)) # Add the final (0,0) or start of local alignment path
    
    # Reverse the path list for display from start to end
    return aligned_seq1, aligned_seq2, path[::-1]

# --- 2. Streamlit Application Layout ---

st.set_page_config(layout="wide", page_title="Sequence Alignment Tool")

st.title("🧬 Dynamic Programming Sequence Alignment")
st.markdown("Implemention of **Needleman-Wunsch** (Global) and **Smith-Waterman** (Local) algorithms.")
 # Trigger diagram for context [cite: 2]

# --- Input and Settings Sidebar ---
with st.sidebar:
    st.header("⚙️ Alignment Settings")
    
    # Sequence Inputs
    st.subheader("Sequences")
    seq_a = st.text_input("Sequence 1", value="ACGT", max_chars=50).upper().strip()
    seq_b = st.text_input("Sequence 2", value="CATG", max_chars=50).upper().strip()

    # Alignment Type
    st.subheader("Alignment Type")
    alignment_type = st.radio("Choose Algorithm:", ["Global Alignment (Needleman-Wunsch)", "Local Alignment (Smith-Waterman)"])
    
    is_local = alignment_type == "Local Alignment (Smith-Waterman)"

    # Scoring Metrics [cite: 34]
    st.subheader("Scoring Metrics")
    match_score = st.number_input("Match Score", value=2, min_value=0, step=1)
    mismatch_score = st.number_input("Mismatch Score (must be < 0 for effective scoring)", value=-1, max_value=0, step=1)
    gap_penalty = st.number_input("Gap Penalty (must be < 0)", value=-2, max_value=0, step=1)

    if not seq_a or not seq_b:
        st.error("Please enter both sequences.")
        st.stop()

# --- 3. Run Alignment and Display Results ---

st.header(f"Results: {alignment_type}")
st.write(f"**Sequence 1:** `{seq_a}`")
st.write(f"**Sequence 2:** `{seq_b}`")
st.write(f"**Scoring:** Match: `{match_score}`, Mismatch: `{mismatch_score}`, Gap: `{gap_penalty}`")

try:
    # 1. Perform Alignment
    score, matrix, traceback, max_pos = perform_alignment(
        seq_a, seq_b, match_score, mismatch_score, gap_penalty, is_local
    )

    # 2. Traceback
    aligned_seq1, aligned_seq2, path = traceback_alignment(
        seq_a, seq_b, matrix, traceback, max_pos, is_local
    )

    # --- Alignment Score and Sequences ---
    st.subheader("✅ Alignment Summary")
    
    # Display the score
    if is_local:
        st.metric("Optimal Local Alignment Score", value=score, help="The highest score found in the matrix (Smith-Waterman).")
    else:
        st.metric("Optimal Global Alignment Score", value=score, help="The score in the bottom-right cell (Needleman-Wunsch).")
    
    # Display the aligned sequences
    st.markdown("#### Aligned Sequences")
    st.code(f"Seq 1: {aligned_seq1}\nSeq 2: {aligned_seq2}")

    st.markdown("---")

    # --- Score Matrix Visualization ---
    st.subheader("📊 Alignment Matrix (Score)")

    # Prepare DataFrame for score matrix visualization
    matrix_df = pd.DataFrame(matrix, index=["-"] + list(seq_a), columns=["-"] + list(seq_b))
    
    # Function to highlight the optimal path cells
    def highlight_path(s):
        """Highlights cells that are part of the optimal alignment path."""
        is_path = [(i, j) in path for i in range(matrix_df.shape[0]) for j in range(matrix_df.shape[1])]
        
        # Reshape the boolean list to match the DataFrame shape
        is_path_df = np.array(is_path).reshape(matrix_df.shape)
        
        # Apply style: yellow background for path, blue for the max score cell (local)
        styles = []
        for i in range(matrix_df.shape[0]):
            row_styles = []
            for j in range(matrix_df.shape[1]):
                style = ''
                if is_path_df[i, j]:
                    style = 'background-color: yellow; color: black;'
                
                # Highlight the maximum score cell differently for local alignment
                if is_local and i == max_pos[0] and j == max_pos[1]:
                    style = 'background-color: blue; color: white; font-weight: bold;'
                
                row_styles.append(style)
            styles.append(row_styles)
        
        return styles

    # Display the styled matrix
    st.dataframe(
        matrix_df.style.apply(highlight_path, axis=None), 
        use_container_width=True
    )
    
    st.markdown("---")

    # --- Traceback Path Visualization ---
    st.subheader("↩️ Traceback Path")
    st.info("The path highlights the cells used to derive the optimal alignment. For Local Alignment, it starts at the maximum score and stops at a cell with score 0.")

    # Prepare DataFrame for traceback visualization (showing direction codes)
    traceback_df = pd.DataFrame(traceback, index=["-"] + list(seq_a), columns=["-"] + list(seq_b))
    
    # Replace codes with directional symbols for better visualization
    # 0: Stop/End, 1: ↖ (Diagonal), 2: ↑ (Up), 3: ← (Left)
    traceback_symbols = traceback_df.replace({
        0: '◼️',  # Stop
        1: '↖️',  # Diagonal
        2: '⬆️',  # Up
        3: '⬅️'   # Left
    })

    st.dataframe(
        traceback_symbols.style.apply(highlight_path, axis=None), # Reuse path highlighting
        use_container_width=True
    )
    
    st.caption("Key: ↖️: Diagonal, ⬆️: Up, ⬅️: Left, ◼️: Stop (Local only)")

except Exception as e:
    st.error(f"An error occurred during alignment: {e}")

# --- End of Streamlit Code ---