#TAY CHING XIAN A23CS0307
import streamlit as st
import numpy as np
import pandas as pd

def initialize_matrix(seq1_len, seq2_len, gap_penalty, is_local=False):
    rows = seq1_len + 1
    cols = seq2_len + 1
    matrix = np.zeros((rows, cols), dtype=int)
    traceback = np.zeros((rows, cols), dtype=int)

    if not is_local:
        for i in range(1, rows):
            matrix[i, 0] = matrix[i-1, 0] + gap_penalty
            traceback[i, 0] = 2
        for j in range(1, cols):
            matrix[0, j] = matrix[0, j-1] + gap_penalty
            traceback[0, j] = 3
    
    return matrix, traceback

def perform_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, is_local=False):
    m, n = len(seq1), len(seq2)
    matrix, traceback = initialize_matrix(m, n, gap_penalty, is_local)
    
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            
            current_score = match_score if seq1[i-1] == seq2[j-1] else mismatch_score
            score_diag = matrix[i-1, j-1] + current_score
            
            score_up = matrix[i-1, j] + gap_penalty
            
            score_left = matrix[i, j-1] + gap_penalty
            
            scores = [score_diag, score_up, score_left]
            
            if is_local:
                scores.append(0) 
            
            max_score_cell = max(scores)
            matrix[i, j] = max_score_cell

            if is_local and max_score_cell == 0:
                traceback[i, j] = 0
            elif max_score_cell == score_diag:
                traceback[i, j] = 1
            elif max_score_cell == score_up:
                traceback[i, j] = 2
            elif max_score_cell == score_left:
                traceback[i, j] = 3
            
            if is_local and matrix[i, j] > max_score:
                max_score = matrix[i, j]
                max_pos = (i, j)

    final_score = max_score if is_local else matrix[m, n]
    return final_score, matrix, traceback, max_pos

def traceback_alignment(seq1, seq2, matrix, traceback, max_pos, is_local=False):
    aligned_seq1 = ""
    aligned_seq2 = ""
    
    i, j = max_pos if is_local else (len(seq1), len(seq2))
    
    path = []
    
    while i > 0 or j > 0:
        path.append((i, j))
        
        if is_local and traceback[i, j] == 0:
            break
        
        move = traceback[i, j]
        
        if move == 1 and i > 0 and j > 0:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif move == 2 and i > 0:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        elif move == 3 and j > 0:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        elif i > 0:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        elif j > 0:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        else:
            break
    
    if i == 0 and j == 0:
        path.append((0, 0))
    
    return aligned_seq1, aligned_seq2, path[::-1]

st.set_page_config(layout="wide", page_title="Sequence Alignment Tool")

st.title("🧬 Dynamic Programming Sequence Alignment")
st.markdown("Implemention of **Needleman-Wunsch** (Global) and **Smith-Waterman** (Local) algorithms.")

with st.sidebar:
    st.header("⚙️ Alignment Settings")
    
    st.subheader("Sequences")
    seq_a = st.text_input("Sequence 1", value="ACGT", max_chars=50).upper().strip()
    seq_b = st.text_input("Sequence 2", value="CATG", max_chars=50).upper().strip()

    st.subheader("Alignment Type")
    alignment_type = st.radio("Choose Algorithm:", ["Global Alignment (Needleman-Wunsch)", "Local Alignment (Smith-Waterman)"])
    
    is_local = alignment_type == "Local Alignment (Smith-Waterman)"

    st.subheader("Scoring Metrics")
    match_score = st.number_input("Match Score", value=2, min_value=0, step=1)
    mismatch_score = st.number_input("Mismatch Score (must be < 0 for effective scoring)", value=-1, max_value=0, step=1)
    gap_penalty = st.number_input("Gap Penalty (must be < 0)", value=-2, max_value=0, step=1)

    if not seq_a or not seq_b:
        st.error("Please enter both sequences.")
        st.stop()

st.header(f"Results: {alignment_type}")
st.write(f"**Sequence 1:** `{seq_a}`")
st.write(f"**Sequence 2:** `{seq_b}`")
st.write(f"**Scoring:** Match: `{match_score}`, Mismatch: `{mismatch_score}`, Gap: `{gap_penalty}`")

try:
    score, matrix, traceback, max_pos = perform_alignment(
        seq_a, seq_b, match_score, mismatch_score, gap_penalty, is_local
    )

    aligned_seq1, aligned_seq2, path = traceback_alignment(
        seq_a, seq_b, matrix, traceback, max_pos, is_local
    )

    st.subheader("✅ Alignment Summary")
    
    if is_local:
        st.metric("Optimal Local Alignment Score", value=score, help="The highest score found in the matrix (Smith-Waterman).")
    else:
        st.metric("Optimal Global Alignment Score", value=score, help="The score in the bottom-right cell (Needleman-Wunsch).")
    
    st.markdown("#### Aligned Sequences")
    st.code(f"Seq 1: {aligned_seq1}\nSeq 2: {aligned_seq2}")

    st.markdown("---")

    st.subheader("📊 Alignment Matrix (Score)")

    matrix_df = pd.DataFrame(matrix, index=["-"] + list(seq_a), columns=["-"] + list(seq_b))
    
    def highlight_path(data):
        styles = pd.DataFrame('', index=data.index, columns=data.columns)
        
        for i in range(len(data.index)):
            for j in range(len(data.columns)):
                if (i, j) in path:
                    styles.iloc[i, j] = 'background-color: yellow; color: black;'
                
                if is_local and i == max_pos[0] and j == max_pos[1]:
                    styles.iloc[i, j] = 'background-color: blue; color: white; font-weight: bold;'
        
        return styles

    st.dataframe(
        matrix_df.style.apply(highlight_path, axis=None), 
        use_container_width=True
    )
    
    st.markdown("---")

    st.subheader("↩️ Traceback Path")
    st.info("The path highlights the cells used to derive the optimal alignment. For Local Alignment, it starts at the maximum score and stops at a cell with score 0.")

    traceback_df = pd.DataFrame(traceback, index=["-"] + list(seq_a), columns=["-"] + list(seq_b))
    
    traceback_symbols = traceback_df.replace({
        0: '◼️',
        1: '↖️',
        2: '⬆️',
        3: '⬅️'
    })

    st.dataframe(
        traceback_symbols.style.apply(highlight_path, axis=None),
        use_container_width=True
    )
    
    st.caption("Key: ↖️: Diagonal, ⬆️: Up, ⬅️: Left, ◼️: Stop (Local only)")

except Exception as e:
    st.error(f"An error occurred during alignment: {e}")
