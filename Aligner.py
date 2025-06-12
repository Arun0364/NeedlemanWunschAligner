# Import necessary libraries
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.Align import substitution_matrices


# Function to save aligned sequences in FASTA format
def save_alignment_fasta(aln1, aln2, name1="Seq1", name2="Seq2", file_name="alignment.fasta"):
    with open(file_name, "w") as f:
        f.write(f">{name1}\n{aln1}\n")
        f.write(f">{name2}\n{aln2}\n")
    print(f"\n✅ Alignment saved to {file_name}")

# Needleman-Wunsch global alignment using BLOSUM62 scoring
def needleman_wunsch_blosum(seq1, seq2, name1="Seq1", name2="Seq2", gap=-4):
	# Convert sequences to uppercase
	seq1 = seq1.upper()
	seq2 = seq2.upper()
	n,m = len(seq1), len(seq2)

	#Initializing scoring and traceback matrices
	score_matrix = np.zeros((n+1, m+1), dtype=int)
	traceback = np.zeros((n+1, m+1), dtype=str)

	#Loading BLOSUM62 from biopython
	blosum62 = substitution_matrices.load("BLOSUM62")	

	# Initialize the first row and column with gap penalties
	for i in range(1, n+1):
		score_matrix[i][0] = score_matrix[i-1][0] + gap
		traceback[i][0] = '↑' # Upward move
	for j in range(1, m+1):
		score_matrix[0][j] = score_matrix[0][j-1] + gap
		traceback[0][j] = '←' # Leftward move

	# Fill scoring and traceback matrices
	for i in range(1, n+1):
		for j in range(1, m+1):
			match = blosum62[seq1[i-1], seq2[j-1]] # Match/mismatch score
			diag = score_matrix[i-1][j-1] + match  # Diagonal move
			up = score_matrix[i-1][j] + gap        # Gap in seq2
			left = score_matrix[i][j-1] + gap      # Gap in seq1

			# Choose maximum score direction
			max_score = max(diag, up, left)
			score_matrix[i][j] = max_score
			
			if max_score == diag:
				traceback[i][j] = '↖'
			elif max_score == up:
				traceback[i][j] = '↑'
			else:
				traceback[i][j] = '←'
	
	# Traceback to build the aligned sequences
	align1 = ""
	align2 = ""
	i,j = n, m
	path_coords = [] # To store traceback path for heatmap
	
	while i > 0 or j > 0:
		path_coords.append((i, j))  # Save path coordinate
		direction = traceback[i][j]
		if direction == '↖':
			align1 = seq1[i-1] + align1
			align2 = seq2[j-1] + align2
			i -= 1
			j -= 1
		elif direction == '↑':
			align1 = seq1[i-1] + align1
			align2 = '-' + align2
			i -= 1
		elif direction == '←':
			align1 = '-' + align1
			align2 = seq2[j-1] + align2
			j -= 1


	path_coords.append((0, 0))     # Include the origin point
	path_coords = set(path_coords) # Convert to set for unique dots
			
	# Create a DataFrame from score matrix for heatmap
	df = pd.DataFrame(score_matrix, index = ["-"] + list(seq1), columns = ["-"] + list(seq2))

	# Plot heatmap of alignment matrix
	plt.figure(figsize = (10,8))
	ax = sns.heatmap(df, annot = True, fmt = "d", cmap = "YlGnBu", linewidth = 0.5, square = True, cbar=False)

	# Overlay red circles on traceback path
	for (i, j) in path_coords:
		ax.add_patch(plt.Circle((j + 0.5, i + 0.5), 0.15, color='red', fill=True))	

	# Customize plot
	plt.title("Needleman-Wunsch Alignment Score Matrix (BLOSUM62)")
	plt.xlabel("Sequence 2")
	plt.ylabel("Sequence 1")
	plt.tight_layout()

	# Save the heatmap to a PNG file
	plt.savefig("alignment_heatmap.png", dpi=300)
	plt.show()
	
	# Print results
	print("\nFinal Alignment Score:", score_matrix[n][m])
	print("\nAligned Sequences:")
	print("Seq1:", align1)
	print("     ", ''.join('|' if a == b else ' ' for a, b in zip(align1, align2)))
	print("Seq2:", align2)

	# Save aligned sequences to a FASTA file
	save_alignment_fasta(align1, align2, name1, name2)


	return align1, align2, score_matrix[n][m]

# ======= MAIN SCRIPT EXECUTION ======= #

# User inputs for sequence names and sequences
name1 = input("Enter name for first sequence: ")
seq1 = input("Enter first sequence: ")
name2 = input("Enter name for second sequence: ")
seq2 = input("Enter second sequence: ")

# Run alignment
needleman_wunsch_blosum(seq1, seq2, name1, name2)