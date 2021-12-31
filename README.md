# Aligning sequences

## Implement the Smith-Waterman algorithm to align two sequences.

For this assignment, implement the `smith_waterman` function in align_sequences.py. The function takes in two nucleotide sequences (seq1 and seq2), and scores/penalty for a match, mismatch, opening of a new gap, and extension of a gap. You should write the `smith_waterman` function to calculate the maximum score of the local alignment and the resulting aligned subsequences from the two input sequences. For example if seq1 is `ACGTTA` and seq2 is `ACTTA`, match=1, mismatch=1, gapopen=1, gapextend=0, then the expected output of the function will be `4, 'ACGTTA','AC-TTA'`

## Testing

You can test your work by running:

```
python3 testdriver
```
