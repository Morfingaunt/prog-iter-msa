# Progressive and Iterative MSA

### Preparing
1. First of all it's required to comment this line in Class TabularMSA (skbio/alignment/_tabular_msa.py):
>self._assert_valid_sequences(sequences)

It's required to add the support of all ASCII symbols.

2. Install IAB library :
>https://github.com/lsl5/An-Introduction-To-Applied-Bioinformatics

### Launching
>./perf.py

or throuth any other Python interpreter:
>python3.x perf.py