% My project involves calculating the folding energies of RNA sequences of
% the a specific part of the trp operon.  I also want to see the folding
% energy when you make any number of changes of nucleic acid base sequence
% of part of the trp operon sequence of E. Coli (Yanofsky et al., 1981). I
% specificall chose to focus on the trpA gene. Potentially, I’d like to show
% the location and position of the mutated structures and energies obtained
% then examine and compare the structures of varying sequences. Then, I plan
% on showing the diagrams of each of the structures obtained and the
% energies. Next, we can see the base changes you made and how such changes
% would affect the transcription of the downstream genes if those changes
% were introduced into the genome of the trp operon of E. coli. We can also
% view the energy of stems in terms and the effect on in vivo function.

% First I compare mutations of the original sequence and a random
% permutation of mutations invoked into the original trpA sequence. The
% enzyme this trpA gene codes for is tryptophan synthase subunit ? from the
% Escherichia coli K-12 substr. MG1655

trpA = fastaread('nuc-seq-file-EG11024.txt.alter.fas'); %read file into system
[trpAHeader, trpASequenceDNA] = fastaread('nuc-seq-file-EG11024.txt.alter.fas') %repeated reading file into system, labels trpA header and trpA Sequence
trpASequenceRNA = dna2rna(trpASequenceDNA) %convert DNA sequence to RNA

trpAseqMut = trpASequenceRNA
point_mutations = 20
length_of_trpA = length(trpASequenceRNA)
s = RandStream('mlfg6331_64')
trpAseqMut(randi(length_of_trpA,1,point_mutations)) = datasample(s,'ACGU',point_mutations,'Weights',[0.25 0.25 0.25 0.25])

% performed an alignment of code sequences for the original trpA RNA code
% and its mutation
[trpA_Score, trpA_Alignment] = swalign(trpASequenceRNA,trpAseqMut,'Alphabet','NT')
trpA_length = length(trpA_Alignment)
trpA_Array = trpA_Alignment(2,:)
trpA_matches = count(trpA_Array,'|')
Average_trpA_Reads = trpA_matches/trpA_length

% Made a scoring matrix with the original trpA RNA code and its mutation to
% compare and observe optimal alignment
matchval = 2;
mismatchval = -1;
ofdiag = ones(4)-eye(4);
S = matchval*eye(4)+mismatchval*ofdiag
figure;
[score, align, start] = swalign(trpASequenceRNA,trpAseqMut,'alphabet','nt','ScoringMatrix',S,'GapOpen',1,'Showscore',true);
align

% Second, I plan on calculating the folding energies of the trp operon

[trpA_bracket, trpA_energy, trpA_matrix] = rnafold(trpASequenceRNA); 
%returns RNAmatrix, a connectivity matrix representing the secondary structure associated with the minimum free energy.
trpA_energy

[trpAMutation_bracket, trpAMutation_energy, trpAMutation_matrix] = rnafold(trpAseqMut);
trpAMutation_energy

% Third, I show the diagram of RNA folding 

trpAfold = rnafold(trpASequenceRNA);
%returns the minimum free energy secondary structure of the normal trpA gene in bracket notation. In the bracket notation, each dot represents an unpaired base, while a pair of equally nested, opening and closing brackets represents a base pair.
trpAMutationfold = rnafold(trpAseqMut);
%returns the minimum free energy secondary structure of the mutated trpA gene

rnaplot(trpAfold, 'sequence', trpASequenceRNA, 'format', 'mountain', 'colorby', 'pair');
hold on;
rnaplot(trpAMutationfold, 'sequence', trpAseqMut, 'format', 'mountain', 'colorby', 'pair'); hold off;


rnaplot(trpAfold, 'sequence', trpASequenceRNA, 'format', 'dotdiagram');
trpA_energy

rnaplot(trpAMutationfold, 'sequence', trpAseqMut, 'format', 'dotdiagram', 'selection', randi(length_of_trpA,1,point_mutations));
trpAMutation_energy

