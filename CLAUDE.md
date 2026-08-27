\# AfroTB Phylogeny-Aware GNN Project



\## Project goal



Develop a phylogeny-aware Graph Neural Network for joint prediction of

tuberculosis drug resistance and Mycobacterium tuberculosis lineage using

the Afro-TB dataset.



The project focuses on African TB genomic diversity and evaluates whether

incorporating evolutionary relationships between isolates improves

drug-resistance prediction compared with conventional machine-learning

approaches.



\## Current development stage



This is a university undergraduate research project.



Review 2 requires approximately 60–75% implementation and a working

demonstration.



The implementation should prioritize:

1\. Correctness

2\. Reproducibility

3\. Biological validity

4\. Clear modular design

5\. Efficient computation

6\. Minimal unnecessary dependencies



\## Team responsibilities



\### Person 1

Owns:

1\. Afro-TB dataset acquisition and inspection

2\. VCF preprocessing

3\. Resistance-associated mutation feature extraction

4\. AMR labels

5\. Lineage labels

6\. Data-quality validation

7\. SNP-distance / phylogenetic representation

8\. Phylogenetic graph construction



\### Person 2

Owns:

1\. Graph Neural Network architecture

2\. Multi-task learning

3\. Model training

4\. Validation

5\. Model evaluation



\### Person 3

Owns:

1\. Random Forest baseline

2\. XGBoost baseline

3\. Global/non-African comparison where data permits

4\. Evaluation metrics

5\. Explainability

6\. Result visualization



\## Important biological rules



Never invent:

1\. Mutation positions

2\. Resistance associations

3\. Lineage assignments

4\. Dataset statistics

5\. Biological conclusions



Use the verified Afro-TB dataset and published literature as the source of

biological information.



If information is missing or ambiguous, report it instead of guessing.



\## Data contract



The final processed dataset should provide:



X\_mutations

N × F mutation feature matrix



y\_amr

N × D multi-label drug-resistance matrix



y\_lineage

N-dimensional lineage labels



edge\_index

2 × E graph connectivity matrix



edge\_weight

E-dimensional graph edge weights



splits

Fixed train/validation/test sample IDs shared by all models.



The exact value of F must be determined from the verified dataset rather

than assumed.



\## Development rules



1\. Do not process all 13,753 genomes during early development.

2\. First test on a very small subset.

3\. Then test on approximately 100 samples.

4\. Scale to the complete dataset only after validation.

5\. Do not modify unrelated files.

6\. Do not rewrite working code unnecessarily.

7\. Do not install packages without explaining why they are required.

8\. Prefer simple and reproducible implementations.

9\. Preserve sample IDs throughout every preprocessing stage.

10\. Keep train/validation/test splits fixed and reproducible.

11\. Never commit raw datasets or large generated files to Git.

12\. Never commit credentials, API keys, or .env files.



\## Person 1 output requirements



Person 1 should eventually provide:



1\. Clean sample metadata

2\. Mutation feature matrix

3\. AMR labels

4\. Lineage labels

5\. Dataset quality report

6\. Fixed train/validation/test splits

7\. Phylogenetic/SNP-distance representation

8\. Graph edge\_index

9\. Graph edge weights

10\. Documentation describing how each was generated



\## Coding style



Use clear descriptive variable and function names.



Prefer readable scientific Python over unnecessarily clever code.



Every important preprocessing function should have a concise docstring.



Avoid unnecessary abstractions.



Before implementing a major component:

1\. Inspect the existing files.

2\. State the implementation plan briefly.

3\. Implement only the requested component.

4\. Test it on a small dataset.

5\. Report the result.



Do not automatically proceed to the next project component.



\## Claude Code behaviour



When asked to modify code:



1\. Inspect only the relevant files first.

2\. Do not scan the entire repository unless necessary.

3\. Do not modify unrelated files.

4\. Do not create duplicate implementations.

5\. Do not generate large amounts of explanatory text.

6\. Keep responses concise and implementation-focused.

7\. Ask for clarification when biological assumptions are uncertain.

8\. Never fabricate missing data.

