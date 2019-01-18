# Hierarchical Bayes SDT models

## Example Data:
[cidrs18_fulldatascreened.csv](cidrs18_fulldatascreened.csv) is an unpublished data set containing trial-level test phase data for participants in a joint priming/recognition/source experiment. At study, participants are shown nouns on a red or blue background, and hear a sentence about the word and the background colour. In one (within-subjects) condition, the sentence describes the object designated by the noun as being red or blue (Condition = 1, Item-condition), in the other the sentence describes the object interacting with a red/blue object (Condition = 2, Context-condition).

## Univariate
### (Univeriate/UVSDT_participant.stan)

Univariate Signal/Noise SDT, with participant effects