# Hierarchical Bayes SDT models

## Example Data:
[cidrs18_fulldatascreened.csv](cidrs18_fulldatascreened.csv) is an unpublished data set containing trial-level test phase data for participants in a joint priming/recognition/source experiment. At study, participants are shown nouns on a red or blue background, and hear a sentence about the word and the background colour. In one (within-subjects) condition, the sentence describes the object designated by the noun as being red or blue (Condition = 1, Item-condition), in the other the sentence describes the object interacting with a red/blue object (Condition = 2, Context-condition).

## Univariate
### [UVSDT_participant](Univariate/UVSDT_participant.stan)

Univariate Signal/Noise UV SDT with participant effects, reasonably uninformative priors

### [EVSDT_participant](Univariate/EVSDT_participant.stan)

Equal variance version of above, assuming 1 participant specific sigma for signal and noise items

### [DP_EVSDT_Source2R_participant](Univariate/DP_EVSDT_Source2R_participant.stan)

Specific version of Yonelinas dual-process model -> for source memory, recollection for signal and noise items (Source A v Source B), and Rs and Rn estimated separately for signal and noise items (written for unitization literature exploration, Diana et al, 2008)), participant effects

## To Play with:

* item effects, especially for unitization, probably less of an issue for standard word list but since the data is there probably good to explore full model
* parameterize as theta_s = theta_n + theta_d' and estimate d' and alpha' instead [https://doi.org/10.1037/a0028144](https://doi.org/10.1037/a0028144)
* Potentially not that trivial re: working around the lack of multi_normal_cdf:
	* bivariate recognition/source
	* multivariate recognition/source/priming (cidrs task) -- maybe start with the non-hierarchical version here
		* bivariate Priming-Recognition or Priming-Source should work more easily out of the box