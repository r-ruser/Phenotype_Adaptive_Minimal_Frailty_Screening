## Methods for manuscript use

We compared the locked SHAP-selected six-item score with two item-selection benchmarks. For the simple-correlation benchmark, the absolute Spearman correlation between each of the 16 candidate items and continuous complete FI was calculated among participants in the four outer-training cohorts using pairwise observed data. The six highest-ranked items were selected, with the candidate-bank order resolving ties. Each selected score was an unweighted mean of its observed items and required at least five observed items. The held-out cohort contributed exclusively to evaluation.

For the random benchmark, a fixed seed (20261001) selected 1,000 distinct combinations uniformly without replacement from the 8,008 possible six-item combinations. These same 1,000 combinations were evaluated in each of the five held-out cohorts. A shared risk set required all original universal and phenotype-specific FI1–FI16 scores to be calculable and at least 15 of the 16 candidate items to be observed. This rule ensured at least five observed items for every possible six-item combination. Outcome-specific comparisons additionally required the relevant outcome to be observed. All three item-selection approaches were evaluated on the same participants within each cohort and outcome.

Concurrent evaluation included AUROC for complete FI ≥0.25, Spearman correlation with continuous complete FI, and MAE and RMSE of the raw six-item mean relative to continuous complete FI. Prospective evaluation included AUROC for follow-up FI ≥0.25 and new ADL limitation. Original stored prospective outcomes retained the original eligibility and ascertainment rules. Correlation and prediction metrics quantify ranking and discrimination; MAE and RMSE quantify numerical agreement on the raw score scales.

SHAP-minus-simple-correlation differences used 500 paired participant-bootstrap samples per cohort and outcome. Identical equal-weight scores yield zero paired differences on every resample. Random-set results are summarized by the median and 2.5th–97.5th percentiles across the 1,000 combinations. These percentiles describe variation across item combinations. SHAP empirical performance percentiles use the fraction of random sets with poorer performance, assigning half weight to ties, and orient MAE/RMSE so that lower errors represent better performance. Cohort-mean summaries average the five cohort-specific metrics with equal cohort weights for each item combination.

The six-item benchmark comparison belongs to the descriptive, exploratory length analysis. Both trained rankings use outer-training data; the prominence of the candidate length reflects review of the held-out length curves.

## Main quantitative findings

The simple-correlation top-six set matched the SHAP top-six set in all five outer folds: lifting, stair climbing, self-rated health, chair rise, stooping/kneeling/crouching, and grocery shopping. Item order differed between approaches, while their equal-weight six-item scores were identical. All 30 evaluated paired metric differences and paired-bootstrap intervals were 0.

The fair comparison samples comprised 34,126 concurrent participants, 21,818 participants with observed follow-up FI outcomes, and 18,472 participants with new-ADL outcomes. Relative to the main analysis risk sets, the candidate-observation rule removed 76, 43, and 31 participants, respectively.

| Metric, equally weighted cohort mean | SHAP and simple correlation | Random-set median | Random-set 2.5th–97.5th percentiles | SHAP better-than-random percentile |
| --- | ---: | ---: | ---: | ---: |
| Concurrent AUROC | 0.9388 | 0.8924 | 0.7611–0.9268 | 100.0% |
| Concurrent Spearman | 0.8692 | 0.7967 | 0.6534–0.8487 | 100.0% |
| Raw-score MAE | 0.1787 | 0.1193 | 0.0869–0.1877 | 5.4% |
| Raw-score RMSE | 0.2293 | 0.1460 | 0.1131–0.2156 | 0.3% |
| Follow-up FI AUROC | 0.8242 | 0.7758 | 0.6598–0.8152 | 100.0% |
| New ADL AUROC | 0.7060 | 0.6646 | 0.5571–0.6972 | 99.9% |

### Interpretation for manuscript use

The trained six-item set showed stronger ranking and discrimination than most sampled six-item combinations. Simple correlation recovered the same item set and the same equal-weight score. Thus, these results support the stability and utility of the six-item candidate set and quantify an observed performance gain of zero for SHAP relative to this simple selection benchmark.

Raw six-item means displayed larger numerical errors than typical random combinations, alongside stronger discrimination. This pattern highlights the distinction between ordering participants and matching the numerical complete-FI scale. The raw six-item mean requires its own scale interpretation and training-based threshold or mapping. The complete-FI 0.25 threshold defines the reference classification in this analysis; the original locked six-item cutoffs govern six-item classification.

