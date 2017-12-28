## Results

### Clinvar counts
* [total missense clinvar](https://github.com/samesense/mahdi_epi/blob/master/notebooks/clinvar-report.ipynb)

### Predict panel variant status

#### Questions
* Does training w/ clinvar do better than MPC>=2?
* Does training with the panel do better than MPC>=2?
* Does training with panel genes from clinvar do better than all of clinvar?
* Does training with the panel do better than training with clinvar?
* Do gene specific models do better than global models?

#### Global cutoff methods
* MPC>2 is pathogenic (baseline MPC cutoff)
* hold out test gene from panel and train
* train using all of clinvar minus testing variants
* train using clinvar panel genes minus testing variants

#### Single gene cutoff
* When genes have at least 5 pathogenic and 5 benign variants, evaluate the gene
* Use all methods from the global evaluation, but evaluation is limited by above statement
* Train using clinvar variants from this gene

#### Results
* Using the baseline MPC>=2 is always worse than?
