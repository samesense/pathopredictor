### 2018_05_18
* single gene fig, method, result
* background
* discussion
* fig1 and fig importance
* I don't see how clinvar single ran w/ the right data. need to account for clinvar set in main prediction data
* I can't seem to parse gnomad this morining.
* Clinvar single vs total was not a problem, but the script should have been clearer. It's under control now.
* No training data after removing vars. Try not using vest cutoff. Fathmm too. That's a lot of work. There are some left after the vest and fathmm filters.
* Also must redo clinvar restrictions for these guys. Do this only if vest and fathmm are removed.
* only epi has enough training vars

### 2018_05_17
* Mehods done except single gene
* single gene: single_gene_collapse_.003_.1.pdf. Need to filter gnomad first. Running. parse_gnomad like clinvar and panel. Don't forget to filter vars after w/ limit_eval_general. remove panel vars from extra gnomad vars, and remove extra gnomad vars from clinvar test
* what about the count figure? How is it made if global eval limits my genes? I should update fig1. It uses results, so it should be fine.

### 2018_05_10
* First tackle background
* Then methods