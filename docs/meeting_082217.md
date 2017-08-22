### 082217

### Updates
* Fixed mutalyzer indels except for complex indels (ins+del)
* burden w/ pathogenic burden evaluation:
    * enrichmed: qval < .01
    * not enriched: qval > .2
    * ignored 'none' domains
    * lof, missense + inframe, total
    * lof is not enriched in pathogenic. Are there some gene/domains I can spot check? relaxing qval and domain does not help either missense nor lof
* can't dl mtr scores, will use missense-badness instead
