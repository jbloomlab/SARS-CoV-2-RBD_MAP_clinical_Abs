# Analyze naturally occurring mutations at sites of strong escape
This Python Jupyter notebook sees how many naturally occuring mutations are observed at each site of strong escape

## Set up analysis
Import Python modules:


```python
import math
import os

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import pandas as pd

from plotnine import *

import yaml
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Read escape profiles config, which tells which sets to make plots for:


```python
with open(config['escape_profiles_config']) as f:
    escape_profiles_config = yaml.safe_load(f)
```

Create output directory:


```python
os.makedirs(config['gisaid_mutations_dir'], exist_ok=True)
```

Read counts of naturally ocurring mutations:


```python
print(f"Reading mutation counts from {config['gisaid_mutation_counts']}")

mut_counts = pd.read_csv(config['gisaid_mutation_counts'])
```

    Reading mutation counts from results/GISAID_mutations/mutation_counts.csv


Read sites of "strong escape" from all antibodies / sera:


```python
print(f"Reading sites of strong escape from {config['strong_escape_sites']}")

strong_sites = pd.read_csv(config['strong_escape_sites'])
```

    Reading sites of strong escape from results/escape_profiles/strong_escape_sites.csv


Read escape fractions for all antibodies / sera:


```python
print(f"Reading escape fractions from {config['escape_fracs']}")

escape_fracs = (
    pd.read_csv(config['escape_fracs'])
    .query('library == "average"')
    .drop(columns='site')
    .rename(columns={'mutation': 'mutant',
                     'label_site': 'site'})
    [['condition', 'site', 'wildtype', 'mutant', config['mut_metric'], config['site_metric']]]
    )

escape_fracs
```

    Reading escape fractions from results/escape_scores/escape_fracs.csv





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>condition</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutant</th>
      <th>mut_escape_frac_epistasis_model</th>
      <th>site_total_escape_frac_epistasis_model</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CB6_400</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>0.001093</td>
      <td>0.04032</td>
    </tr>
    <tr>
      <th>1</th>
      <td>CB6_400</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>0.002187</td>
      <td>0.04032</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CB6_400</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>0.001093</td>
      <td>0.04032</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CB6_400</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>0.001093</td>
      <td>0.04032</td>
    </tr>
    <tr>
      <th>4</th>
      <td>CB6_400</td>
      <td>331</td>
      <td>N</td>
      <td>G</td>
      <td>0.005493</td>
      <td>0.04032</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>7787</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>R</td>
      <td>0.002079</td>
      <td>0.03943</td>
    </tr>
    <tr>
      <th>7788</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>S</td>
      <td>0.002079</td>
      <td>0.03943</td>
    </tr>
    <tr>
      <th>7789</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>V</td>
      <td>0.002076</td>
      <td>0.03943</td>
    </tr>
    <tr>
      <th>7790</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>W</td>
      <td>0.002072</td>
      <td>0.03943</td>
    </tr>
    <tr>
      <th>7791</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>Y</td>
      <td>0.002079</td>
      <td>0.03943</td>
    </tr>
  </tbody>
</table>
<p>7792 rows × 6 columns</p>
</div>



## Counts of mutations at sites of escape
Get counts of naturally occurring mutations at sites of escape, along with the actual escape values:

First get mutation-level counts:


```python
mutcounts_strong_sites = (
    strong_sites[['condition', 'threshold', 'site']]
    .merge(mut_counts, how='inner', on='site')
    .merge(escape_fracs[['condition', 'site', 'wildtype', config['site_metric']]].drop_duplicates(),
           on=['condition', 'site', 'wildtype'],
           validate='many_to_one')
    .assign(mutation=lambda x: x['wildtype'] + x['site'].astype(str) + x['mutant'])
    .sort_values('count', ascending=False)
    )
```

Now get site-level counts (aggregating all mutations at a site):


```python
sitecounts_strong_sites = (
    mutcounts_strong_sites
    .assign(mut_count=lambda x: x['mutation'] + ' (' + x['count'].astype(str) + ')')
    .groupby(['condition', 'threshold', 'site', 'wildtype', config['site_metric']])
    .aggregate({'count': 'sum', 'mut_count': ', '.join})
    .rename(columns={'mut_count': 'counts_by_mutation'})
    .reset_index()
    .sort_values('count', ascending=False)
    )

print(f"Here are first few lines showing the most frequently mutated sites of escape:")
display(HTML(sitecounts_strong_sites.head(n=20).to_html(index=False)))
```

    Here are first few lines showing the most frequently mutated sites of escape:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>site_total_escape_frac_epistasis_model</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>CB6_400</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>1.4950</td>
      <td>11434</td>
      <td>N501Y (11216), N501T (209), N501S (6), N501I (3)</td>
    </tr>
    <tr>
      <td>REGN10987_400</td>
      <td>sensitive</td>
      <td>439</td>
      <td>N</td>
      <td>2.2040</td>
      <td>5262</td>
      <td>N439K (5262)</td>
    </tr>
    <tr>
      <td>REGN10987_400</td>
      <td>default</td>
      <td>439</td>
      <td>N</td>
      <td>2.2040</td>
      <td>5262</td>
      <td>N439K (5262)</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>default</td>
      <td>453</td>
      <td>Y</td>
      <td>4.0980</td>
      <td>994</td>
      <td>Y453F (994)</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>sensitive</td>
      <td>453</td>
      <td>Y</td>
      <td>4.0980</td>
      <td>994</td>
      <td>Y453F (994)</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>sensitive</td>
      <td>484</td>
      <td>E</td>
      <td>3.0380</td>
      <td>615</td>
      <td>E484K (536), E484Q (54), E484A (8), E484G (7), E484D (6), E484R (3), E484V (1)</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>3.0380</td>
      <td>615</td>
      <td>E484K (536), E484Q (54), E484A (8), E484G (7), E484D (6), E484R (3), E484V (1)</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>sensitive</td>
      <td>484</td>
      <td>E</td>
      <td>0.9369</td>
      <td>615</td>
      <td>E484K (536), E484Q (54), E484A (8), E484G (7), E484D (6), E484R (3), E484V (1)</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>sensitive</td>
      <td>417</td>
      <td>K</td>
      <td>17.8300</td>
      <td>349</td>
      <td>K417N (347), K417R (2)</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>default</td>
      <td>417</td>
      <td>K</td>
      <td>4.5550</td>
      <td>349</td>
      <td>K417N (347), K417R (2)</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>default</td>
      <td>417</td>
      <td>K</td>
      <td>17.8300</td>
      <td>349</td>
      <td>K417N (347), K417R (2)</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>sensitive</td>
      <td>417</td>
      <td>K</td>
      <td>4.5550</td>
      <td>349</td>
      <td>K417N (347), K417R (2)</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>default</td>
      <td>475</td>
      <td>A</td>
      <td>10.2900</td>
      <td>110</td>
      <td>A475S (59), A475V (50), A475T (1)</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>sensitive</td>
      <td>475</td>
      <td>A</td>
      <td>10.2900</td>
      <td>110</td>
      <td>A475S (59), A475V (50), A475T (1)</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>default</td>
      <td>475</td>
      <td>A</td>
      <td>1.8420</td>
      <td>110</td>
      <td>A475S (59), A475V (50), A475T (1)</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>sensitive</td>
      <td>475</td>
      <td>A</td>
      <td>1.8420</td>
      <td>110</td>
      <td>A475S (59), A475V (50), A475T (1)</td>
    </tr>
    <tr>
      <td>REGN10933+REGN10987_400</td>
      <td>sensitive</td>
      <td>446</td>
      <td>G</td>
      <td>0.1137</td>
      <td>101</td>
      <td>G446V (90), G446S (8), G446A (2), G446D (1)</td>
    </tr>
    <tr>
      <td>REGN10987_400</td>
      <td>default</td>
      <td>446</td>
      <td>G</td>
      <td>17.1200</td>
      <td>101</td>
      <td>G446V (90), G446S (8), G446A (2), G446D (1)</td>
    </tr>
    <tr>
      <td>REGN10987_400</td>
      <td>sensitive</td>
      <td>446</td>
      <td>G</td>
      <td>17.1200</td>
      <td>101</td>
      <td>G446V (90), G446S (8), G446A (2), G446D (1)</td>
    </tr>
    <tr>
      <td>REGN10987_400</td>
      <td>sensitive</td>
      <td>440</td>
      <td>N</td>
      <td>6.6100</td>
      <td>70</td>
      <td>N440K (53), N440Y (11), N440D (3), N440S (2), N440T (1)</td>
    </tr>
  </tbody>
</table>


Now plot mutation counts (any mutation) at each site of escape for each antibody / sera:


```python
nconditions = sitecounts_strong_sites['condition'].nunique() * sitecounts_strong_sites['threshold'].nunique()
ncol = 8
nrow = math.ceil(nconditions / ncol)

p = (ggplot(sitecounts_strong_sites) +
     aes(config['site_metric'], 'count') +
     geom_point(alpha=0.5) +
     facet_wrap('~ condition + threshold', ncol=ncol) +
     theme(figure_size=(2 * ncol, 2 * nrow)) +
     xlab('site-level escape') +
     ylab('sequences with mutation at site')
     )

_ = p.draw()
```


    
![png](natural_mutations_files/natural_mutations_20_0.png)
    


## Perform analyses on subsets
We perform analyses on all subsets in the escape profiles config for which this is specified:


```python
for name, specs in escape_profiles_config.items():
    if 'analyze_natural_mutations' not in specs or not specs['analyze_natural_mutations']:
        continue
    print(f"\nAnalyzing natural mutations for {name}")
    
    conditions = specs['conditions']
    
    threshold = specs['plot_auto_identified_sites']
    if threshold not in sitecounts_strong_sites['threshold'].unique():
        raise ValueError(f"invalid threshold {threshold} for {name}")
    
    # get count for conditions of interest for this subset
    df = (sitecounts_strong_sites
          .query('condition in @conditions')
          .query('threshold == @threshold')
          .assign(condition=lambda x: x['condition'].map(conditions))
          .drop(columns=config['site_metric'])
          )
    countsfile = os.path.join(config['gisaid_mutations_dir'], f"{name}_mutation_counts.csv")
    print(f"Writing counts of mutations at sites of strong escape to {countsfile}. First few lines:")
    display(HTML(df.head(n=10).to_html(index=False)))
    df.to_csv(countsfile, index=False)
    
    # make plot showing escape sites with more than mincounts mutations
    if 'natural_mutations_mincounts' in specs:
        mincounts = specs['natural_mutations_mincounts']
    else:
        mincounts = 5
    plotsfile = os.path.join(config['gisaid_mutations_dir'], f"{name}_mutation_counts.pdf")
    print('Plotting which antibodies / sera are escaped by mutations at all sites of '
          f"escape with at least {mincounts} mutation counts and saving to {plotsfile}.")
    plot_df = (
        # data frame with all combinations of conditions and sites
        pd.DataFrame.from_records([(condition, site) for condition in conditions.values()
                                   for site in df['site'].unique()],
                                  columns=['condition', 'site'])
        # annotate sites of escape
        .merge(df.assign(escape=lambda x: x['count'] >= mincounts)[['condition', 'site', 'escape']],
               how='left',
               validate='one_to_one',
               on=['condition', 'site'])
        .assign(escape=lambda x: x['escape'].fillna(False))
        # add wildtype and counts of mutations at each site
        .merge(sitecounts_strong_sites[['site', 'wildtype', 'count']].drop_duplicates(),
               how='left',
               validate='many_to_one',
               on='site')
        # get only sites with sufficient mutation counts
        .query('count > @mincounts')
        # only get sites where at least one antibody escapes
        .assign(n_escape=lambda x: x.groupby('site')['escape'].transform('sum'))
        .query('n_escape > 0')
        # order conditions, and order sites by count after making nice label
        .assign(site_label=lambda x: x['wildtype'] + x['site'].astype(str) + ' (' + x['count'].astype(str) + ')')
        .sort_values('count')
        .assign(condition=lambda x: pd.Categorical(x['condition'], list(conditions.values()), ordered=True),
                site_label=lambda x: pd.Categorical(x['site_label'], x['site_label'].unique(), ordered=True)
                )
        )
    p = (ggplot(plot_df) +
         aes('condition', 'site_label', fill='escape') +
         geom_tile(color='black', size=0.3) +
         theme(axis_text_x=element_text(angle=90),
               figure_size=(0.3 * plot_df['condition'].nunique(), 0.3 * plot_df['site_label'].nunique()),
               panel_background=element_blank(),
               ) +
         xlab('') +
         ylab('') +
         scale_fill_manual(values=['white', 'dimgray'])
         )
    p.save(plotsfile, verbose=False)
    fig = p.draw()
    display(fig)
    plt.close(fig)
```

    
    Analyzing natural mutations for REGN_and_LY-CoV016
    Writing counts of mutations at sites of strong escape to results/GISAID_mutations/REGN_and_LY-CoV016_mutation_counts.csv. First few lines:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>LY-CoV016</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>11434</td>
      <td>N501Y (11216), N501T (209), N501S (6), N501I (3)</td>
    </tr>
    <tr>
      <td>REGN10987</td>
      <td>sensitive</td>
      <td>439</td>
      <td>N</td>
      <td>5262</td>
      <td>N439K (5262)</td>
    </tr>
    <tr>
      <td>REGN10933</td>
      <td>sensitive</td>
      <td>453</td>
      <td>Y</td>
      <td>994</td>
      <td>Y453F (994)</td>
    </tr>
    <tr>
      <td>REGN10933</td>
      <td>sensitive</td>
      <td>484</td>
      <td>E</td>
      <td>615</td>
      <td>E484K (536), E484Q (54), E484A (8), E484G (7), E484D (6), E484R (3), E484V (1)</td>
    </tr>
    <tr>
      <td>LY-CoV016</td>
      <td>sensitive</td>
      <td>484</td>
      <td>E</td>
      <td>615</td>
      <td>E484K (536), E484Q (54), E484A (8), E484G (7), E484D (6), E484R (3), E484V (1)</td>
    </tr>
    <tr>
      <td>LY-CoV016</td>
      <td>sensitive</td>
      <td>417</td>
      <td>K</td>
      <td>349</td>
      <td>K417N (347), K417R (2)</td>
    </tr>
    <tr>
      <td>REGN10933</td>
      <td>sensitive</td>
      <td>417</td>
      <td>K</td>
      <td>349</td>
      <td>K417N (347), K417R (2)</td>
    </tr>
    <tr>
      <td>LY-CoV016</td>
      <td>sensitive</td>
      <td>475</td>
      <td>A</td>
      <td>110</td>
      <td>A475S (59), A475V (50), A475T (1)</td>
    </tr>
    <tr>
      <td>REGN10933</td>
      <td>sensitive</td>
      <td>475</td>
      <td>A</td>
      <td>110</td>
      <td>A475S (59), A475V (50), A475T (1)</td>
    </tr>
    <tr>
      <td>REGN10933 + REGN10987</td>
      <td>sensitive</td>
      <td>446</td>
      <td>G</td>
      <td>101</td>
      <td>G446V (90), G446S (8), G446A (2), G446D (1)</td>
    </tr>
  </tbody>
</table>


    Plotting which antibodies / sera are escaped by mutations at all sites of escape with at least 5 mutation counts and saving to results/GISAID_mutations/REGN_and_LY-CoV016_mutation_counts.pdf.



    
![png](natural_mutations_files/natural_mutations_22_3.png)
    



```python

```
