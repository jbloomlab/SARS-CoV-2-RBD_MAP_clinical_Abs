# Output PDBs with escape scores as B factors
This Python Jupyter notebook outputs PDBs with the escape scores as B factors.

Though we will want more elaborate series of commands to codify our visualization of these RBD structures colored by escape, the series of commands below, when executed in a `PyMol` session with one of these PDBs open, will color the RBD surface according to escape scores.

For example, to normalize each structure colored by the max mut effect, we might want to have a white to red scale from 0 to 1:

     create RBD, chain E
     hide all; show cartoon, chain A; color gray20, chain A
     show surface, RBD; spectrum b, white red, RBD, minimum=0, maximum=1
     
For something like total escape, maybe we want each structure normalized to the maximum total escape in that structure, in which case we can just leave the maximum argument empty:

     create RBD, chain E
     hide all; show cartoon, chain A; color gray20, chain A
     show surface, RBD; spectrum b, white red, RBD, minimum=0
     
We write PDBs with B factors indicating the total site escape and maximum mutation escape at each site, and the same with these values normalized to the maximum for the full structure (the latter are easier to process in `Chimera`).

First, import Python modules:


```python
import collections
import copy
import os
import warnings

import Bio.PDB

import dms_variants.pdb_utils

from IPython.display import display, HTML

import pandas as pd

import yaml
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Read configuration for outputting PDBs:


```python
print(f"Reading PDB output configuration from {config['output_pdbs_config']}")
with open(config['output_pdbs_config']) as f:
    output_pdbs_config = yaml.safe_load(f)
```

    Reading PDB output configuration from data/output_pdbs_config.yaml


Make output directory:


```python
os.makedirs(config['pdb_outputs_dir'], exist_ok=True)
```

Read escape fractions and compute **total** and **maximum** escape at each site, and also the total and maximum escape at each site normalized to be between 0 and 1 for each selection:


```python
print(f"Reading escape fractions from {config['escape_fracs']}")

escape_fracs = (
    pd.read_csv(config['escape_fracs'])
    .query('library == "average"')
    .assign(site=lambda x: x['label_site'])
    .groupby(['selection', 'site'])
    .aggregate(total_escape=pd.NamedAgg(config['mut_metric'], 'sum'),
               max_escape=pd.NamedAgg(config['mut_metric'], 'max')
               )
    .reset_index()
    .assign(max_total_escape=lambda x: x.groupby('selection')['total_escape'].transform('max'),
            max_max_escape=lambda x: x.groupby('selection')['max_escape'].transform('max'),
            norm_total_escape=lambda x: x['total_escape'] / x['max_total_escape'],
            norm_max_escape=lambda x: x['max_escape'] / x['max_max_escape'])
    )

display(HTML(escape_fracs.head().to_html(index=False)))
```

    Reading escape fractions from results/escape_scores/escape_fracs.csv



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>selection</th>
      <th>site</th>
      <th>total_escape</th>
      <th>max_escape</th>
      <th>max_total_escape</th>
      <th>max_max_escape</th>
      <th>norm_total_escape</th>
      <th>norm_max_escape</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>CB6_400</td>
      <td>331</td>
      <td>0.040314</td>
      <td>0.013650</td>
      <td>17.83399</td>
      <td>0.9972</td>
      <td>0.002261</td>
      <td>0.013688</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>332</td>
      <td>0.085961</td>
      <td>0.026530</td>
      <td>17.83399</td>
      <td>0.9972</td>
      <td>0.004820</td>
      <td>0.026604</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>333</td>
      <td>0.048161</td>
      <td>0.013420</td>
      <td>17.83399</td>
      <td>0.9972</td>
      <td>0.002701</td>
      <td>0.013458</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>334</td>
      <td>0.043309</td>
      <td>0.006111</td>
      <td>17.83399</td>
      <td>0.9972</td>
      <td>0.002428</td>
      <td>0.006128</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>335</td>
      <td>0.032980</td>
      <td>0.004959</td>
      <td>17.83399</td>
      <td>0.9972</td>
      <td>0.001849</td>
      <td>0.004973</td>
    </tr>
  </tbody>
</table>


Now map the escape metrics to the B-factors.
For sites where no mutations have escape scores:
 - In the RBD chain(s) fill the B-factor for non-normalized scores to -1 to enable collapsing to zero or callout as a a separate class, depending how we choose to color sites for different visualizations. For normalized scores, fill to 0.
 - In other chains, always fill missing B factors to 0.  


```python
for name, specs in output_pdbs_config.items():
    print(f"\nMaking PDB mappings for {name} to {specs['pdbfile']}")
    assert os.path.isfile(specs['pdbfile'])
    
    # get escape fracs just for conditions of interest
    if isinstance(specs['conditions'], str) and specs['conditions'].upper() == 'ALL':
        conditions = escape_fracs['selection'].unique().tolist()
    else:
        assert isinstance(specs['conditions'], list)
        conditions = specs['conditions']
    print(f"Making mappings for {len(conditions)} conditions.")
    df = escape_fracs.query('selection in @conditions')
    
    # get chains
    assert isinstance(specs['chains'], list)
    print('Mapping to the following chains: ' + ', '.join(specs['chains']))
    df = pd.concat([df.assign(chain=chain) for chain in specs['chains']], ignore_index=True)
    
    # make mappings for each condition and metric
    for condition, df in df.groupby('selection'):
        print(f"  Writing B-factor re-assigned PDBs for {condition} to:")
    
        for metric in ['total_escape', 'max_escape', 'norm_total_escape', 'norm_max_escape']:
        
            # what do we assign to missing sites?
            missing_metric = collections.defaultdict(lambda: 0)  # non-RBD chains always fill to zero
            for chain in specs['chains']:
                if 'norm' in metric:
                    missing_metric[chain] = 0  # missing sites in RBD are 0 for normalized metric PDBs
                else:
                    missing_metric[chain] = -1  # missing sites in RBD are -1 for non-normalized metric PDBs
        
            fname = os.path.join(config['pdb_outputs_dir'], f"{condition}_{name}_{metric}.pdb")
            print(f"    {fname}")
            
            dms_variants.pdb_utils.reassign_b_factor(input_pdbfile=specs['pdbfile'],
                                                     output_pdbfile=fname,
                                                     df=df,
                                                     metric_col=metric,
                                                     missing_metric=missing_metric)
```

    
    Making PDB mappings for 6m0j to data/pdbs/6M0J.pdb
    Making mappings for 4 conditions.
    Mapping to the following chains: E
      Writing B-factor re-assigned PDBs for CB6_400 to:
        results/pdb_outputs/CB6_400_6m0j_total_escape.pdb
        results/pdb_outputs/CB6_400_6m0j_max_escape.pdb
        results/pdb_outputs/CB6_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/CB6_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for REGN10933+REGN10987_400 to:
        results/pdb_outputs/REGN10933+REGN10987_400_6m0j_total_escape.pdb
        results/pdb_outputs/REGN10933+REGN10987_400_6m0j_max_escape.pdb
        results/pdb_outputs/REGN10933+REGN10987_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/REGN10933+REGN10987_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for REGN10933_400 to:
        results/pdb_outputs/REGN10933_400_6m0j_total_escape.pdb
        results/pdb_outputs/REGN10933_400_6m0j_max_escape.pdb
        results/pdb_outputs/REGN10933_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/REGN10933_400_6m0j_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for REGN10987_400 to:
        results/pdb_outputs/REGN10987_400_6m0j_total_escape.pdb
        results/pdb_outputs/REGN10987_400_6m0j_max_escape.pdb
        results/pdb_outputs/REGN10987_400_6m0j_norm_total_escape.pdb
        results/pdb_outputs/REGN10987_400_6m0j_norm_max_escape.pdb
    
    Making PDB mappings for 6xdg to data/pdbs/6xdg.pdb
    Making mappings for 3 conditions.
    Mapping to the following chains: E
      Writing B-factor re-assigned PDBs for REGN10933+REGN10987_400 to:
        results/pdb_outputs/REGN10933+REGN10987_400_6xdg_total_escape.pdb
        results/pdb_outputs/REGN10933+REGN10987_400_6xdg_max_escape.pdb
        results/pdb_outputs/REGN10933+REGN10987_400_6xdg_norm_total_escape.pdb
        results/pdb_outputs/REGN10933+REGN10987_400_6xdg_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for REGN10933_400 to:
        results/pdb_outputs/REGN10933_400_6xdg_total_escape.pdb
        results/pdb_outputs/REGN10933_400_6xdg_max_escape.pdb
        results/pdb_outputs/REGN10933_400_6xdg_norm_total_escape.pdb
        results/pdb_outputs/REGN10933_400_6xdg_norm_max_escape.pdb
      Writing B-factor re-assigned PDBs for REGN10987_400 to:
        results/pdb_outputs/REGN10987_400_6xdg_total_escape.pdb
        results/pdb_outputs/REGN10987_400_6xdg_max_escape.pdb
        results/pdb_outputs/REGN10987_400_6xdg_norm_total_escape.pdb
        results/pdb_outputs/REGN10987_400_6xdg_norm_max_escape.pdb
    
    Making PDB mappings for 7c01 to data/pdbs/7c01_single.pdb
    Making mappings for 1 conditions.
    Mapping to the following chains: A
      Writing B-factor re-assigned PDBs for CB6_400 to:
        results/pdb_outputs/CB6_400_7c01_total_escape.pdb
        results/pdb_outputs/CB6_400_7c01_max_escape.pdb
        results/pdb_outputs/CB6_400_7c01_norm_total_escape.pdb
        results/pdb_outputs/CB6_400_7c01_norm_max_escape.pdb

