# Analyze viral escape-mutant selections
Analyze results from viral escape-mutant selections.

Import Python modules:


```python
import collections
import math
import os

import Bio.SeqIO

import dms_variants.constants
from dms_variants.constants import CBPALETTE
from dms_variants.utils import single_nt_accessible

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import mizani

import numpy

import pandas as pd

from plotnine import *

import yaml
```

Read in configuration and then escape-mutant selection results:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
    
print(f"Reading escape-mutant selection results from {config['escape_selection_results']}")
with open(config['escape_selection_results']) as f:
    selection_results = yaml.safe_load(f)
```

    Reading escape-mutant selection results from data/escape_selection_results.yaml


Make output directory:


```python
os.makedirs(config['escape_selections_dir'], exist_ok=True)
```

Read escape-mutation mapping and deep mutational scanning results, and then merge them:


```python
# read escape fractions
escape_fracs = (
    pd.read_csv(config['escape_fracs'])
    .query('library == "average"')
    .rename(columns={config['mut_metric']: 'mutation_escape',
                     config['site_metric']: 'site_escape',
                     'selection': 'antibody'})
    .assign(site=lambda x: x['label_site'])
    [['antibody', 'site', 'wildtype', 'mutation', 'mutation_escape', 'site_escape']]
    )

# read DMS data
bind_expr = (
    pd.read_csv(config['mut_bind_expr'])
    .rename(columns={'site_SARS2': 'site',
                     'bind_avg': 'ACE2 binding',
                     'expr_avg': 'RBD expression',
                     })
    .assign(mutation=lambda x: x['mutant'])
    [['site', 'mutation', 'ACE2 binding', 'RBD expression']]
    )

# merge escape and DMS data
escape_dms = (
    escape_fracs
    .merge(bind_expr,
           on=['site', 'mutation'],
           how='left',
           validate='many_to_one',
           )
    )

# first few lines of data frame
display(HTML(escape_dms.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>antibody</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>mutation_escape</th>
      <th>site_escape</th>
      <th>ACE2 binding</th>
      <th>RBD expression</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>CB6_400</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>0.001044</td>
      <td>0.03952</td>
      <td>-0.03</td>
      <td>-0.11</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>0.002087</td>
      <td>0.03952</td>
      <td>0.03</td>
      <td>-0.44</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>0.001044</td>
      <td>0.03952</td>
      <td>0.00</td>
      <td>-0.31</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>0.001044</td>
      <td>0.03952</td>
      <td>-0.10</td>
      <td>-0.70</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>331</td>
      <td>N</td>
      <td>G</td>
      <td>0.005013</td>
      <td>0.03952</td>
      <td>-0.04</td>
      <td>-0.25</td>
    </tr>
  </tbody>
</table>


Get data frame with just escape-selection counts and then add to the main data frame of escape / DMS data:


```python
records = []
antibody_order = {}
spike_nt_seq = {}
for selection_set, selection_set_d in selection_results.items():
    antibody_order[selection_set] = {}
    spike_nt_seq[selection_set] = str(Bio.SeqIO.read(selection_set_d['spike_sequence'],
                                                     'genbank',
                                                     ).seq)
    for antibody, antibody_d in selection_set_d['antibodies'].items():
        antibody_order[selection_set][antibody] = antibody_d['display_name']
        if 'mutations' in antibody_d:
            for mutation_str, n in antibody_d['mutations'].items():
                wt = mutation_str[0]
                site = int(mutation_str[1: -1])
                mutation = mutation_str[-1]
                if ('label_mutations' in antibody_d) and mutation_str in antibody_d['label_mutations']:
                    mutation_label = mutation_str
                else:
                    mutation_label = None
                records.append((selection_set, antibody, site, wt, mutation, n, mutation_label))
                assert 1 == len(escape_dms.query('antibody == @antibody')
                                          .query('wildtype == @wt')
                                          .query('site == @site')
                                          .query('mutation == @mutation')
                                          ), f"{mutation_str} not in `escape_dms` once for {antibody}"
        if 'label_mutations' in antibody_d:
            for mutation_str in antibody_d['label_mutations']:
                if 'mutations' not in antibody_d or mutation_str not in antibody_d['mutations']:
                    wt = mutation_str[0]
                    site = int(mutation_str[1: -1])
                    mutation = mutation_str[-1]
                    records.append((selection_set, antibody, site, wt, mutation, 0, mutation_str))
                    assert 1 == len(escape_dms.query('antibody == @antibody')
                                              .query('wildtype == @wt')
                                              .query('site == @site')
                                              .query('mutation == @mutation')
                                              ), f"{mutation_str} not in `escape_dms` once for {antibody}"
            
selection_df = pd.DataFrame.from_records(
                records,
                columns=['selection_set', 'antibody', 'site',
                         'wildtype', 'mutation', 'n_selected', 'mutation_label'],
                )

display(HTML(selection_df.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>selection_set</th>
      <th>antibody</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>n_selected</th>
      <th>mutation_label</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>vir_whelan_VSV_S2K146</td>
      <td>S2K146_63</td>
      <td>489</td>
      <td>Y</td>
      <td>H</td>
      <td>36</td>
      <td>Y489H</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2K146_63</td>
      <td>489</td>
      <td>Y</td>
      <td>H</td>
      <td>36</td>
      <td>Y489H</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2E12_56</td>
      <td>476</td>
      <td>G</td>
      <td>S</td>
      <td>1</td>
      <td>G476S</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2E12_56</td>
      <td>475</td>
      <td>A</td>
      <td>D</td>
      <td>1</td>
      <td>A475D</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2E12_56</td>
      <td>476</td>
      <td>G</td>
      <td>D</td>
      <td>8</td>
      <td>G476D</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2E12_56</td>
      <td>485</td>
      <td>G</td>
      <td>D</td>
      <td>6</td>
      <td>G485D</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2D106_68</td>
      <td>484</td>
      <td>E</td>
      <td>A</td>
      <td>1</td>
      <td>E484A</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2D106_68</td>
      <td>490</td>
      <td>F</td>
      <td>C</td>
      <td>2</td>
      <td>F490C</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2D106_68</td>
      <td>484</td>
      <td>E</td>
      <td>K</td>
      <td>11</td>
      <td>E484K</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2D106_68</td>
      <td>484</td>
      <td>E</td>
      <td>V</td>
      <td>1</td>
      <td>E484V</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>346</td>
      <td>R</td>
      <td>G</td>
      <td>2</td>
      <td>R346G</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>484</td>
      <td>E</td>
      <td>G</td>
      <td>1</td>
      <td>E484G</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>490</td>
      <td>F</td>
      <td>S</td>
      <td>2</td>
      <td>F490S</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>449</td>
      <td>Y</td>
      <td>C</td>
      <td>1</td>
      <td>Y449C</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>448</td>
      <td>N</td>
      <td>D</td>
      <td>1</td>
      <td>N448D</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>450</td>
      <td>N</td>
      <td>D</td>
      <td>1</td>
      <td>N450D</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>450</td>
      <td>N</td>
      <td>K</td>
      <td>2</td>
      <td>N450K</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>452</td>
      <td>L</td>
      <td>R</td>
      <td>2</td>
      <td>L452R</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>490</td>
      <td>F</td>
      <td>V</td>
      <td>1</td>
      <td>F490V</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>346</td>
      <td>R</td>
      <td>I</td>
      <td>1</td>
      <td>R346I</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>490</td>
      <td>F</td>
      <td>L</td>
      <td>6</td>
      <td>F490L</td>
    </tr>
    <tr>
      <td>vir_whelan_VSV_main</td>
      <td>S2X58_18</td>
      <td>450</td>
      <td>N</td>
      <td>H</td>
      <td>3</td>
      <td>N450H</td>
    </tr>
  </tbody>
</table>


Add escape-selection counts and labels to main data frame with escape scores, and just keep antibodies of interest for each selection:


```python
escape_dms_selection = (
    pd.concat([escape_dms.merge(df,
                                on=['antibody', 'site', 'wildtype', 'mutation'],
                                how='left',
                                validate='one_to_one',
                                )
                          .assign(selection_set=selection_set,
                                  antibody_name=lambda x: x['antibody'].map(antibody_order[selection_set])
                                  )
                          .query('antibody_name.notnull()')
               for selection_set, df in selection_df.groupby('selection_set')
               ])
    .assign(n_selected=lambda x: x['n_selected'].fillna(0).astype(int),
            n_selected_total=lambda x: x.groupby(['selection_set', 'antibody'])
                                        ['n_selected'].transform('sum'),
            any_selected=lambda x: x['n_selected'] > 0,
            max_at_site_escape=lambda x: x.groupby(['antibody', 'site'])['mutation_escape'].transform('max'),
            mean_at_site_escape=lambda x: x.groupby(['antibody', 'site'])['mutation_escape'].transform('mean'),
            )
    )

display(HTML(escape_dms_selection.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>antibody</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>mutation_escape</th>
      <th>site_escape</th>
      <th>ACE2 binding</th>
      <th>RBD expression</th>
      <th>selection_set</th>
      <th>n_selected</th>
      <th>mutation_label</th>
      <th>antibody_name</th>
      <th>n_selected_total</th>
      <th>any_selected</th>
      <th>max_at_site_escape</th>
      <th>mean_at_site_escape</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>S2K146_63</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>0.004645</td>
      <td>0.026</td>
      <td>-0.03</td>
      <td>-0.11</td>
      <td>vir_whelan_VSV_S2K146</td>
      <td>0</td>
      <td>NaN</td>
      <td>S2K146</td>
      <td>36</td>
      <td>False</td>
      <td>0.005989</td>
      <td>0.001625</td>
    </tr>
    <tr>
      <td>S2K146_63</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>0.005989</td>
      <td>0.026</td>
      <td>0.03</td>
      <td>-0.44</td>
      <td>vir_whelan_VSV_S2K146</td>
      <td>0</td>
      <td>NaN</td>
      <td>S2K146</td>
      <td>36</td>
      <td>False</td>
      <td>0.005989</td>
      <td>0.001625</td>
    </tr>
    <tr>
      <td>S2K146_63</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>0.000793</td>
      <td>0.026</td>
      <td>0.00</td>
      <td>-0.31</td>
      <td>vir_whelan_VSV_S2K146</td>
      <td>0</td>
      <td>NaN</td>
      <td>S2K146</td>
      <td>36</td>
      <td>False</td>
      <td>0.005989</td>
      <td>0.001625</td>
    </tr>
    <tr>
      <td>S2K146_63</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>0.003720</td>
      <td>0.026</td>
      <td>-0.10</td>
      <td>-0.70</td>
      <td>vir_whelan_VSV_S2K146</td>
      <td>0</td>
      <td>NaN</td>
      <td>S2K146</td>
      <td>36</td>
      <td>False</td>
      <td>0.005989</td>
      <td>0.001625</td>
    </tr>
    <tr>
      <td>S2K146_63</td>
      <td>331</td>
      <td>N</td>
      <td>G</td>
      <td>0.002445</td>
      <td>0.026</td>
      <td>-0.04</td>
      <td>-0.25</td>
      <td>vir_whelan_VSV_S2K146</td>
      <td>0</td>
      <td>NaN</td>
      <td>S2K146</td>
      <td>36</td>
      <td>False</td>
      <td>0.005989</td>
      <td>0.001625</td>
    </tr>
  </tbody>
</table>


Now add in the viral codon sequence and determine what amino-acid mutations are single-nucleotide accessible:


```python
codon_df = pd.DataFrame()
sites = escape_dms_selection['site'].unique()
for selection_set, spike in spike_nt_seq.items():
    codon_df = codon_df.append(
        pd.DataFrame({'selection_set': selection_set, 'site': sites})
          .assign(viral_codon=lambda x: x['site'].map(lambda r: spike[3 * (r - 1): 3 * r]),
                  viral_aa=lambda x: x['viral_codon'].map(dms_variants.constants.CODON_TO_AA))
        )
    
escape_dms_selection = (
    escape_dms_selection
    .drop(columns=['viral_codon', 'viral_aa'], errors='ignore')
    .merge(codon_df,
           how='left',
           on=['selection_set', 'site'],
           validate='many_to_one',
           )
    .assign(single_nt_accessible=lambda x: x.apply(lambda row: single_nt_accessible(row['viral_codon'],
                                                                                    row['mutation']),
                                                   axis=1)
            )
    )

# check viral spike is same one used in our DMS. If not, we need to somehow deal with that...
aa_mismatch = escape_dms_selection.query('wildtype != viral_aa')
if len(aa_mismatch):
    raise ValueError('mismatches at the following amino acids: ' + aa_mismatch)

display(HTML(escape_dms_selection.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>antibody</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>mutation_escape</th>
      <th>site_escape</th>
      <th>ACE2 binding</th>
      <th>RBD expression</th>
      <th>selection_set</th>
      <th>n_selected</th>
      <th>mutation_label</th>
      <th>antibody_name</th>
      <th>n_selected_total</th>
      <th>any_selected</th>
      <th>max_at_site_escape</th>
      <th>mean_at_site_escape</th>
      <th>viral_codon</th>
      <th>viral_aa</th>
      <th>single_nt_accessible</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>S2K146_63</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>0.004645</td>
      <td>0.026</td>
      <td>-0.03</td>
      <td>-0.11</td>
      <td>vir_whelan_VSV_S2K146</td>
      <td>0</td>
      <td>NaN</td>
      <td>S2K146</td>
      <td>36</td>
      <td>False</td>
      <td>0.005989</td>
      <td>0.001625</td>
      <td>AAT</td>
      <td>N</td>
      <td>False</td>
    </tr>
    <tr>
      <td>S2K146_63</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>0.005989</td>
      <td>0.026</td>
      <td>0.03</td>
      <td>-0.44</td>
      <td>vir_whelan_VSV_S2K146</td>
      <td>0</td>
      <td>NaN</td>
      <td>S2K146</td>
      <td>36</td>
      <td>False</td>
      <td>0.005989</td>
      <td>0.001625</td>
      <td>AAT</td>
      <td>N</td>
      <td>True</td>
    </tr>
    <tr>
      <td>S2K146_63</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>0.000793</td>
      <td>0.026</td>
      <td>0.00</td>
      <td>-0.31</td>
      <td>vir_whelan_VSV_S2K146</td>
      <td>0</td>
      <td>NaN</td>
      <td>S2K146</td>
      <td>36</td>
      <td>False</td>
      <td>0.005989</td>
      <td>0.001625</td>
      <td>AAT</td>
      <td>N</td>
      <td>False</td>
    </tr>
    <tr>
      <td>S2K146_63</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>0.003720</td>
      <td>0.026</td>
      <td>-0.10</td>
      <td>-0.70</td>
      <td>vir_whelan_VSV_S2K146</td>
      <td>0</td>
      <td>NaN</td>
      <td>S2K146</td>
      <td>36</td>
      <td>False</td>
      <td>0.005989</td>
      <td>0.001625</td>
      <td>AAT</td>
      <td>N</td>
      <td>False</td>
    </tr>
    <tr>
      <td>S2K146_63</td>
      <td>331</td>
      <td>N</td>
      <td>G</td>
      <td>0.002445</td>
      <td>0.026</td>
      <td>-0.04</td>
      <td>-0.25</td>
      <td>vir_whelan_VSV_S2K146</td>
      <td>0</td>
      <td>NaN</td>
      <td>S2K146</td>
      <td>36</td>
      <td>False</td>
      <td>0.005989</td>
      <td>0.001625</td>
      <td>AAT</td>
      <td>N</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


Make plots showing effects of all mutations stratified by how many times observed in selections, and whether single-nucleotide accessible or not from the viral spike.
For mutations observed in selections point size, is proportional to times observed:


```python
for selection_set, df in escape_dms_selection.groupby('selection_set'):
    
    min_size = selection_results[selection_set]['min_size']
    max_size = selection_results[selection_set]['max_size']
    size_scale = selection_results[selection_set]['max_size']
                                                  
    if 'escape_metric' in selection_results[selection_set]:
        escape_metric = selection_results[selection_set]['escape_metric']
    else:
        escape_metric = 'mutation_escape'
    assert escape_metric in df.columns
                                                  
    if 'ncol' in selection_results[selection_set]:
        ncol = min(df['antibody'].nunique(), selection_results[selection_set]['ncol'])
    else:
        ncol = df['antibody'].nunique()
    nrow = math.ceil(df['antibody'].nunique() / ncol)
    
    print(f"\nMaking plot for {selection_set}")
    
    selected_not_accessible = len(df.query('any_selected and not single_nt_accessible'))
    
    if 'custom_categories' in selection_results[selection_set]:
        custom_cats = selection_results[selection_set]['custom_categories']
        addtl_cats = list({cat: True for a_cats in custom_cats.values() for cat in a_cats.values()})
    else:
        custom_cats = {}
        addtl_cats = []
    
    def get_point_category(row):
        if row['antibody_name'] in custom_cats and (row['mutation_str'] in 
                                                    custom_cats[row['antibody_name']]):
            return custom_cats[row['antibody_name']][row['mutation_str']]
        elif row['any_selected'] and row['single_nt_accessible']:
            return 'selected (single-nucleotide)'
        elif row['any_selected']:
            return 'selected (multi-nucleotide)'
        elif row['single_nt_accessible']:
            return 'single-nucleotide'
        else:
            return 'multi-nucleotide'
        
    df = df.assign(
            mutation_str=lambda x: x['wildtype'] + x['site'].astype(str) + x['mutation'],
            antibody_name=lambda x: pd.Categorical(x['antibody_name'],
                                                   antibody_order[selection_set].values(),
                                                   ordered=True),
            point_area=lambda x: 0.5 * size_scale * numpy.clip(x['n_selected'].fillna(0),
                                                               min_size,
                                                               max_size),
            point_category=lambda x: x.apply(get_point_category, axis=1),
            )
    possible_point_categories = ['single-nucleotide', 'multi-nucleotide',
                                 'selected (single-nucleotide)',
                                 'selected (multi-nucleotide)'] + addtl_cats
    observed_point_categories = [cat for cat in possible_point_categories
                                 if cat in df['point_category'].unique()]
    df = df.assign(point_category=lambda x: pd.Categorical(x['point_category'],
                                                           observed_point_categories,
                                                           ordered=True)
                   )
    
    if 'shapes' in selection_results[selection_set]:
        cat_shapes = selection_results[selection_set]['shapes']
    else:
        cat_shapes = ['o', 'x', 'D', '^']
    assert len(cat_shapes) >= len(observed_point_categories), 'not enough shapes'
    if 'colors' in selection_results[selection_set]:
        cat_colors = selection_results[selection_set]['colors']
    else: 
        cat_colors = [CBPALETTE[2], CBPALETTE[0], CBPALETTE[1], CBPALETTE[1]]
    assert len(cat_colors) >= len(observed_point_categories), 'not enough colors'
    if 'alphas' in selection_results[selection_set]:
        cat_alphas = selection_results[selection_set]['alphas']
    else:  
        cat_alphas = [0.4, 0.3, 0.95, 0.95]
    assert len(cat_colors) >= len(observed_point_categories), 'not enough alphas'
    
    p = (ggplot(df) +
         aes(escape_metric, 'ACE2 binding',
             color='point_category', alpha='point_category',
             size='point_area', shape='point_category', label='mutation_label') +
         geom_point() +
         geom_text(data=df.query('mutation_label.notnull()'),
                   size=(selection_results[selection_set]['label_fontsize'] if
                         'label_fontsize' in selection_results[selection_set]
                         else 8),
                   va='top', ha='right', alpha=1, nudge_x=-0.025, nudge_y=-0.025,
                   show_legend=False,
                   # see here for adjust_text: https://stackoverflow.com/a/57717833
                   adjust_text={'avoid_text': True,
                                'avoid_points': False,
                                'avoid_self': True,
                                'expand_text': [1.05, 1.2]},
                   ) +
         facet_wrap('~ antibody_name',
                    scales=(selection_results[selection_set]['facet_scales']
                            if 'facet_scales' in selection_results[selection_set]
                            else 'fixed'),
                    ncol=ncol) +
         scale_color_manual(values=cat_colors) +
         scale_alpha_manual(values=cat_alphas) +
         scale_size_area(limits=(size_scale * min_size, size_scale * max_size),
                                 max_size=4) +
         scale_shape_manual(values=cat_shapes) +
         theme_classic() +
         scale_x_continuous(name={'mutation_escape': 'mutation escape fraction',
                                  'site_escape': 'total escape at site',
                                  'max_at_site_escape': 'max escape at site',
                                  'mean_at_site_escape': 'mean escape at site',
                                  }[escape_metric],
                            breaks=mizani.breaks.mpl_breaks(nbins=3),
                            
                            ) +
         scale_y_continuous(expand=(0.07, 0)) +
         guides(alpha=False, size=False,
                shape=guide_legend(override_aes={'size': 3},
                                   title='mutation type'),
                color=guide_legend(title='mutation type'),
                ) +
         theme(figure_size=(2.3 * ncol, 2.3 * nrow),
               legend_position='top' if 'legend_position' not in selection_results[selection_set]
                                else selection_results[selection_set]['legend_position'])
         )
    
    # ad hoc code to put desired points on top: those with only a few of that color
    fig = p.draw()
    for child in fig.get_children():
        if 'AxesSubplot' in str(child):
            for c in child.collections:
                if len(c._linewidths) < 20:
                    c.zorder = 2
    
    plotfile = os.path.join(config['escape_selections_dir'], f"{selection_set}.pdf")
    svgfile = os.path.splitext(plotfile)[0] + '.svg'
    print(f"Saving plot to {plotfile} and {svgfile}")
    fig.savefig(plotfile, bbox_inches='tight')
    fig.savefig(svgfile, bbox_inches='tight')
    display(fig)
    plt.close(fig)
```

    
    Making plot for vir_whelan_VSV_S2K146
    Saving plot to results/escape_selections/vir_whelan_VSV_S2K146.pdf and results/escape_selections/vir_whelan_VSV_S2K146.svg



    
![png](escape_selections_files/escape_selections_15_1.png)
    


    
    Making plot for vir_whelan_VSV_main


    /fh/fast/bloom_j/software/miniconda3/envs/SARS-CoV-2-RBD_MAP/lib/python3.7/site-packages/plotnine/guides/guide_legend.py:126: PlotnineWarning: Duplicated override_aes is ignored.


    Saving plot to results/escape_selections/vir_whelan_VSV_main.pdf and results/escape_selections/vir_whelan_VSV_main.svg



    
![png](escape_selections_files/escape_selections_15_5.png)
    



```python

```
