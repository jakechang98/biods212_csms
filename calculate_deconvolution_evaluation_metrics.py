import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr, spearmanr
import os
import scanpy as sc
import squidpy as sq
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns


dataset_names = ['SN048_A121573_Rep1', 'SN048_A121573_Rep2', 'SN048_A416371_Rep1', 'SN048_A416371_Rep2','SN123_A551763_Rep1', 'SN84_A120838_Rep1', 'SN84_A120838_Rep2']

def read_csv_to_df(file_path):
    return pd.read_csv(file_path)


def calculate_deconvolution_validation_metrics(cside_output_directory, cytospace_output_directory):
  # Initialize dictionaries to store DataFrames
  cside_dfs = {}
  cytospace_dfs = {}
  location_dfs = {}

  total_rmse_values = {}
  total_pearson_correlations = {}
  total_spearman_correlations = {}
  adatas = {}

  for name in dataset_names:
      results_file = os.path.join(cside_output_directory, f"{name}_results.csv")
      if os.path.exists(results_file):
          cside_dfs[name] = read_csv_to_df(results_file)
          cside_dfs[name].rename(columns={'Unnamed: 0':'SpotID'}, inplace=True)
          cside_dfs[name]['SpotID'] = cside_dfs[name]['SpotID'].str.replace('.', '-', regex=False)
          print(f"Read results DataFrame for {name}")
          print(len(cside_dfs[name].index))

      cell_type_file = os.path.join(cytospace_output_directory, f"{name}/fractional_abundances_by_spot.csv")
      locations_file = os.path.join(cytospace_output_directory, f"{name}/assigned_locations.csv")
      if os.path.exists(cell_type_file):
          cytospace_dfs[name] = read_csv_to_df(cell_type_file)
          location_dfs[name] = read_csv_to_df(locations_file)
          print(f"Read cell types DataFrame for {name}")
          print(len(cytospace_dfs[name].index))
      
      print("Before Intersection")
      print(len(cytospace_dfs[name].index))
      print(len(cside_dfs[name].index))

      cytospace_dfs[name] = cytospace_dfs[name][cytospace_dfs[name]['SpotID'].isin(cside_dfs[name]['SpotID'])]
      cside_dfs[name] = cside_dfs[name][cside_dfs[name]['SpotID'].isin(cytospace_dfs[name]['SpotID'])]
      location_dfs[name] = location_dfs[name][location_dfs[name]['SpotID'].isin(cytospace_dfs[name]['SpotID'])]
      location_dfs[name] = location_dfs[name][['SpotID', 'row', 'col']].copy()
      location_dfs[name] = location_dfs[name].drop_duplicates()
      
      cside_dfs[name] = cside_dfs[name].sort_values(by='SpotID')
      cytospace_dfs[name] = cytospace_dfs[name].sort_values(by='SpotID')
      location_dfs[name] = location_dfs[name].sort_values(by='SpotID')
      
      print("After Intersection")
      print(len(cytospace_dfs[name].index))
      print(len(cside_dfs[name].index))

      cside_dfs[name].rename(columns={
      'B cells_1': 'B cells 1', 'B cells_2': 'B cells 2', 'B cells_3': 'B cells 3',
      'Epithelial cells_1': 'Epithelial cells 1', 'Epithelial cells_2': 'Epithelial cells 2',
      'Myeloids_1': 'Myeloids 1', 'Myeloids_2': 'Myeloids 2',
      'Stromal cells_1': 'Stromal cells 1', 'Stromal cells_2': 'Stromal cells 2', 'Stromal cells_3': 'Stromal cells 3',
      'T cells_1': 'T cells 1', 'T cells_2': 'T cells 2'}, inplace=True)

      cell_types = list(cside_dfs[name].columns)[1:]

      rmse_values = {}
      pearson_correlations = {}
      spearman_correlations = {}

      # make adatas
      cside_annData = ad.AnnData(X=cside_dfs[name].iloc[:, 1:].values, 
                                obs=location_dfs[name][['row', 'col']], 
                                var=pd.DataFrame(index=cside_dfs[name].columns[1:]))

      cytospace_annData = ad.AnnData(X=cytospace_dfs[name].iloc[:, 1:].values, 
                                    obs=location_dfs[name][['row', 'col']], 
                                    var=pd.DataFrame(index=cytospace_dfs[name].columns[1:]))

      cside_annData.obsm['spatial'] = location_dfs[name][['row', 'col']].values
      cytospace_annData.obsm['spatial'] = location_dfs[name][['row', 'col']].values
      adatas[name] = {'cside': cside_annData, 'cytospace': cytospace_annData}

      # calcualte morans index
      sq.gr.spatial_neighbors(adatas[name]['cside'], coord_type='grid')
      sq.gr.spatial_autocorr(adatas[name]['cside'], mode="moran")
      sq.gr.spatial_neighbors(adatas[name]['cytospace'], coord_type='grid')
      sq.gr.spatial_autocorr(adatas[name]['cytospace'], mode="moran")

      for cell_type in cell_types:
          if cell_type in cside_dfs[name].columns and cell_type in cytospace_dfs[name].columns:
              # RMSE
              rmse = np.sqrt(mean_squared_error(cside_dfs[name][cell_type].values, cytospace_dfs[name][cell_type]))
              rmse_values[cell_type] = rmse
              
              # Pearson
              corr, _ = pearsonr(cside_dfs[name][cell_type], cytospace_dfs[name][cell_type])
              pearson_correlations[cell_type] = corr

              # Spearman
              corr, _ = spearmanr(cside_dfs[name][cell_type], cytospace_dfs[name][cell_type])
              spearman_correlations[cell_type] = corr

                  
          total_rmse_values[name] = rmse_values
          print(f"Calcualted RMSE value for {name}")

          total_pearson_correlations[name] = pearson_correlations
          print(f"Calcualted Pearson Correlation value for {name}")
          
          total_spearman_correlations[name] = spearman_correlations
          print(f"Calcualted Spearman Correlation value for {name}")

  rmse_df = pd.DataFrame.from_dict(total_rmse_values, orient='index')
  pearson_df = pd.DataFrame.from_dict(total_pearson_correlations, orient='index')
  spearman_df = pd.DataFrame.from_dict(total_spearman_correlations, orient='index')
  rmse_df.to_csv('./cytpsoace_cside_rmse.csv')
  pearson_df.to_csv('./cytpsoace_cside_pearson.csv')
  spearman_df.to_csv('./cytpsoace_cside_spearman.csv')

def make_moransI_plots(adatas):
  morans_cside = {}
  morans_cytospace = {}

  for name in dataset_names:
    morans_i_results = adatas[name]['cytospace'].uns['moranI']
    morans_i_df = pd.DataFrame(morans_i_results)
    morans_i_df['Dataset'] = name  # Add a column to identify the dataset
    morans_cytospace[name] = morans_i_df.sort_values('I')

  combined_morans_cytospace_df = pd.concat(morans_cytospace.values())
  combined_morans_cytospace_df['Gene'] = combined_morans_cytospace_df.index
    
  for name in dataset_names:
      morans_i_results = adatas[name]['cside'].uns['moranI']
      morans_i_df = pd.DataFrame(morans_i_results)
      morans_i_df['Dataset'] = name  # Add a column to identify the dataset
      morans_cside[name] = morans_i_df.sort_values('I')

  combined_morans_cside_df = pd.concat(morans_cside.values())
  combined_morans_cside_df['Gene'] = combined_morans_cside_df.index

  all_genes = sorted(set(combined_morans_cytospace_df['Gene']).union(set(combined_morans_cside_df['Gene'])))
  combined_morans_cytospace_df['Gene'] = pd.Categorical(combined_morans_cytospace_df['Gene'], categories=all_genes, ordered=True)
  combined_morans_cside_df['Gene'] = pd.Categorical(combined_morans_cside_df['Gene'], categories=all_genes, ordered=True)

  # Plot the combined Moran's I results side by side
  fig, axes = plt.subplots(1, 2, figsize=(24, 8), sharey=True)

  # Cytospace Plot
  sns.scatterplot(data=combined_morans_cytospace_df, x='Gene', y='I', hue='Dataset', palette='tab10', ax=axes[0], s=100)
  axes[0].set_xlabel('Cell Type', size=15)
  axes[0].set_ylabel('Moran\'s I', size=15)
  axes[0].set_title('CytoSPACE Moran\'s I', size=20)
  axes[0].set_xticklabels(combined_morans_cytospace_df['Gene'], size=15)
  axes[0].tick_params(axis='x', rotation=45, size=30)
  axes[0].legend(title='Sample')

  # Cside Plot
  sns.scatterplot(data=combined_morans_cside_df, x='Gene', y='I', hue='Dataset', palette='tab10', ax=axes[1], s=100)
  axes[1].set_xlabel('Cell Type', size=15)
  axes[1].set_ylabel('Moran\'s I', size=15)
  axes[1].set_title('C-SIDE Moran\'s I', size=20)
  axes[1].set_xticklabels(combined_morans_cytospace_df['Gene'], size=15)
  axes[1].tick_params(axis='x', rotation=45, size=30)
  axes[1].legend(title='Sample')

  plt.tight_layout()
  plt.show()
        

