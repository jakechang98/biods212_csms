import pandas as pd
import matplotlib.pyplot as plt
import requests



def get_pathways(gene):
    user_list = {'list': '\n'.join(gene)}
    response = requests.post('https://maayanlab.cloud/Enrichr/addList', files=user_list)
    if response.status_code == 200:
        user_list_id = response.json()['userListId']
        enrichment_response = requests.get(f'https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType=KEGG_2021_Human')
        if enrichment_response.status_code == 200:
            return enrichment_response.json().get('KEGG_2021_Human', [])
    return []


def plot_enrichement_pathways(svgs):
  # Get the pathways for the genes
  genes = svgs['gene_name'].unique()
  all_enriched_pathways = []
  for gene in genes:
      pathways = get_pathways([gene])
      for pathway in pathways:
          all_enriched_pathways.append({
              'Gene': gene,
              'Pathway': pathway[1],
              'P-value': pathway[2],
              'Z-score': pathway[3]
          })

  # Create a DataFrame from the results
  plot_df = pd.DataFrame(all_enriched_pathways)
  
  # Plot the stacked bar plot
  # Pivot the data for plotting
  pivot_df = plot_df.pivot_table(index='Pathway', columns='Cell Type', values='P-value', fill_value=0)
  pivot_df.plot(kind='bar', stacked=True, figsize=(14, 8), colormap='tab10')
  plt.title('Relative Contribution of Pathways in Each Cell Type')
  plt.xlabel('Pathway')
  plt.ylabel('Z-score')
  plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
  plt.tight_layout()
  plt.tick_params(axis='x', rotation=90, size=15)
  plt.show()

