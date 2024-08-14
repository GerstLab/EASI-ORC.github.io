import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from ipywidgets import widgets
import os


# variables
save_btn = widgets.Button(description='Save Figure')
save_btn.style.button_color = 'GhostWhite'
graph_styles = ['default', '_mpl-gallery', '_mpl-gallery-nogrid', 'bmh', 'ggplot', 'seaborn-v0_8', 'seaborn-v0_8-bright', 'seaborn-v0_8-colorblind', 'seaborn-v0_8-dark', 'seaborn-v0_8-dark-palette', 'seaborn-v0_8-darkgrid', 'seaborn-v0_8-deep', 'seaborn-v0_8-muted', 'seaborn-v0_8-notebook', 'seaborn-v0_8-paper', 'seaborn-v0_8-pastel', 'seaborn-v0_8-poster', 'seaborn-v0_8-talk', 'seaborn-v0_8-ticks', 'seaborn-v0_8-white', 'seaborn-v0_8-whitegrid', 'tableau-colorblind10']
file_type=['png', 'jpg', 'tiff', 'pdf', 'ps', 'eps', 'svg']
axis_ticks=list(range(10,20))
DPI=list(range(500, 2501, 500)) 
main_title=list(range(8,33))
axis_labels=list(range(8,33))
legend=list(range(8,33))
plot_style=graph_styles
bins=(1, 200, 1)


# functions
def process_coords_column(coord_string):
    coords_list = []
    curr_coords = coord_string[1:-2]
    curr_coords = "]," + curr_coords + ",["
    curr_coords = curr_coords.split("],[")
    curr_coords = curr_coords[1:-1]
    for item in curr_coords:
        coords_list.append(item.split(","))
    return coords_list


def prepare_dataframes(df):
    df_file_convert_to_lists = df.dropna(axis=1, how='all').drop(columns=['Total mRNA per Cell', 'Total Colocolized With Organelle Near', 
                                                                        'Total Colocolized With Organelle Far', 'Total Not Colocolized with Organelle']).copy()
    df_file_convert_to_lists.loc[:, 'Spots Coordinates Intensity and Colocalization (Far Near or Not Colocolized)'] = \
                                (df_file_convert_to_lists['Spots Coordinates Intensity and Colocalization (Far Near or Not Colocolized)'].apply(process_coords_column))
    return df_file_convert_to_lists


def find_rna_z_plane(row, coverages_df, rna_df, idx):
    rna_z = str(rna_df.loc[row, 'z'])
    coverage_column = [col for col in coverages_df.columns if rna_z in col][0]
    rna_df.loc[row, 'organelle_coverage'] = coverages_df.loc[idx, coverage_column]
    return rna_df


def rna_status(row):
    if row['status'] == 'nc':
        row['not_colocalized'] = 1
    elif row['status'] == 'organelle_near':
        row['organelle_near_colocalized'] = 1
        row['colocalized'] = 1
    elif row['status'] == 'organelle_far':
        row['organelle_far_colocalized'] = 1
        row['colocalized'] = 1
    return row
   

def rna_df_function(idx, df):
    rna = df['Spots Coordinates Intensity and Colocalization (Far Near or Not Colocolized)'].iloc[idx]
    rna_df = pd.DataFrame(rna, columns=['x', 'y', 'z', 'intensity', 'status'])
    rna_df[['x', 'y', 'z', 'intensity']] = rna_df[['x', 'y', 'z', 'intensity']].apply(pd.to_numeric)
    rna_df['cell_id'] = df.index[idx]
    rna_df[['colocalized', 'not_colocalized', 'organelle_near_colocalized', 'organelle_far_colocalized']] = 0
    if len(rna_df) >0:
        for row in np.arange(len(rna_df)):
            rna_df = find_rna_z_plane(row, rna_df=rna_df, coverages_df=df, idx=idx)
            rna_df = rna_df.apply(lambda x: rna_status(row=x), axis=1).fillna(0, inplace=False)       
        return rna_df
    else:
        rna_df['organelle_coverage'] = None
        return rna_df
    

def rna_to_cell_df(df):
    copy_df = df.copy() 
    copy_df['full_id'] = copy_df['treatment'] + '_' + copy_df['cell_id'].astype(str) 
    rna_counts = copy_df['full_id'].value_counts()
    rna_counts = pd.DataFrame(rna_counts)
    rna_counts['full_id'] = rna_counts.index
    rna_counts.rename(columns={'count': 'rna_count'}, inplace=True)
    rna_counts = pd.DataFrame(rna_counts).reset_index(drop=True)
    copy_df = pd.merge(copy_df, rna_counts, on='full_id', how='left')
    numeric_cols = copy_df.select_dtypes(include='number').columns
    grouped_df = copy_df.groupby('full_id')[numeric_cols].mean().reset_index()
    grouped_df['cell_id'] = grouped_df['full_id'].str.split('_').str[1]
    grouped_df['treatment'] = grouped_df['full_id'].str.split('_').str[0]
    grouped_df.drop(columns=['x', 'y', 'z', 'intensity', 'full_id', 'organelle_coverage'], axis=1, inplace=True)
    return grouped_df


def plot_assist(axis_ticks, axis_labels, legend, main_title, x_axis_label, y_axis_label, main_title_label):
    matplotlib.rcParams["axes.spines.right"] = False
    matplotlib.rcParams["axes.spines.top"] = False
    plt.xticks(fontsize=axis_ticks)
    plt.yticks(fontsize=axis_ticks)
    plt.xlabel(x_axis_label, fontsize=axis_labels)
    plt.ylabel(y_axis_label, fontsize=axis_labels)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize =legend)
    plt.suptitle(main_title_label, fontsize=main_title, y = 1.05)


def filtered_mean_coverages(row, coverage_threshold, cols):
    row['filtered_mean_organelle_coverage'] = np.mean(list(filter(lambda cov: (cov > coverage_threshold[0]) and cov < coverage_threshold[1], np.array(row[cols]))))
    return row


def on_button_clicked(b):
    out = widgets.Output()
    os.makedirs(save_path, exist_ok=True)
    fig.savefig(plt_save_name, bbox_inches = 'tight', dpi = 1000)
    with out:
        out.clear_output()


def global_assist(min=None, max=None, DPI=None, fig_name=None, file_type=None, filter=True, intensity_threshold=None, coverage_threshold=None, rna_num_threshold=None, path=None):
    global fig, plt_save_name, dpi, save_btn, save_path
    fig = plt.gcf()
    dpi = DPI
    save_path = os.path.join(path, 'Figures')
    save_btn.on_click(on_button_clicked)
    if filter:
        global threshold
        threshold = [min, max]
        plt_save_name = rf"{save_path}{fig_name} ({round(min,2)} to {round(max,2)}).{str(file_type)}"
    else:
        plt_save_name = rf'{save_path}{fig_name}, intensity {intensity_threshold}, coverage {coverage_threshold}, rna_num {rna_num_threshold}.{str(file_type)}' 