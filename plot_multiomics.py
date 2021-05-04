import numpy as np
import matplotlib.pyplot as plt

n_modifications = 3


def plot_distribution_of_designs(df): 
    bar_height = 1
    labels = ['KO', 'NoMod', 'UP']
    colors = ['#019600', 'grey', '#219AD8']
        
    plt.style.use('seaborn-white')
    
    dataframe = df.copy()
    reactions = dataframe.columns
    
    n_rec = len(dataframe)
    dataframe.loc[n_rec] = [[list(dataframe[reaction]).count(int(i))/n_rec*100 
                             for i in range(n_modifications)]  for reaction in reactions]
    
    data = [ [dataframe.iloc[-1][r][num] for r in reactions] 
            for num in range(n_modifications)]
    
    y_pos = np.arange(len(reactions))

    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111)

    # Remove frame
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    patch_handles = []
    # left alignment of data starts at zero
    left = np.zeros(len(reactions)) 
    for i, d in enumerate(data):
        patch_handles.append(ax.barh(y_pos, d, 
                                     color=colors[i%len(colors)], edgecolor='white',
                                     height=bar_height, align='center', 
                                     left=left, label=labels[i]))
        left += d

    # search all of the bar segments and annotate
    for j in range(n_modifications):
        for i, patch in enumerate(patch_handles[j].get_children()):
            bl = patch.get_xy()
            x = 0.5*patch.get_width() + bl[0]
            y = 0.5*patch.get_height() + bl[1]
            ax.text(x,y, "%d%%" % (data[j][i]), ha='center')

    ax.set_title('Distribution of modifications')
    plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=True, 
                    labelbottom=False)
    plt.yticks(y_pos, reactions)
    ax.invert_yaxis()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    
    
def plot_DO_extmets(od,ext_metabolites):
    fig, ax = plt.subplots(figsize=(12,4), ncols=2, nrows=1)
    od.plot(ax=ax[0], style='s-', title='Cell', label='dcw', legend=True)
    ax[0].set_xlabel("Hour")
    ax[0].set_ylabel("Concentration [gDW/L]")
    ext_metabolites.plot(ax=ax[1], style='o-', title='External Metabolites')
    ax[1].set_xlabel("Hour")
    ax[1].set_ylabel("Concentration [mM]")
    
    
def pred_vs_actual(df):
    """Plots the predictions of a machine learning model.

    Create a bar plot of machine learning model predictions vs.
    actual values from the data set along with a 95% credible interval.
    """

    plt.style.use('seaborn-darkgrid')

    fontsize = 16

    predicted_mean = df['Mean predicted Isoprenol [mM]'][0]
    predicted_std = df['SD Isoprenol [mM]'][0]

    observed = df['Actual Isoprenol [mM]'][0]
    
    x_label = ['predicted', 'actual']
    x_pos = np.arange(len(x_label))
    width = 0.6  # the width of the bars

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.bar(0, predicted_mean,  yerr=1.96*predicted_std, width=width, 
           color='#397479', alpha=0.8, ecolor='#303030', capsize=4)
    ax.bar(1, observed, width, color='grey', alpha=0.8)

    plt.xticks(x_pos, x_label)

    ax.set_ylim([0, 0.7])

    ax.set_ylabel('Isoprenol [mM]', fontsize=fontsize)
    ax.set_title('Best recommendation', fontsize=fontsize)

    plt.tick_params(axis='both', which='major', labelsize=fontsize)

    # Save the figure and show
    plt.show()

    fig.savefig('../data/ART_prediction_vs_actual_recommendation.png',
                    bbox_inches='tight', transparent=False, dpi=300)
    
        