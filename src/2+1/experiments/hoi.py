import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os


k0_values = [1.0, 5.0, 7.0]
T_values = [3, 5, 10]
target_volume = 3000
ksteps = target_volume * 100

fig, axes = plt.subplots(len(k0_values), len(T_values), figsize=(10, 10))

for i, k0 in enumerate(k0_values):
    for j, T in enumerate(T_values):
        path = f'measurements_thermal/k0={k0}'
        filepaths_ar = glob.glob(f'{path}/T{T}*tswps=1000_*kstps={ksteps}*acceptance_ratios.npy')
        ar_add = []
        ar_delete = []
        ar_flip = []
        ar_shift = []
        ar_ishift = []

        for filepath in filepaths_ar:
            ar = np.load(filepath)
            ar_add.append(ar[:, 0])
            ar_delete.append(ar[:, 1])
            ar_flip.append(ar[:, 2])
            ar_shift.append(ar[:, 3])
            ar_ishift.append(ar[:, 4])

        df_ar_add = pd.DataFrame(ar_add)
        df_ar_delete = pd.DataFrame(ar_delete)
        df_ar_flip = pd.DataFrame(ar_flip)
        df_ar_shift = pd.DataFrame(ar_shift)
        df_ar_ishift = pd.DataFrame(ar_ishift)

        df_ar_add = df_ar_add.melt(var_name='sweep', value_name='ar', ignore_index=False).reset_index()
        df_ar_delete = df_ar_delete.melt(var_name='sweep', value_name='ar', ignore_index=False).reset_index()
        df_ar_flip = df_ar_flip.melt(var_name='sweep', value_name='ar', ignore_index=False).reset_index()
        df_ar_shift = df_ar_shift.melt(var_name='sweep', value_name='ar', ignore_index=False).reset_index()
        df_ar_ishift = df_ar_ishift.melt(var_name='sweep', value_name='ar', ignore_index=False).reset_index()

        df_ar_add['sweep'] += 1
        df_ar_delete['sweep'] += 1
        df_ar_flip['sweep'] += 1
        df_ar_shift['sweep'] += 1
        df_ar_ishift['sweep'] += 1

        ax = axes[i, j]
        sns.lineplot(data=df_ar_add, x='sweep', y='ar', errorbar='sd', label='Add', color='royalblue', ax=ax)
        sns.lineplot(data=df_ar_delete, x='sweep', y='ar', errorbar='sd', label='Delete', color='orangered', ax=ax)
        sns.lineplot(data=df_ar_flip, x='sweep', y='ar', errorbar='sd', label='Flip', color='seagreen', ax=ax)
        sns.lineplot(data=df_ar_shift, x='sweep', y='ar', errorbar='sd', label='Shift', color='red', ax=ax)
        sns.lineplot(data=df_ar_ishift, x='sweep', y='ar', errorbar='sd', label='Inverse Shift', color='purple', ax=ax)
        ax.set_title(f'$T$={T}, $k_0$={k0}', fontsize=14)

        # Remove x and y labels
        ax.set_xlabel('')
        ax.set_ylabel('')

        # Only keep axis ticks on the left and bottom
        if i != len(k0_values) - 1:
            ax.set_xticks([])
        if j != 0:
            ax.set_yticks([])

        ax.get_legend().remove()

# Add general labels to the figure
fig.suptitle('Acceptance probabilities for different $k_0$ and $T$ values', fontsize=16)
fig.text(0.51, 0.2, 'Sweep', ha='center', fontsize=16)
fig.text(0, 0.59, 'Acceptance probability', va='center', rotation='vertical', fontsize=16)

# Remove rendundant legend labels
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='lower right', bbox_to_anchor=(1, -1.1), fontsize=14)

# Put the legend out of the figure
plt.subplots_adjust(right=0.85)
plt.tight_layout()
plt.savefig('plots/acceptance_ratios/acceptance_probabilities_T.png', dpi=400, bbox_inches='tight')
plt.show()