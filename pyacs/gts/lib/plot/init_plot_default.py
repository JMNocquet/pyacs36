

###############################################################################
# GLOBAL PARAMETERS FOR PLOTS         
###############################################################################

# FIGURE SIZE

plot_size_1_component=(10,3.5)
plot_size_2_component=(8,8)
plot_size_3_component=(7,8)

# SPACE BETWEEN SUBPLOTS
hspace=0.3

# PLOT SETTINGS FOR MAIN MAKERS,ERRORBARS AND SUPERIMPOSED  
marker_main_symbol='o'
marker_main_color='b'
marker_main_size=4
marker_main_colorlw=0.0
marker_main_linewidth=0
marker_main_linecolor='k'

error_bar_color='grey'
error_bar_linewidth=1
error_bar_capsize=0

line_2_color='r'
line_2_width=.5
line_2_style='-'
marker_2_symbol=None
marker_2_color='r'
marker_2_size=None
marker_2_colorlw=0.0
marker_2_linewidth=1
marker_2_linecolor='r'

H_kwargs_plot_date={}
H_kwargs_plot_date_superimposed={}
H_kwargs_errorbar={}
H_kwargs_errorbar_superimposed={}

###############################################################################
# INIT PLOT SETTINGS         
###############################################################################

def __init_kwargs():

    global H_kwargs_errorbar
    global H_kwargs_errorbar_superimposed
    
    H_plot_settings = {'plot_size_1_component': plot_size_1_component, \
                       'plot_size_2_component': plot_size_2_component, \
                       'plot_size_3_component': plot_size_3_component, \
                       'hspace' : hspace, \
                        }
    
    H_kwargs_errorbar={'linewidth' : 0,\
                       'marker' : marker_main_symbol ,\
                       'markersize' : marker_main_size, \
                       'markerfacecolor' : marker_main_color,\
                       'markeredgecolor' : marker_main_color,\
                       'markeredgewidth' : marker_main_colorlw,\
                       'ecolor' : error_bar_color,\
                       'elinewidth' : error_bar_linewidth,\
                       'capsize' : error_bar_capsize }
    
    H_kwargs_errorbar_superimposed={ \
                        'color' : line_2_color,\
                        'linewidth' : line_2_width,\
                        'linestyle' : line_2_style,\
                        'marker' : marker_2_symbol,\
                        'markersize' : marker_2_size, \
                        'markerfacecolor' : marker_2_color,\
                        'markeredgecolor' : marker_2_color,\
                        'markeredgewidth' : marker_2_colorlw, \
                        'ecolor' : error_bar_color,\
                        'elinewidth' : error_bar_linewidth,\
                        'capsize' : error_bar_capsize }

    return H_plot_settings , H_kwargs_errorbar , H_kwargs_errorbar_superimposed


