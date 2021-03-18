import pandas as pd
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource

df = pd.read_csv('kde_data.csv')

source = ColumnDataSource(data=df)

p = figure(title = "Your Complex Spin Splitting Energy Compared to Training Data", 
               sizing_mode="scale_both", plot_width=400, plot_height=250)
p.xaxis.axis_label = 'Spin Splitting Energy (kcal/mol)'
p.yaxis.axis_label = 'KDE of Training Data'
p.line('sse', 'kde', source=source, line_width=2)

show(p)
