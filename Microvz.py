#This program generates visualizations of microbiome data in the form of simplified conceptual models.
#The program is designed to be run from the command line, with the following arguments:
#Use argparse to parse command line arguments
import argparse
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("--parse_only", help="Flag to skip running any programs and only parse their output using selected parameters", default=False, type=bool)
args = parser.parse_args()
import plotly.graph_objects as go
import numpy as np

cube_fig = go.Figure(data=[
    go.Mesh3d(
        golden_ratio = (1 + 5 ** 0.5) / 2
        print(golde)
        a = [[0,1,-1,golden_ratio,-golden_ratio],[1,-1,golden_ratio,-golden_ratio,0],[golden_ratio,-golden_ratio,0,1,-1]]
        #12 vertices of an icosahedron
        ico_vertices = list(itertools.product(*a))
        print(ico_vertices)
        x = [0]
        y = [1]
        z = [golden_ratio]
        colorbar_title='z',
        colorscale=[[0, 'green'],
                    [1, 'blue']],
        # Intensity of each vertex, which will be interpolated and color-coded
        intensity = [1,1,1,1,0,0,0,0],
        # i, j and k give the vertices of triangles
        alphahull=0,
        opacity=0.8,
        name='y',
        showscale=False
    )
])

#save figure
cube_fig.write_html("cube.html")