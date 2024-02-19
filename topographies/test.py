import plotly.graph_objects as go
import numpy as np

x = np.linspace(0, 1e3)
y = np.linspace(0, 1e3)

X, Y = np.meshgrid(x, y)

def h(x, y):
    return -2250+y*0.005

Z = h(X,Y)

fig = go.Figure(data=[go.Surface(x=X, y=Y,z=Z)])
fig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                  highlightcolor="limegreen", project_z=True))


fig.update_layout(title='Test', autosize=True,
                  #scene_camera_eye=dict(x=1.87, y=0.88, z=-0.64),
                  #width=500, height=500,
                  #margin=dict(l=65, r=50, b=65, t=90)
)


fig.show()