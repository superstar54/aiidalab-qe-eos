"""Results view widgets (MOVE TO OTHER MODULE!)

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""

from aiidalab_qe.panel import ResultPanel
import ipywidgets as ipw
import nglview
import traitlets


class EOSResults(ResultPanel):
    def __init__(self, wc_node, **kwargs):
        self.wc_node = wc_node

        self.summary_view = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>EOS.</h4></div>"""
        )
        #
        children = [self.summary_view, g]
        super().__init__(
            children=children,
            **kwargs,
        )

    def udpate(self):
        import plotly.graph_objects as go

        eos = self.wc_node.outputs.eos.get_dict()
        # init figure
        g = go.FigureWidget(
            layout=go.Layout(
                title=dict(text="XPS"),
                barmode="overlay",
            )
        )
        g.layout.xaxis.title = "Volume (A^3)"
        g.layout.xaxis.autorange = "Energy (eV)"
        for site, d in data.items():
            g.add_scatter(x=eos["vs"], y=eos["es"])
        self.children = [self.summary_view, g]
        
        