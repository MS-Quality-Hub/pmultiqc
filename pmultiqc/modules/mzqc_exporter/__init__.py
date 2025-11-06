from .mzqc_exporter import MzQCExporterModule

from typing import Any, Dict, Optional, Sequence, Union

from mzqc import MZQCFile as qc

from multiqc.plots import bargraph, linegraph, box, scatter, table
from multiqc.plots.bargraph import InputDatasetT, InputCategoriesT, BarPlot, BarPlotConfig


# initialize the module
__all__ = ["MzQCExporterModule"]


# some static / special functions to wrap the plots and add metrices

def plot_bargraph_and_add_mzqc( 
        data: Union[InputDatasetT, Sequence[InputDatasetT]],
        cats: Optional[Union[InputCategoriesT, Sequence[InputCategoriesT]]] = None,
        pconfig: Optional[Union[Dict[str, Any], BarPlotConfig]] = None,
        mzqcexporter: MzQCExporterModule = None,
        accession: str = None,
)  -> Union["BarPlot", str, None]:
    # just forward the execution of the plot
    bar_html = bargraph.plot(
        data=data,
        cats=cats,
        pconfig=pconfig,
    )

    # add the metric for this exporter
    if mzqcexporter is not None:
        mzqcexporter.add_metric(data, accession)

    return bar_html
