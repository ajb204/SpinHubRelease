WORKFLOW_REGISTRY = [
    {
        "key": "prepare",
        "title": "Prepare raw data",
        "description": "Select the working directory, choose the nmrPipe input file, and generate processed spectra.",
    },
    {
        "key": "inspect",
        "title": "Inspect projections",
        "description": "Open 2D projections to check data quality, contours, and spectral regions.",
    },
    {
        "key": "decon",
        "title": "Detect peaks / deconvolve",
        "description": "Run deconvolution and peak detection, then review the fitted outputs.",
    },
    {
        "key": "slices",
        "title": "Explore slices",
        "description": "Open 1D and 2D slice viewers for detailed local inspection.",
    },
    {
        "key": "special",
        "title": "Specialized analyses",
        "description": "Launch pseudo-2D, pseudo-3D, phasing, uSTA, or full 3D analysis as needed.",
    },
]

WORKFLOW_BY_KEY = {item["key"]: item for item in WORKFLOW_REGISTRY}
