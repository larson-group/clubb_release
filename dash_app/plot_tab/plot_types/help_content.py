"""User-facing interpretation help for every Plot-tab plot family."""


PLOT_HELP = {
    "profile": {
        "title": "Profile plot",
        "overview": (
            "Shows the selected variable as a function of height. CLUBB profiles are "
            "averaged over the active global time window, unless point playback is active."
        ),
        "sections": (
            {
                "heading": "How to read it",
                "bullets": (
                    "The horizontal axis is the selected variable in the units stored in NetCDF; the vertical axis is height.",
                    "A line to the right or left means a larger or smaller value at that height. It does not imply a vertical flux unless the selected variable is itself a flux.",
                    "Enabled SAM and COAMPS profiles use their standard benchmark colors. In directory-compare mode, line style identifies the output directory.",
                ),
            },
            {
                "heading": "Time, columns, and error",
                "bullets": (
                    "Range mode averages all selected CLUBB records; point mode displays the record associated with the playback time.",
                    "Single-column mode uses the chosen CLUBB column. Overplot mode draws every enabled column rather than averaging them.",
                    "When a benchmark error is shown, it is a scaled RMSE over overlapping plotted heights. The scale is the benchmark profile range, with protective floors for cloud fraction and mean cloud water.",
                ),
            },
            {
                "heading": "Nuances",
                "bullets": (
                    "CLUBB and LES grids are interpolated only for the error calculation; their plotted lines retain their native vertical grids.",
                    "A good time-mean profile can hide timing or regime errors. Use time-series or time-height plots to check evolution.",
                    "Zooming changes what you see but does not change the selected averaging or error height range; use the global height control for that.",
                ),
            },
        ),
        "size": "medium",
    },
    "budget": {
        "title": "Budget plot",
        "overview": (
            "Shows the vertical profiles of the diagnosed terms in one CLUBB prognostic "
            "equation, all expressed as tendencies in the same units."
        ),
        "sections": (
            {
                "heading": "What the lines mean",
                "bullets": (
                    "The term ending in _bt is the diagnosed total time tendency of the budgeted variable.",
                    "The remaining named terms are individual physical or numerical contributions such as advection, buoyancy, pressure, turbulent production, dissipation, forcing, microphysics, or clipping. Their exact availability depends on the equation and stats file.",
                    "The residual is the diagnosed closure mismatch after the reported contributions are combined. A small residual relative to the dominant terms indicates that the written budget closes.",
                    "Positive terms increase the variable locally; negative terms decrease it. Large opposing terms may nearly cancel while still being dynamically important.",
                ),
            },
            {
                "heading": "Averaging and interpretation",
                "bullets": (
                    "Every line is averaged over the active global time selection before plotting.",
                    "Budget terms are rates, not accumulated changes. Multiplying a steady tendency by a time interval is only an approximation when terms evolve during that interval.",
                    "A time-mean residual can look small because positive and negative instantaneous residuals cancel. Point playback is useful for checking that possibility.",
                    "Terms absent from the selected NetCDF file are skipped; an incomplete legend can therefore mean incomplete stats output rather than a zero tendency.",
                ),
            },
            {
                "heading": "Current limitations",
                "bullets": (
                    "Budget plots use one CLUBB column and are disabled for output-directory comparison and all-column overplot mode.",
                    "Some operations outside a particular advance routine can appear in clipping, surface, forcing, or residual terms rather than the process term you first expect.",
                    "Budget closure does not by itself establish physical correctness: an equation can close exactly while one of its parameterized terms is biased.",
                ),
            },
        ),
        "size": "large",
    },
    "pdf_contour": {
        "title": "PDF contour plot",
        "overview": (
            "Shows CLUBB's diagnosed two-dimensional PDF contours over a probability-and-transport "
            "background. For ARM with SAM enabled, the background is a real raw 3-D SAM plane; "
            "otherwise it is CLUBB's analytic PDF."
        ),
        "sections": (
            {
                "heading": "Probability and signed transport",
                "bullets": (
                    "Gold marks high probability. The card-level SAM color switch selects signed local w′r꜀′, w′rₜ′, w′χ′, w′θₗ′, or the rtp2-style third-order kernel w′rₜ′². Blue marks positive contribution, red negative, and pale colors mean probability and the selected signal are both strong.",
                    "For the cloud-water default, positive covariance does not simply mean upward motion: both an unusually cloudy updraft and an unusually dry downdraft contribute positively. Cloudy downdrafts and dry updrafts contribute negatively. The rtp2 option is a third moment, not a transport flux.",
                    "With Probability only selected, purple boundaries enclose the highest-density 50%, 80%, and 99.5% of the CLUBB PDF. With a weighted signal selected, they instead enclose the corresponding fractions of CLUBB's absolute local contribution, matching the raw-SAM contour convention. Their line strength fades outward, so the 50% core is emphasized and the 99.5% envelope remains visible but quieter.",
                    "For a selected time range, the displayed density is the average of the PDFs at the included output times, not a PDF constructed from simply averaged parameters.",
                    "Probability and signed-transport magnitude use separate robust log normalizations before they are combined into color. Raw color intensity therefore supports geometric interpretation, not direct comparison of peak density values between projections.",
                    "All five signal choices are exact when the background is a raw SAM plane. For an analytic CLUBB-only background, only the cloud-water default is colored; the other choices deliberately display probability without an approximate red/blue field, while their purple contours still use the analytic conditional local contribution.",
                ),
            },
            {
                "heading": "CLUBB Gaussian components",
                "bullets": (
                    "Red, cyan, and—when present—green crosses mark CLUBB's Gaussian mixture components. The plot detects a nonzero component-3 weight and applies CLUBB's nested 3G weighting: components 1 and 2 split the remaining mass after component 3.",
                    "Dashed 2× RMS shapes show each component's width and covariance tilt after moment-matching the selected time window. Two-Gaussian output remains a two-component plot.",
                    "For an exactly bivariate Gaussian, a 2× RMS (Mahalanobis-radius 2) ellipse encloses about 86.5% of that component. It is a component-shape guide, not one of the total-PDF probability boundaries.",
                    "A component can be difficult to see when its weight or variance is very small, or when it lies inside the dense part of the other component. Hovering the shape reports its mean weight, centers, standard deviations, and correlation.",
                ),
            },
            {
                "heading": "Raw LES background",
                "bullets": (
                    "The former Gaussian LES moment proxy has been removed because it could look like sampled LES structure that the available profile moments do not contain.",
                    "When SAM is enabled and a matching <case>_3d directory exists, the app finds the nearest time and height in the shared JULY_2017_3D_RECREATIONS raw run. It bins the actual horizontal plane with the shared signed-transport convention and repeats that same background behind every CLUBB output panel.",
                    "With a raw plane present, CLUBB's analytic heatmap is disabled so the sampled LES structure remains visible; purple CLUBB lines and component shapes still come from CLUBB. Both purple and gold contours follow probability in Probability only mode, or absolute local contribution for the selected weighted signal. Gold raw-SAM lines use a lightly smoothed rank field, but their 50%, 80%, and 99.5% cutoffs are calibrated to the original weights. If no matching raw run exists, the background falls back to CLUBB's analytic probability and signed transport.",
                ),
            },
            {
                "heading": "Projection-specific nuances",
                "bullets": (
                    "w–rt probability uses total-water moments directly. Cloud transport is obtained by analytically conditioning each trivariate Gaussian on (w,rt) and evaluating its expected positive χ/cloud water.",
                    "w–χ uses CLUBB's linearized thermodynamic transform from rt and thl moments. The vertical χ=0 line is the saturation/cloud threshold.",
                    "Cloud water is rc=max(χ,0). In w–rc, the heatmap shows only the continuous cloudy density for rc>0; clear-air probability is concentrated at rc=0 and is not a two-dimensional contour. Axes and moment ellipses are clipped at that physical boundary.",
                    "For that reason, the signed w–rc integral omits covariance carried by the clear-air boundary mass at rc=0. Use signed w–χ or w–rt when you want the full mapped covariance accounting.",
                    "The annotation decomposes diagnosed w′rc′ into cloudy-mean transport and within-cloud covariance when all required fields are available.",
                    "Multiple CLUBB output directories are displayed in separate panels with shared axes. When raw SAM is available, each panel receives the identical sampled background so only the closure contours change.",
                ),
            },
        ),
        "size": "large",
    },
    "custom": {
        "title": "Custom profile expression",
        "overview": (
            "Evaluates an arithmetic expression from available profile fields and plots the "
            "result against height using the active time and column selection."
        ),
        "sections": (
            {
                "heading": "Expression syntax",
                "bullets": (
                    "Use field names directly, for example wpthlp + wprtp or clip(wprcp / les.wprcp, -2, 2).",
                    "Prefix a field with les. to use the configured LES benchmark, for example les.wprcp. Bare fields refer to CLUBB output.",
                    "Supported operations include ordinary arithmetic and a restricted set of safe functions such as abs, clip, sqrt, log, exp, erf, minimum, maximum, where, and trigonometric functions. np.function and numpy.function aliases are accepted for supported functions.",
                    "Keyword arguments, arbitrary Python, indexing, and unapproved attributes are intentionally rejected.",
                ),
            },
            {
                "heading": "Alignment and units",
                "bullets": (
                    "Fields are time-averaged first and then evaluated on compatible CLUBB profile/column selections. A single line may be broadcast across multiple selected lines when unambiguous.",
                    "LES fields are interpolated to the CLUBB height grid before participating in an expression.",
                    "The app cannot infer algebraic output units. Ratios involving small denominators should use clip or another explicit guard to avoid spikes and misleading axis ranges.",
                ),
            },
            {
                "heading": "Nuances",
                "bullets": (
                    "A ratio of time-averaged quantities is generally not equal to the time average of their instantaneous ratio.",
                    "In multi-directory mode the same expression is evaluated independently for each output. An LES-prefixed field remains the common benchmark.",
                    "Non-finite results are omitted from plotting; an empty plot usually means missing fields, incompatible shapes, or invalid arithmetic.",
                ),
            },
        ),
        "size": "large",
    },
    "timeseries": {
        "title": "Time-series plot",
        "overview": (
            "Shows how a scalar or column-resolved NetCDF variable evolves through the full "
            "time extent of the selected output file."
        ),
        "sections": (
            {
                "heading": "How it differs from other plots",
                "bullets": (
                    "The horizontal axis is file time and the vertical axis is the selected variable.",
                    "Unlike profile and budget plots, this plot intentionally uses the full file time extent and is not cropped by the global profile averaging window.",
                    "It is best suited to true time-series variables or fields whose non-time dimensions can be reduced to selected lines.",
                ),
            },
            {
                "heading": "Columns and comparison",
                "bullets": (
                    "Single-column mode shows the selected column. Overplot mode shows each enabled column separately.",
                    "With multiple output directories, line style identifies the directory. In directory-plus-all-column mode, color distinguishes columns while the directory legend describes line style.",
                    "Similar time series do not guarantee similar vertical structure; pair this view with profiles or time-height plots.",
                ),
            },
        ),
        "size": "medium",
    },
    "timeheight": {
        "title": "Time–height plot",
        "overview": (
            "Shows a variable as a two-dimensional field with time on the horizontal axis, "
            "height on the vertical axis, and color representing the variable value."
        ),
        "sections": (
            {
                "heading": "Color and time averaging",
                "bullets": (
                    "Read the colorbar for the variable's NetCDF units. Color represents the field value, not probability density unless the chosen field is explicitly a density.",
                    "The global start time sets the beginning of the displayed sequence. The selected average length forms consecutive time bins; each column of pixels is the average in one bin and is positioned at that bin's end time.",
                    "Small-scale extrema can weaken or disappear as the averaging interval grows.",
                ),
            },
            {
                "heading": "Columns and benchmarks",
                "bullets": (
                    "Single-column mode uses the selected CLUBB column. All-column mode averages enabled CLUBB columns before plotting rather than drawing overlapping heatmaps.",
                    "When compatible LES panels are enabled, CLUBB and benchmark panels use common color limits so colors have the same meaning. The card is locked to large size to keep stacked panels readable.",
                    "Time–height plots are currently disabled for multi-directory comparison; use separate cards or runs for that comparison.",
                ),
            },
        ),
        "size": "medium",
    },
    "subcolumn": {
        "title": "Subcolumn plot",
        "overview": (
            "Shows the vertical profiles of CLUBB subcolumn samples for a selected variable, "
            "revealing variability hidden by a grid-mean profile."
        ),
        "sections": (
            {
                "heading": "What each line represents",
                "bullets": (
                    "When an explicit subcolumn dimension exists, each line is one sampled subcolumn profile after the active time averaging.",
                    "If no explicit subcolumn dimension exists, the plot falls back to available column-resolved profiles. The note above the graph states which interpretation is active.",
                    "The spread and clustering of lines describe sampled heterogeneity; they are not uncertainty bars and should not be interpreted as independent observations.",
                ),
            },
            {
                "heading": "Time and columns",
                "bullets": (
                    "Range mode averages each subcolumn line over the selected time records; point mode follows playback.",
                    "Single-column mode uses the chosen CLUBB column. All-column mode includes lines from every enabled column and can produce a dense plot.",
                    "Subcolumn plots are disabled in multi-directory compare mode and require compatible subcolumn or column dimensions in the output.",
                ),
            },
        ),
        "size": "medium",
    },
}
